/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.BarcodeEditDistanceQuery;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.text.NumberFormat;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * Determine the barcode for each read in an Illumina lane.
 * For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt.
 * An output file contains a line for each read in the tile, aligned with the regular basecall output
 * The output file contains the following tab-separated columns:
 * - read subsequence at barcode position
 * - Y or N indicating if there was a barcode match
 * - matched barcode sequence (empty if read did not match one of the barcodes).  If there is no match
 * but we're close to the threshold of calling it a match we output the barcode that would have been
 * matched but in lower case
 * - distance to best matching barcode, "mismatches" (*)
 * - distance to second-best matching barcode, "mismatchesToSecondBest" (*)
 *
 * NOTE (*): Due to an optimization the reported mismatches & mismatchesToSecondBest values may be inaccurate as long as
 * the conclusion (match vs. no-match) isn't affected. For example, reported mismatches and
 * mismatchesToSecondBest may be smaller than their true value if mismatches is truly larger than MAX_MISMATCHES.
 * Also, mismatchesToSecondBest might be smaller than its true value if its true value is greater than
 * mismatches + MIN_MISMATCH_DELTA.
 */
@CommandLineProgramProperties(

        summary = ExtractIlluminaBarcodes.USAGE_SUMMARY + ExtractIlluminaBarcodes.USAGE_DETAILS,
        oneLineSummary = ExtractIlluminaBarcodes.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class ExtractIlluminaBarcodes extends ExtractBarcodesProgram {
    static final String USAGE_SUMMARY = "Tool determines the barcode for each read in an Illumina lane.  ";
    static final String USAGE_DETAILS = "<p>This tool determines the numbers of reads containing barcode-matching sequences and provides " +
            "statistics on the quality of these barcode matches.</p> " +
            "<p>Illumina sequences can contain at least two types of barcodes, sample and molecular (index).  Sample barcodes " +
            "(B in the read structure) are used to demultiplex pooled samples while index barcodes (M in the read structure) are used " +
            "to differentiate multiple reads of a template when carrying out paired-end sequencing.  Note that this tool only extracts " +
            "sample (B) and not molecular barcodes (M).</p>" +
            "" +
            "<p>Barcodes can be provided in the form of a list (BARCODE_FILE) or a string representing the barcode (BARCODE).  " +
            "The BARCODE_FILE contains multiple fields including '" + SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "' (or '" + SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "_1'), '"
            + SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "_2' (optional), '" + SampleInputParameters.BARCODE_NAME_COLUMN + "', and '" + SampleInputParameters.LIBRARY_NAME_COLUMN + "'. " +
            "In contrast, the BARCODE argument is used for runs with reads containing a single " +
            "barcode (nonmultiplexed) and can be added directly as a string of text e.g. BARCODE=CAATAGCG.</p>" +
            "" +
            "<p>Data is output per lane/tile within the BaseCalls directory with the file name format of 's_{lane}_{tile}_barcode.txt'.  " +
            "These files contain the following tab-separated columns:" +
            "<ul> " +
            "<li>Read subsequence at barcode position</li>" +
            "<li>Y or N indicating if there was a barcode match</li>" +
            "<li>Matched barcode sequence (empty if read did not match one of the barcodes)</li>  " +
            "<li>The number of mismatches if there was a barcode match</li>  " +
            "<li>The number of mismatches to the second best barcode if there was a barcode match</li>  " +
            "</ul>" +
            "If there is no match but we're close to the threshold of calling it a match, we output the barcode that would have been " +
            "matched but in lower case.  Threshold values can be adjusted to accommodate barcode sequence mismatches from the reads." +
            "  The metrics file produced by the ExtractIlluminaBarcodes program indicates the number of matches (and mismatches)" +
            " between the barcode reads and the actual barcodes.  These metrics are provided both per-barcode and per lane and can be " +
            "found in the BaseCalls directory.</p>" +
            "<p>For poorly matching barcodes, the order of specification of barcodes can cause arbitrary output differences.</p>" +
            "" +
            "<h4>Usage example:</h4> " +
            "<pre>" +
            "java -jar picard.jar ExtractIlluminaBarcodes \\<br />" +
            "              BASECALLS_DIR=/BaseCalls/ \\<br />" +
            "              LANE=1 \\<br />" +
            "          READ_STRUCTURE=25T8B25T \\<br />" +
            "              BARCODE_FILE=barcodes.txt \\<br />" +
            "              METRICS_FILE=metrics_output.txt " +
            "</pre>" +
            "" +
            "Please see the ExtractIlluminaBarcodes.BarcodeMetric " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric'>definitions</a> " +
            "for a complete description of the metrics produced by this tool.</p>" +
            "" +
            "<hr />";

    // The following attributes define the command-line arguments

    @Argument(doc = "Tab-delimited file of barcode sequences, barcode name and, optionally, library name.  " +
            "Barcodes must be unique and all the same length.  Column headers must be '" +  SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "' (or '" + SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "_1'), '"
            + SampleInputParameters.BARCODE_SEQUENCE_COLUMN + "_2' (optional), '"+ SampleInputParameters.BARCODE_NAME_COLUMN + "', and '" + SampleInputParameters.LIBRARY_NAME_COLUMN + "'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Argument(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Argument(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    private static final Log LOG = Log.getInstance(ExtractIlluminaBarcodes.class);

    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();

    public ExtractIlluminaBarcodes() {
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(METRICS_FILE);
        if (OUTPUT_DIR == null) {
            OUTPUT_DIR = BASECALLS_DIR;
        }
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIR);

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = Runtime.getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            numProcessors = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            numProcessors = NUM_PROCESSORS;
        }

        LOG.info("Processing with " + numProcessors + " PerTileBarcodeExtractor(s).");
        final ThreadPoolExecutor pool = new ThreadPoolExecutorWithExceptions(numProcessors);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }

        final BarcodeMetric noMatchMetric =  new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        final BarcodeExtractor barcodeExtractor = createBarcodeExtractor();

        final List<PerTileBarcodeExtractor> extractors = new ArrayList<>(factory.getAvailableTiles().size());
        // TODO: This is terribly inefficient; we're opening a huge number of files via the extractor constructor and we never close them.
        for (final int tile : factory.getAvailableTiles()) {
            final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                    tile,
                    getBarcodeFile(tile),
                    factory,
                    barcodeExtractor);

            extractors.add(extractor);
        }

        for (final PerTileBarcodeExtractor extractor : extractors) {
            pool.submit(extractor);
        }
        pool.shutdown();
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Per tile extractor executor", pool, Duration.ofMinutes(5));

        LOG.info("Processed " + extractors.size() + " tiles.");
        LOG.info("Cache grew to " + this.barcodeLookupMap.size() + " entries.");
        for (final PerTileBarcodeExtractor extractor : extractors) {
            for (final String key : barcodeToMetrics.keySet()) {
                barcodeToMetrics.get(key).merge(extractor.getMetrics().get(key));
            }
            noMatchMetric.merge(extractor.getNoMatchMetric());
            if (extractor.getException() != null) {
                LOG.error("Abandoning metrics calculation because one or more PerTileBarcodeExtractors failed.");
                return 4;
            }
        }

        // Finish metrics tallying.
        finalizeMetrics(barcodeToMetrics, noMatchMetric);

        // Warn about minimum qualities and assert that we've achieved the minimum.
        for (final Map.Entry<Byte, Integer> entry : bclQualityEvaluationStrategy.getPoorQualityFrequencies().entrySet()) {
            LOG.warn(String.format("Observed low quality of %s %s times.", entry.getKey(), entry.getValue()));
        }
        bclQualityEvaluationStrategy.assertMinimumQualities();

        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        INPUT_FILE = BARCODE_FILE;
        readStructure = new ReadStructure(READ_STRUCTURE.replaceAll("[TM]", "S"));
        return super.customCommandLineValidation();
    }

    public static void finalizeMetrics(final Map<String, BarcodeMetric> barcodeToMetrics,
                                       final BarcodeMetric noMatchMetric) {
        // Finish metrics tallying.
        long totalReads = noMatchMetric.READS;
        long totalPfReads = noMatchMetric.PF_READS;
        long totalPfReadsAssigned = 0;
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            totalReads += barcodeMetric.READS;
            totalPfReads += barcodeMetric.PF_READS;
            totalPfReadsAssigned += barcodeMetric.PF_READS;
        }

        if (totalReads > 0) {
            noMatchMetric.PCT_MATCHES = noMatchMetric.READS / (double) totalReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PCT_MATCHES = barcodeMetric.READS / (double) totalReads;
                if (barcodeMetric.PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        if (totalPfReads > 0) {
            noMatchMetric.PF_PCT_MATCHES = noMatchMetric.PF_READS / (double) totalPfReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PF_PCT_MATCHES = barcodeMetric.PF_READS / (double) totalPfReads;
                if (barcodeMetric.PF_PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PF_PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        // Calculate the normalized matches
        if (totalPfReadsAssigned > 0) {
            final double mean = (double) totalPfReadsAssigned / (double) barcodeToMetrics.values().size();
            for (final BarcodeMetric m : barcodeToMetrics.values()) {
                m.PF_NORMALIZED_MATCHES = m.PF_READS / mean;
            }
        }
    }

    /**
     * Create a barcode filename corresponding to the given tile qseq file.
     */
    private File getBarcodeFile(final int tile) {
        return new File(OUTPUT_DIR,
                "s_" + LANE + "_" + tileNumberFormatter.format(tile) + "_barcode.txt" + (COMPRESS_OUTPUTS ? ".gz" : ""));
    }

    /**
     * Extracts barcodes and accumulates metrics for an entire tile.
     */
    public static class PerTileBarcodeExtractor implements Runnable {
        private final int tile;
        private final File barcodeFile;
        private final Map<String, BarcodeMetric> metrics;
        private final BarcodeMetric noMatch;
        private Exception exception = null;
        private final boolean usingQualityScores;
        private BaseIlluminaDataProvider provider;
        private final ReadStructure outputReadStructure;
        private final BarcodeExtractor barcodeExtractor;

        /**
         * Constructor
         *  @param tile               The number of the tile being processed; used for logging only.
         *  @param barcodeFile        The file to write the barcodes to
         */
        public PerTileBarcodeExtractor(
                final int tile,
                final File barcodeFile,
                final IlluminaDataProviderFactory factory,
                final BarcodeExtractor extractor
        ) {
            this.barcodeExtractor = extractor;
            this.tile = tile;
            this.barcodeFile = barcodeFile;
            this.usingQualityScores = barcodeExtractor.minimumBaseQuality > 0;
            this.metrics = new LinkedHashMap<>(barcodeExtractor.metrics.size());
            for (final String key : barcodeExtractor.metrics.keySet()) {
                this.metrics.put(key, BarcodeMetric.copy(barcodeExtractor.metrics.get(key)));
            }

            this.noMatch = BarcodeMetric.copy(barcodeExtractor.noMatch);
            this.provider = factory.makeDataProvider(tile);
            this.outputReadStructure = factory.getOutputReadStructure();
        }

        // These methods return the results of the extraction
        public synchronized Map<String, BarcodeMetric> getMetrics() {
            return this.metrics;
        }

        public synchronized BarcodeMetric getNoMatchMetric() {
            return this.noMatch;
        }

        public synchronized Exception getException() {
            return this.exception;
        }

        /**
         * run method which extracts barcodes and accumulates metrics for an entire tile
         */
        public synchronized void run() {
            try {
                LOG.info("Extracting barcodes for tile " + tile);

                // Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                // is non-overlapping sets of files so make the  data providers in the individual threads for PerTileBarcodeExtractors
                // so they are not all waiting for each others file operations

                // Most likely we have SKIPS in our read structure since we replace all template reads with skips in the input data structure
                // (see customCommnandLineValidation), therefore we must use the outputReadStructure to index into the output cluster data
                final int[] barcodeIndices = outputReadStructure.sampleBarcodes.getIndices();
                final BufferedWriter writer = IOUtil.openFileForBufferedWriting(barcodeFile);
                final byte[][] barcodeSubsequences = new byte[barcodeIndices.length][];
                final byte[][] qualityScores = usingQualityScores ? new byte[barcodeIndices.length][] : null;
                while (provider.hasNext()) {
                    // Extract the barcode from the cluster and write it to the file for the tile
                    final ClusterData cluster = provider.next();

                    for (int i = 0; i < barcodeIndices.length; i++) {
                        barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
                        if (usingQualityScores) {
                            qualityScores[i] = cluster.getRead(barcodeIndices[i]).getQualities();
                        }
                    }
                    final boolean passingFilter = cluster.isPf();
                    final BarcodeExtractor.BarcodeMatch match = barcodeExtractor.findBestBarcode(barcodeSubsequences, qualityScores);

                    BarcodeExtractor.updateMetrics(match, passingFilter, metrics, noMatch);

                    final String yOrN = (match.matched ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(StringUtil.bytesToString(bc));
                    }
                    writer.write("\t" + yOrN + "\t" + match.barcode + "\t" + match.mismatches + "\t" + match.mismatchesToSecondBest);
                    writer.newLine();
                }
                writer.close();
            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                CloserUtil.close(provider);
                provider = null;
            }
        }

    }
}
