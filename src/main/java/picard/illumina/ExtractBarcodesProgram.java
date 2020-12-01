package picard.illumina;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.illumina.parser.*;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;

import java.io.File;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public abstract class ExtractBarcodesProgram extends CommandLineProgram {
    @Argument(doc = "The distance metric that should be used to compare the barcode-reads and the provided barcodes for finding the best and second-best assignments.")
    public DistanceMetric DISTANCE_MODE = DistanceMetric.HAMMING;

    @Argument(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Argument(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Argument(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Argument(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even if the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(shortName = "GZIP", doc = "Compress output FASTQ files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    @Argument(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    @Argument(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true)
    public File METRICS_FILE;

    @Argument
    public Boolean BARCODE_EXTRACT = false;

    @Argument(optional = true)
    public File INPUT_PARAMS_FILE;

    protected IlluminaDataProviderFactory factory;
    protected final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
    protected final static int MAX_LOOKUP_SIZE = 4096;
    protected final ConcurrentHashMap<BarcodeExtractor.ByteString, BarcodeExtractor.BarcodeMatch> barcodeLookupMap = new ConcurrentHashMap<>(MAX_LOOKUP_SIZE);
    protected final Set<BarcodeExtractor.ByteString> barcodeByteStrings = new HashSet<>();
    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);
    protected BarcodeMetric noMatchMetric;
    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();

    /**
     * The read structure of the actual Illumina Run, i.e. the readStructure of the input data
     */
    protected ReadStructure readStructure;

    protected BarcodeExtractor createBarcodeExtractor() {
        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.sampleBarcodes.length()];
        final byte[][] perfectScores = new byte[readStructure.sampleBarcodes.length()][];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                perfectScores[index] = new byte[d.length];
                Arrays.fill(perfectScores[index], (byte) 60);
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }

        this.noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        for (final BarcodeMetric metric : barcodeToMetrics.values()) {
            this.barcodeByteStrings.add(new BarcodeExtractor.ByteString(metric.barcodeBytes));
        }

        final BarcodeExtractor barcodeExtractor = new BarcodeExtractor(barcodeToMetrics,
                barcodeLookupMap,
                barcodeByteStrings,
                noMatchMetric,
                MAX_NO_CALLS,
                MAX_MISMATCHES,
                MIN_MISMATCH_DELTA,
                MINIMUM_BASE_QUALITY,
                DISTANCE_MODE);

        // Prepopulate the lookup map with all perfect barcodes
        for(BarcodeMetric metric : barcodeToMetrics.values()){
            BarcodeExtractor.BarcodeMatch match = barcodeExtractor.calculateBarcodeMatch(metric.barcodeBytes,
                    perfectScores);
            barcodeLookupMap.put(new BarcodeExtractor.ByteString(metric.barcodeBytes), match);
        }

        // Prepopulate all no call barcode match
        BarcodeExtractor.BarcodeMatch noCallMatch = barcodeExtractor.calculateBarcodeMatch(noMatchMetric.barcodeBytes,
                perfectScores);

        barcodeLookupMap.put(new BarcodeExtractor.ByteString(noMatchMetric.barcodeBytes), noCallMatch);
        return barcodeExtractor;
    }

    /**
     * Validate that POSITION >= 1, and that all BARCODEs are the same length and unique
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        readStructure = new ReadStructure(READ_STRUCTURE);
        final List<String> messages = new ArrayList<>();
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);

        final Set<IlluminaDataType> datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores)) :
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF));
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy, datatypes);

        if (INPUT_PARAMS_FILE != null) {
            SampleInputParameters.parseInputFile(INPUT_PARAMS_FILE, readStructure, barcodeToMetrics, messages);
            if (barcodeToMetrics.keySet().isEmpty()) {
                messages.add("No barcodes have been specified.");
            }
        }

       return messages.toArray(new String[0]);
    }

    protected String[] collectErrorMessages(List<String> messages, String[] superErrors) {
        if(superErrors != null && superErrors.length > 0) {
            messages.addAll(Arrays.asList(superErrors));
        }

        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[0]);
    }

    protected void outputMetrics() {
        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
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
}
