package picard.illumina;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.illumina.parser.*;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public abstract class ExtractBarcodesProgram extends CommandLineProgram {

    @Argument(doc = "Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"INPUT_FILE"})
    public List<String> BARCODE = new ArrayList<>();

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

    @Argument(optional = true)
    public File INPUT_FILE;

    protected IlluminaDataProviderFactory factory;
    protected final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();

    protected final static int MAX_LOOKUP_SIZE = 4096;
    protected final ConcurrentHashMap<BarcodeExtractor.ByteString, BarcodeExtractor.BarcodeMatch> barcodeLookupMap = new ConcurrentHashMap<>(MAX_LOOKUP_SIZE);
    protected final Set<BarcodeExtractor.ByteString> barcodeByteStrings = new HashSet<>();
    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);
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

        BarcodeMetric noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

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
        final ArrayList<String> messages = new ArrayList<>();
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);
        /*
          In extract illumina barcodes we NEVER want to look at the template reads nor the molecular barcodes, therefore replace them with
          skips because IlluminaDataProvider and its factory will neither open these nor produce ClusterData with the template reads in them,
          thus reducing the file IO and value copying done by the data provider
         */

        final Set<IlluminaDataType> datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores)) :
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF));
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy, datatypes);

        if (INPUT_FILE != null) {
            SampleInputParameters.parseInputFile(INPUT_FILE, readStructure, barcodeToMetrics, messages);
        } else {
            final int numBarcodes = readStructure.sampleBarcodes.length();
            final Set<String> barcodes = new HashSet<>();

            for (final String barcode : BARCODE) {
                if (barcodes.contains(barcode)) {
                    messages.add("Barcode " + barcode + " specified more than once.");
                }
                barcodes.add(barcode);
                int barcodeNum = 0;
                int pos = 0;
                final String[] bcStrings = new String[numBarcodes];
                for (final ReadDescriptor rd : readStructure.descriptors) {
                    if (rd.type != ReadType.Barcode) {
                        continue;
                    }
                    bcStrings[barcodeNum] = barcode.substring(pos, pos + rd.length);
                    pos += rd.length;
                    ++barcodeNum;
                }

                final BarcodeMetric metric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(bcStrings), bcStrings);
                barcodeToMetrics.put(barcode, metric);
            }
        }
        if ((INPUT_FILE != null || !BARCODE.isEmpty()) && barcodeToMetrics.keySet().isEmpty()) {
            messages.add("No barcodes have been specified.");
        }
        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }
}
