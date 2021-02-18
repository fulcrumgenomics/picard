package picard.illumina;

import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.Tuple;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * This class encloses sample input parameters used for any ExtractBarcodesProgram (ExtractIlluminaBarcodes,
 * IlluminaBasecallsToFastq, IlluminaBasecallsToSam). This is used to parse tab delimited sample information including
 * barcodes, library names and output prefixes.
 */
public class SampleInputParameters {
    /**
     * Column header for the first barcode sequence (preferred).
     */

    public static final Set<String> BARCODE_PREFIXES = new HashSet<>(
            Arrays.asList(SampleInputParameters.BARCODE_SEQUENCE_COLUMN, SampleInputParameters.BARCODE_COLUMN)
    );

    public static final String BARCODE_COLUMN = "barcode";
    public static final String BARCODE_SEQUENCE_COLUMN = "barcode_sequence";

    /**
     * Column header for the barcode name.
     */
    public static final String BARCODE_NAME_COLUMN = "barcode_name";
    /**
     * Column header for the library name.
     */
    public static final String LIBRARY_NAME_COLUMN = "library_name";


    /**
     * Parses any one of the following types of files:
     *
     * ExtractIlluminaBarcodes      BARCODE_FILE
     * IlluminaBasecallsToFastq     MULTIPLEX_PARAMS
     * IlluminaBasecallsToSam       LIBRARY_PARAMS
     *
     * This will validate to file format as well as populate a Map of barcodes to metrics.
     *
     * @param inputFile         The input file that is being parsed
     * @param readStructure     The read structure for the reads of the run
     */
    protected static Tuple<Map<String, BarcodeMetric>, List<String>> parseInputFile(final File inputFile,
                                          final ReadStructure readStructure) {
        List<String> messages = new ArrayList<>();
        Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
        try (final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(inputFile)) {
            List<String> validBarcodeColumns = barcodesParser.columnLabels().stream().filter(name -> {
                boolean isValidPrefix = false;
                for (String columnPrefix : BARCODE_PREFIXES) {
                    isValidPrefix |= name.toUpperCase().startsWith(columnPrefix.toUpperCase());
                }
                return isValidPrefix;
            }).collect(Collectors.toList());

            if (readStructure.sampleBarcodes.length() != validBarcodeColumns.size()) {
                messages.add("Expected " + readStructure.sampleBarcodes.length() + " valid barcode columns, but only found " +
                        String.join(",", validBarcodeColumns));
            }

            Matcher matcher = Pattern.compile("^(.*)_\\d").matcher(validBarcodeColumns.get(0));

            final String sequenceColumn;
            boolean hasMultipleNumberedBarcodeColumns = matcher.matches();
            if (hasMultipleNumberedBarcodeColumns) {
                sequenceColumn = matcher.group(1);
            } else {
                sequenceColumn = validBarcodeColumns.get(0);
            }
            final boolean hasBarcodeName = barcodesParser.hasColumn(BARCODE_NAME_COLUMN);
            final boolean hasLibraryName = barcodesParser.hasColumn(LIBRARY_NAME_COLUMN);
            final int numBarcodes = readStructure.sampleBarcodes.length();
            final Set<String> barcodes = new HashSet<>();
            for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
                final String[] bcStrings = new String[numBarcodes];
                int barcodeNum = 0;
                for (final ReadDescriptor rd : readStructure.descriptors) {
                    if (rd.type != ReadType.Barcode) {
                        continue;
                    }
                    final String header = hasMultipleNumberedBarcodeColumns ? sequenceColumn + "_" + (1 + barcodeNum) : sequenceColumn;
                    final String field = row.getField(header);
                    if (field == null) {
                        messages.add(String.format("Null barcode in column %s of row: %s", header, row.getCurrentLine()));
                        bcStrings[barcodeNum] = "";
                    } else {
                        bcStrings[barcodeNum] = field;
                    }
                    ++barcodeNum;
                }
                final String bcStr = IlluminaUtil.barcodeSeqsToString(bcStrings);
                // if the barcode is all Ns don't add it to metrics (we add noCallMetric separately)
                if(bcStr.contains("N") || bcStr.contains("n")) continue;
                if (barcodes.contains(bcStr)) {
                    messages.add("Barcode " + bcStr + " specified more than once in " + inputFile);
                }
                barcodes.add(bcStr);
                final String barcodeName = (hasBarcodeName ? row.getField(BARCODE_NAME_COLUMN) : "");
                final String libraryName = (hasLibraryName ? row.getField(LIBRARY_NAME_COLUMN) : "");
                final BarcodeMetric metric = new BarcodeMetric(barcodeName, libraryName, bcStr, bcStrings);
                barcodeToMetrics.put(StringUtil.join("", bcStrings), metric);
            }
        }
        return new Tuple<>(barcodeToMetrics, messages);
    }
}
