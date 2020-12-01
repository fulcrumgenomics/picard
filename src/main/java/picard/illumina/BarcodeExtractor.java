package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.SequenceUtil;
import picard.util.BarcodeEditDistanceQuery;
import picard.util.IlluminaUtil;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

public class BarcodeExtractor {
    private final Map<String, BarcodeMetric> metrics;
    private final BarcodeMetric noMatch;
    final Set<ByteString> barcodesBytes;
    final int maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality;
    final DistanceMetric distanceMode;
    final ConcurrentHashMap<ByteString, BarcodeMatch> barcodeLookupMap;

    public BarcodeExtractor(final Map<String, BarcodeMetric> barcodeToMetrics,
                            final ConcurrentHashMap<ByteString, BarcodeMatch> barcodeLookupMap,
                            Set<ByteString> barcodeByteStrings, final BarcodeMetric noMatchMetric,
                            final int maxNoCalls,
                            final int maxMismatches,
                            final int minMismatchDelta,
                            final int minimumBaseQuality,
                            final DistanceMetric distanceMode) {
        this.metrics = barcodeToMetrics;
        this.barcodeLookupMap = barcodeLookupMap;
        this.barcodesBytes = barcodeByteStrings;
        this.noMatch = noMatchMetric;
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        this.distanceMode = distanceMode;
    }

    public Map<String, BarcodeMetric> getMetrics() {
        return this.metrics;
    }

    public BarcodeMetric getNoMatchMetric() {
        return this.noMatch;
    }
    /**
     * Find the best barcode match for the given read sequence, and accumulate metrics
     *
     * NOTE: the returned BarcodeMatch object will contain mismatches mismatchesToSecondBest values that may be
     * inaccurate as long as the conclusion match/no-match isn't affected. for example, mismatches and mismatchesToSecondBest
     * may be smaller than their true value if mismatches is truly larger than maxMismatches.
     * Also, mismatchesToSecondBest might be smaller than its true value if its true value is greater than
     * mismatches + minMismatchDelta. This is due to an optimization which allows the distance calculation to stop once
     * the conclusion (Match or no-Match) can be reached.
     *
     * @param readSubsequences portion of read containing barcode
     * @return perfect barcode string, if there was a match within tolerance, or null if not.
     */
    BarcodeMatch findBestBarcode(final byte[][] readSubsequences,
                                 final byte[][] qualityScores) {
        final boolean canUseLookupTable = areAllQualitiesAboveMinimum(qualityScores, minimumBaseQuality);
        final ByteString barcodesAsString = new ByteString(readSubsequences);
        BarcodeMatch match = canUseLookupTable ? barcodeLookupMap.get(barcodesAsString) : null;
        if (match == null) {
            match = calculateBarcodeMatch(readSubsequences, qualityScores);

            if (canUseLookupTable && match.isMatched()) {
                barcodeLookupMap.put(barcodesAsString, match);
            }
        }

        return match;
    }

    BarcodeMatch calculateBarcodeMatch(final byte[][] readSubsequences,
                                       final byte[][] qualityScores) {
        final BarcodeMatch match;
        String bestBarcode = null;
        match = new BarcodeMatch();

        int totalBarcodeReadBases = 0;
        int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

        for (final byte[] bc : readSubsequences) {
            totalBarcodeReadBases += bc.length;
            for (final byte b : bc) {
                if (SequenceUtil.isNoCall(b)) {
                    ++numNoCalls;
                }
            }
        }

        // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
        // base is a mismatch.
        int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
        int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;

        for (final ByteString barcodeBytes : barcodesBytes) {
            // need to add maxMismatches + minMismatchDelta together since the result might get used as numMismatchesInSecondBestBarcode
            final BarcodeEditDistanceQuery barcodeEditDistanceQuery = new BarcodeEditDistanceQuery(barcodeBytes.bytes, readSubsequences, qualityScores,
                    minimumBaseQuality, Math.min(maxMismatches, numMismatchesInBestBarcode) + minMismatchDelta);
            final int numMismatches = distanceMode.distance(barcodeEditDistanceQuery);

            if (numMismatches < numMismatchesInBestBarcode) {
                if (bestBarcode != null) {
                    numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
                }
                numMismatchesInBestBarcode = numMismatches;
                bestBarcode = barcodeBytes.toString();
            } else if (numMismatches < numMismatchesInSecondBestBarcode) {
                numMismatchesInSecondBestBarcode = numMismatches;
            }
        }

        match.matched = bestBarcode != null &&
                numNoCalls <= maxNoCalls &&
                numMismatchesInBestBarcode <= maxMismatches &&
                numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMismatchDelta;

        match.mismatches = numMismatchesInBestBarcode;
        match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;

        if (match.matched) {
            match.barcode = bestBarcode;
        } else {
            // If we have something that's not a "match" but matches one barcode
            // slightly, we output that matching barcode in lower case
            if (numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases && bestBarcode != null) {
                match.barcode = bestBarcode.toLowerCase();
            } else {
                match.mismatches = totalBarcodeReadBases;
                match.barcode = "";
            }
        }

        return match;
    }

    static void updateMetrics(final BarcodeMatch match, final boolean passingFilter,
                                      final Map<String, BarcodeMetric> metrics, final BarcodeMetric noMatchBarcodeMetric) {
        if (match.matched) {
            final BarcodeMetric matchMetric = metrics.get(match.barcode);
            ++matchMetric.READS;
            if (passingFilter) {
                ++matchMetric.PF_READS;
            }
            if (match.mismatches == 0) {
                ++matchMetric.PERFECT_MATCHES;
                if (passingFilter) {
                    ++matchMetric.PF_PERFECT_MATCHES;
                }
            } else if (match.mismatches == 1) {
                ++matchMetric.ONE_MISMATCH_MATCHES;
                if (passingFilter) {
                    ++matchMetric.PF_ONE_MISMATCH_MATCHES;
                }
            }
        } else {
            ++noMatchBarcodeMetric.READS;
            if (passingFilter) {
                ++noMatchBarcodeMetric.PF_READS;
            }
        }
    }

    /**
     * Checks to ensure that all quality values are greater that the given cutoff such that comparisons can
     * be done just using the bases without further reference to quality scores.
     */
    private static boolean areAllQualitiesAboveMinimum(final byte[][] qualityScores, final int minimumBaseQuality) {
        if (qualityScores == null) return true;

        for (final byte[] qs : qualityScores) {
            for (final byte q : qs) {
                if (q < minimumBaseQuality) {
                    return false;
                }
            }
        }

        return true;
    }
    /**
     * Utility class to hang onto data about the best match for a given barcode
     */
    public static class BarcodeMatch {
        boolean matched;
        String barcode;
        int mismatches;
        int mismatchesToSecondBest;

        public boolean isMatched() {
            return matched;
        }

        public String getBarcode() {
            return barcode;
        }
    }

    /**
     * Class to give a byte[][] a hashcode and equals without copying the whole contents into a String.
     */
    static final class ByteString {
        private final byte[][] bytes;
        private final int hash;

        public ByteString(byte[][] bytes) {
            this.bytes = new byte[bytes.length][];
            System.arraycopy(bytes, 0, this.bytes, 0, bytes.length);

            // Pre-compute the hash-code
            int h = 0;
            for (final byte[] bs : this.bytes) {
                for (final byte b : bs) {
                    h = 31 * h + b;
                }
            }

            this.hash = h;

        }

        @Override
        public final int hashCode() {
            return this.hash;
        }

        @Override
        public boolean equals(final Object obj) {
            try {
                final ByteString that = (ByteString) obj;
                if (this.hash != that.hash) return false;
                if (this.bytes.length != that.bytes.length) return false;
                for (int i=0; i<this.bytes.length; ++i) {
                    if (!Arrays.equals(this.bytes[i], that.bytes[i])) return false;
                }
                return true;
            }
            catch (final Exception e) {
                return false;
            }
        }

        @Override
        public String toString() {
            StringBuilder barcodeBuilder = new StringBuilder();
            for(byte[] barcode : bytes){
                barcodeBuilder.append(new String(barcode, StandardCharsets.UTF_8));
            }
            return barcodeBuilder.toString();
        }
    }
}
