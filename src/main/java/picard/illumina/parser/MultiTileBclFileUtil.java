package picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import picard.illumina.parser.fakers.BclIndexFaker;
import picard.illumina.parser.fakers.MultiTileBclFileFaker;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * NextSeq-style bcl's have all tiles for a cycle in a single file.
 */
public class MultiTileBclFileUtil extends ParameterizedFileUtil {
    final File basecallLaneDir;
    final File bci;
    final TileIndex tileIndex;
    final CycleIlluminaFileMap cycleFileMap = new CycleIlluminaFileMap();

    MultiTileBclFileUtil(final File basecallLaneDir, final int lane) {
        // Since these file names do not contain lane number, first two args to ctor are the same.
        super("^(\\d{4}).bcl.bgzf$", ".bcl.bgzf", basecallLaneDir,
                new MultiTileBclFileFaker(), lane);
        this.basecallLaneDir = basecallLaneDir;
        bci = new File(basecallLaneDir, "s_" + lane + ".bci");
        // Do this once rather than when deciding if these files exist and again later.
        final File[] cycleFiles = IOUtil.getFilesMatchingRegexp(base, matchPattern);
        if (bci.exists()) {
            tileIndex = new TileIndex(bci);
            if (cycleFiles != null) {
                for (final File file : cycleFiles) {
                    final String fileName = file.getName();
                    final String cycleNum = fileName.substring(0, fileName.indexOf('.'));
                    final IlluminaFileMap fileMap = new IlluminaFileMap();
                    for (final Integer tile : tileIndex.getTiles()) {
                        fileMap.put(tile, file);
                    }
                    cycleFileMap.put(Integer.valueOf(cycleNum), fileMap);
                }
            }
        } else {
            tileIndex = null;
        }

    }

    public CycleIlluminaFileMap getFiles(final List<Integer> tiles, final int[] cycles) {
        // Filter input list of cycles according to which actually exist
        final ArrayList<Integer> goodCycleList = new ArrayList<Integer>(cycles.length);
        for (final int cycle : cycles) {
            if (cycleFileMap.containsKey(cycle)) {
                goodCycleList.add(cycle);
            }
        }
        // Ensure cycles are sorted.
        Collections.sort(goodCycleList);
        final int[] goodCycles = new int[goodCycleList.size()];
        for (int i = 0; i < goodCycles.length; ++i) {
            goodCycles[i] = goodCycleList.get(i);
        }

        // Create the map.
        final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();
        if (goodCycles.length > 0) {
            for (final int cycle : goodCycles) {
                final IlluminaFileMap fileMap = cycleFileMap.get(cycle).keep(tiles);
                cycledMap.put(cycle, fileMap);
            }
        }
        return cycledMap;
    }

    @Override
    public boolean filesAvailable() {
        return bci.exists() && !cycleFileMap.isEmpty();
    }

    @Override
    public List<Integer> getTiles() {
        if (tileIndex == null) {
            return Collections.emptyList();
        }
        return tileIndex.getTiles();
    }

    @Override
    public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
        if (tileIndex == null) {
            return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
        }
        final List<String> ret = tileIndex.verify(expectedTiles);
        for (final int expectedCycle : expectedCycles) {
            if (!cycleFileMap.containsKey(expectedCycle)) {
                ret.add(expectedCycle + ".bcl.bgzf not found in " + base);
            }
        }
        return ret;
    }

    @Override
    public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        if (tileIndex == null) {
            return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
        }
        final List<String> ret = tileIndex.verify(expectedTiles);
        for (final int expectedCycle : expectedCycles) {
            if (!cycleFileMap.containsKey(expectedCycle)) {
                //we need to fake a bci file for the tile index
                final BclIndexFaker bciFileFaker = new BclIndexFaker();
                final MultiTileBclFileFaker bclFileFaker = new MultiTileBclFileFaker();
                try {
                    bciFileFaker.fakeBciFile(new File(basecallLaneDir, String.format("%04d", expectedCycle) + ".bcl.bgzf.bci"), tileIndex);
                    bclFileFaker.fakeMultiTileBclFile(new File(basecallLaneDir, String.format("%04d", expectedCycle) + ".bcl.bgzf"), tileIndex);
                } catch (final IOException e) {
                    return Collections.singletonList("Could not create tile index file: " + bci.getAbsolutePath());
                }
                ret.add(expectedCycle + ".bcl.bgzf not found in " + base);
            }
        }
        return ret;
    }
}

