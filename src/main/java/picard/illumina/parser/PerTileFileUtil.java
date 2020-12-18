package picard.illumina.parser;

import picard.illumina.parser.fakers.FileFaker;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PerTileFileUtil extends ParameterizedFileUtil {
    final IlluminaFileMap fileMap;

    public PerTileFileUtil(final String extension, final File base,
                           final FileFaker faker, final int lane) {
        this(extension, base, faker, lane, DefaultSkipEmptyFiles);
    }

    public PerTileFileUtil(final String extension, final File base,
                           final FileFaker faker, final int lane, final boolean skipEmptyFiles) {
        super(true, extension, base, faker, lane, skipEmptyFiles);
        this.fileMap = getTiledFiles(base, matchPattern);
        if (!fileMap.isEmpty()) {
            this.tiles = new int[this.fileMap.keySet().size()];
            int i = 0;
            for(int tile: this.fileMap.keySet()) {
                this.tiles[i] = tile;
                i++;
            }
        } else {
            this.tiles = new int[0];
        }
    }

    @Override
    public boolean filesAvailable() {
        return !fileMap.isEmpty();
    }

    public IlluminaFileMap getFiles() {
        return fileMap;
    }

    public IlluminaFileMap getFiles(final int[]  tiles) {
        return fileMap.keep(tiles);
    }

    @Override
    public List<String> verify(final int[]  expectedTiles, final int[] expectedCycles) {
        return verifyPerTile(this.base, expectedTiles);
    }

    List<String> verifyPerTile(File baseDir, int[]  expectedTiles) {
        final List<String> failures = new LinkedList<>();
        if (!baseDir.exists()) {
            failures.add("Base directory(" + baseDir.getAbsolutePath() + ") does not exist!");
        } else {
            Set<Integer> expectedSet = IntStream.of(expectedTiles).boxed().collect(Collectors.toCollection(HashSet::new));
            int[] missing = IntStream.of(tiles).filter(val -> !expectedSet.contains(val)).toArray();
            if(missing.length > 0){
                failures.add("Missing tile " + Arrays.toString(missing) + " for file type " + extension + ".");
            }
        }
        return failures;
    }

    @Override
    public List<String> fakeFiles(int[] expectedTiles, final int[] cycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        final List<String> failures = new LinkedList<>();
        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
            for (final Integer tile : expectedTiles) {
                if (Arrays.stream(tiles).noneMatch(query ->query == tile) || fileMap.get(tile).length() == 0) {
                    //create a new file of this type
                    try {
                        faker.fakeFile(base, tile, lane, extension);
                    } catch (final IOException e) {
                        failures.add(String.format("Could not create fake file %s: %s", fileMap.get(tile),
                                e.getMessage()));
                    }

                }
            }
        }
        return failures;
    }

}
