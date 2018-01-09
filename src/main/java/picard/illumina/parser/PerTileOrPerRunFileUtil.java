package picard.illumina.parser;

import picard.illumina.parser.fakers.FileFaker;

import java.io.File;
import java.util.List;
import java.util.regex.Pattern;

public class PerTileOrPerRunFileUtil extends PerTileFileUtil {
    private final File runFile;

    public PerTileOrPerRunFileUtil(String extension, File base, FileFaker faker, int lane) {
        super(extension, base, faker, lane);
        Pattern runFileMatchPattern = Pattern.compile("^s" + extension + "$");
        runFile = getRunFile(base.getParentFile(), runFileMatchPattern);
    }

    @Override
    public boolean filesAvailable() {
        return super.filesAvailable() || runFile != null;
    }

    public File getRunFile() {
        return runFile;
    }

    @Override
    public void setTiles(List<Integer> tiles){
        super.setTiles(tiles);
        //now fill in the file map
        if(fileMap.isEmpty()){
            tiles.forEach(i -> fileMap.put(i, runFile));
        }
    }
}
