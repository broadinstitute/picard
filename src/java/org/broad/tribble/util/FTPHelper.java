package org.broad.tribble.util;

import net.sf.samtools.util.ftp.FTPClient;
import net.sf.samtools.util.ftp.FTPStream;
import net.sf.samtools.util.ftp.FTPUtils;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * @author jacob
 * @date 2013-Apr-11
 */
public class FTPHelper implements URLHelper {

    private URL url;

    public FTPHelper(URL url){
        if(!url.getProtocol().toLowerCase().equals("ftp")){
            throw new IllegalArgumentException("FTPHelper can only be used with ftp protocol, not " + url.getProtocol());
        }
        this.url = url;
    }

    @Override
    public URL getUrl() {
        return this.url;
    }

    @Override
    public long getContentLength() throws IOException {
        return FTPUtils.getContentLength(this.url);
    }

    @Override
    public InputStream openInputStream() throws IOException {
        String file = url.getPath();
        FTPClient ftp = FTPUtils.connect(url.getHost(), url.getUserInfo(), null);
        ftp.pasv();
        ftp.retr(file);
        return new FTPStream(ftp);
    }

    @Override
    public InputStream openInputStreamForRange(long start, long end) throws IOException {
        throw new UnsupportedOperationException("Cannot perform range operations over FTP");
    }

    @Override
    public boolean exists() throws IOException {
        return FTPUtils.resourceAvailable(this.url);
    }
}
