package org.broad.tribble.util;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Implementation of URLHelper designed for remote resources.
 * Determines appropriate helper based on protocol (ftp/http),
 * and then delegates.
 *
 * @author jacob
 * @date 2013-Apr-11
 */
public class RemoteURLHelper implements URLHelper {

    private URLHelper wrappedHelper;

    public RemoteURLHelper(URL url){
        String protocol = url.getProtocol().toLowerCase();
        if(protocol.startsWith("http")){
            this.wrappedHelper = new HTTPHelper(url);
        }else if(protocol.startsWith("ftp")){
            this.wrappedHelper = new FTPHelper(url);
        }else{
            throw new IllegalArgumentException("Unable to create helper for url with protocol " + protocol);
        }
    }

    @Override
    public URL getUrl() {
        return this.wrappedHelper.getUrl();
    }

    @Override
    public long getContentLength() throws IOException {
        return this.wrappedHelper.getContentLength();
    }

    @Override
    public InputStream openInputStream() throws IOException {
        return this.wrappedHelper.openInputStream();
    }

    @Override
    public InputStream openInputStreamForRange(long start, long end) throws IOException {
        return this.wrappedHelper.openInputStreamForRange(start, end);
    }

    @Override
    public boolean exists() throws IOException {
        return this.wrappedHelper.exists();
    }
}
