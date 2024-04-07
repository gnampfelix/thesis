package org.husonlab.fmhdist.util;

import java.io.IOException;

public interface FileConsumer {
    public boolean isReady();
    public String getLine() throws IOException;
}
