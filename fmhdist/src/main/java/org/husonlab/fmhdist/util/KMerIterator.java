package org.husonlab.fmhdist.util;

import java.io.Closeable;
import java.util.Iterator;

public interface KMerIterator extends Iterator<byte[]>, Closeable {
    public int getK();
    public byte[] getReverseComplement();
    public KMerCoordinates getCoordinates();
}
