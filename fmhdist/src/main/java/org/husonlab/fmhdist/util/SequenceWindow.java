package org.husonlab.fmhdist.util;

import java.util.ArrayList;
import java.util.List;

public class SequenceWindow {
    private int windowSize;
    private int startPosition;
    private List<byte[]> kmers;
    private int uniqueKmersInWindow;

    public SequenceWindow(int startPosition, int windowSize, List<byte[]> kmers) {
        this.windowSize = windowSize;
        this.startPosition = startPosition;
        this.kmers = kmers;
        this.uniqueKmersInWindow = 0;
        for (byte[] kmer : kmers) {
            if (isUnique(kmer)) {
                uniqueKmersInWindow++;
            }
        }
    }

    private boolean isUnique(byte[] kmer) {
        String kmerString = new String(kmer);
        for (byte[] compare : this.kmers) {
            if (compare == kmer) {
                continue;
            }
            if (new String(compare).equals(kmerString)) {
                return false;
            }
        }
        return true;
    }

    public int getWindowSize() {
        return this.windowSize;
    }

    public int getStartPosition() {
        return this.startPosition;
    }

    public List<byte[]> getKmers() {
        return new ArrayList<>(this.kmers);
    }

    public int getUniqueKmersInWindow() {
        return this.uniqueKmersInWindow;
    }

    public String toString() {
        return String.format("%d,%d,%d", this.startPosition, this.kmers.size(), this.uniqueKmersInWindow);
    }
}
