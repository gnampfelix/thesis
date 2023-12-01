package org.husonlab.util;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;

import jloda.util.FileUtils;


public class FastKMerIterator implements Closeable, Iterator<byte[]> {
    private static boolean[] isLineContainingSkippableChar = new boolean[128];
    static {
        isLineContainingSkippableChar['\t'] = true;
        isLineContainingSkippableChar['>'] = true;
        isLineContainingSkippableChar['\n'] = true;
        isLineContainingSkippableChar['\r'] = true;
        isLineContainingSkippableChar[' '] = true;
    }

    private final byte[] kmer;
    private final byte[] preloaded_kmer;
    private final InputStreamReader reader;
    private final int k;
    
    private byte nextByte; 
    private boolean isNewSequence;

    public FastKMerIterator(int k, String fileName) throws IOException {
        this.k = k;
        this.kmer = new byte[k];
        this.preloaded_kmer = new byte[k];

        this.reader = new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName));
        this.readUntilSequenceStart();
        this.preloadKmer();
    }

    private void preloadKmer() throws IOException {
        // only prepare the first k-1 letters, user will call hasNext() and
        // next(). Assume: We are at the beginning of the actual sequence,
        // headers are already skipped.

        int i = 0;
        while (i < this.k-1 && this.hasNext()) {
            if(isLineContainingSkippableChar[this.nextByte]) {
                this.skipToNextLine();
                if (this.nextByte == '>') {
                    this.readUntilSequenceStart();
                    i = 0; // the current sequence is not long enough to support the k-mer
                }           
            }
            this.preloaded_kmer[i++] = this.nextByte;
            this.isNewSequence = true;
            this.nextByte = (byte) this.reader.read();
        }
    }

    private void readUntilSequenceStart() throws IOException {  
        if (!this.hasNext())
            return;

        this.nextByte = (byte) this.reader.read();
        while (this.hasNext() && isLineContainingSkippableChar[this.nextByte]) {
            this.skipToNextLine();
        }
    }

    private void skipToNextLine() throws IOException {
        while (this.hasNext() && this.nextByte != '\n') {
            this.nextByte = (byte) this.reader.read();
        }
        if (this.hasNext())
            this.nextByte = (byte) this.reader.read();
    }

    public int getK() {
        return this.k;
    }

    @Override
    public void close() throws IOException {
        this.reader.close();
    }

    @Override
    public boolean hasNext() {
        return this.nextByte != -1;
    }

    @Override
    public byte[] next() {
        if (this.isNewSequence) {
            System.arraycopy(this.preloaded_kmer, 0, this.kmer, 0, this.k - 1);
            this.isNewSequence = false;
        } else {
            System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
        }
        this.kmer[this.k-1] = this.nextByte;
        try {
            this.nextByte = (byte) reader.read();
            while(this.hasNext() && isLineContainingSkippableChar[this.nextByte]) {
                this.isNewSequence = this.nextByte == '>' || this.isNewSequence;
                skipToNextLine();
            }
            if(this.isNewSequence) {
                this.preloadKmer();
            }
            
        } catch (IOException e) {
            this.nextByte = -1;
        }
        
        return this.kmer;
    }
    
}
