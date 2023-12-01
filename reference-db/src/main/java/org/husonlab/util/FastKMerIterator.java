package org.husonlab.util;

import java.io.Closeable;
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

    private static final byte[] complementTable = new byte[128];
    static {
        complementTable['A'] = 'T';
        complementTable['T'] = 'A';
        complementTable['G'] = 'C';
        complementTable['C'] = 'G';
        complementTable['a'] = 't';
        complementTable['t'] = 'a';
        complementTable['g'] = 'c';
        complementTable['c'] = 'g';
    }

    private final boolean[] isAmbiguousChar = new boolean[128];    
    private final byte[] kmer;
    private final byte[] complement;
    private final byte[] preloaded_kmer;
    private final byte[] preloaded_complement;
    private final InputStreamReader reader;
    private final int k;
    
    private byte nextByte; 
    private boolean isNewSequence;

    public FastKMerIterator(int k, String fileName, boolean skipN) throws IOException {
        this.k = k;
        this.kmer = new byte[k];
        this.preloaded_kmer = new byte[k];
        this.complement = new byte[k];
        this.preloaded_complement = new byte[k];

        this.isAmbiguousChar['N'] = skipN;
        this.isAmbiguousChar['n'] = skipN;
    
        this.reader = new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName));
        this.readUntilSequenceStart();
        this.preloadKmer();
    }


    private void preloadKmer() throws IOException {
        // only prepare the first k-1 letters, user will call hasNext() and
        // next(). Assume: We are at the beginning of the actual sequence,
        // headers are already skipped.
        // Also assume: this.nextByte is at the first character of the new k-mer

        int i = 0;
        while (i < this.k-1 && this.hasNext()) {
            if(isLineContainingSkippableChar[this.nextByte]) {
                this.skipToNextLine();
                if (this.nextByte == '>') {
                    this.readUntilSequenceStart();
                    i = 0; // the current sequence is not long enough to support the k-mer
                }           
            }

            if(isAmbiguousChar[this.nextByte]) {
                this.nextByte = (byte) this.reader.read();
                i = 0; // we need to start over again!
                continue;
            }

            this.preloaded_kmer[i] = this.nextByte;
            this.preloaded_complement[this.k-1-i] = complementTable[this.nextByte];
            i++;
            this.nextByte = (byte) this.reader.read();
        }
        this.isNewSequence = true;
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
            System.arraycopy(this.preloaded_complement, 1, this.complement, 1, this.k-1);
            this.isNewSequence = false;
        } else {
            System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
            System.arraycopy(this.complement, 0, this.complement, 1, this.k - 1);
        }
        this.kmer[this.k-1] = this.nextByte;
        this.complement[0] = complementTable[this.nextByte];
        try {
            this.nextByte = (byte) reader.read();
            // Skip ambiguous characters ("N"/"n"), dicard all k-mers
            while(this.hasNext() && isAmbiguousChar[this.nextByte]) {
                this.nextByte = (byte) reader.read();
                this.isNewSequence = true;
            }
            while(this.hasNext() && isLineContainingSkippableChar[this.nextByte]) {
                this.isNewSequence = this.nextByte == '>' || this.isNewSequence;
                skipToNextLine();
            }

            // If we start a new sequence (i.e. new entry in fasta file or after
            // an amb. character), we need to preload all of the next k-1
            // characters
            if(this.isNewSequence) {
                this.preloadKmer();
            }
            
        } catch (IOException e) {
            this.nextByte = -1;
        }
        
        return this.kmer;
    }

    /**
     * Returns the reverse complement to the kmer returned by next(). Only correct if called after next()!
     * @return
     */
    public byte[] getComplement() {
        return this.complement;
    }
    
}
