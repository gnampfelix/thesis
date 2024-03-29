package org.husonlab.fmhdist.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import jloda.util.FileUtils;

public class LineKMerIterator implements KMerIterator {
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

    private static final byte[] toUpperTable = new byte[128];
    static {
        toUpperTable['A'] = 'A';
        toUpperTable['T'] = 'T';
        toUpperTable['G'] = 'G';
        toUpperTable['C'] = 'C';
        toUpperTable['a'] = 'A';
        toUpperTable['t'] = 'T';
        toUpperTable['g'] = 'G';
        toUpperTable['c'] = 'C';
        toUpperTable['n'] = 'N';
        toUpperTable['N'] = 'N';
    }

    private final boolean[] isAmbiguousChar = new boolean[128];

    private int k;
    private byte[] kmer;
    private byte[] kmerReverseComplement;
    private byte[] preloadedKmer;
    private byte[] preloadedKmerReverseComplement;

    private byte[] nextLine;
    private byte[] currentLine;
    private int linePointer;

    private boolean isEOF;
    private boolean isPreloaded;
    private BufferedReader reader;

    public LineKMerIterator(int k, BufferedReader reader, boolean skipN) throws IOException {
        this.k = k;
        this.kmer = new byte[k];
        this.kmerReverseComplement = new byte[k];
        this.preloadedKmer = new byte[k];
        this.preloadedKmerReverseComplement = new byte[k];
        
        this.isEOF = false;
        this.isPreloaded = false;

        this.isAmbiguousChar['N'] = skipN;
        this.isAmbiguousChar['n'] = skipN;

        this.linePointer = 0;
        this.reader = reader;
        String current = reader.readLine();
        String next  = reader.readLine();
        
        if (current == null || next == null) {
            throw new IOException("file is too short, valid FASTA files have at least two lines");
        }
        this.currentLine = current.getBytes();
        this.nextLine = next.getBytes();

        this.preload();
    }

    public LineKMerIterator(int k, String fileName, boolean skipN) throws IOException {
       this(k, new BufferedReader(new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName))), skipN);
    }

    private boolean isSequenceChar() {
        return !isLineContainingSkippableChar[this.currentLine[linePointer]];
    }

    private void feedLine() throws IOException {
        if (this.isEOF) {
            return;
        }
        this.currentLine = this.nextLine;
        this.linePointer = 0;
        String next = reader.readLine();
        if (next == null) {
            this.isEOF = true;
        } else {
            this.nextLine = next.getBytes();
        }
    }

    private void moveCursor() throws IOException {
        if (++this.linePointer >= this.currentLine.length) {
            this.feedLine();
        }
    }

    private void preload() throws IOException{
        int i = 0;
        while(this.hasNext() && i < this.k - 1) {
            if(isSequenceChar()) {
                if (isAmbiguousChar[this.currentLine[linePointer]]) {
                    i = 0;
                    this.moveCursor();
                    continue;
                }
                this.preloadedKmer[i] = toUpperTable[this.currentLine[linePointer]];
                this.preloadedKmerReverseComplement[this.k - i - 1] = toUpperTable[complementTable[this.currentLine[linePointer]]];
                this.moveCursor();
                i++;
            } else {
                // no need to check header start explicitely - if NOT a valid sequence character, will skip the current line.
                this.feedLine();
                i = 0;
            }
        }
        this.isPreloaded = true;
    }

    @Override
    public boolean hasNext() {
        return !this.isEOF || this.linePointer < this.currentLine.length;
    }

    @Override
    public byte[] next() {
        if (this.isPreloaded) {
            System.arraycopy(this.preloadedKmer, 0, this.kmer, 0, this.k - 1);
            System.arraycopy(this.preloadedKmerReverseComplement, 1, this.kmerReverseComplement, 1, this.k-1);
            this.isPreloaded = false;
        } else {
            System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
            System.arraycopy(this.kmerReverseComplement, 0, this.kmerReverseComplement, 1, this.k - 1);
        }
        this.kmer[this.k-1] = toUpperTable[this.currentLine[linePointer]];
        this.kmerReverseComplement[0] = toUpperTable[complementTable[this.currentLine[linePointer]]];

        try {
            this.moveCursor();
            if (this.hasNext()) {
                if (this.isSequenceChar() && isAmbiguousChar[this.currentLine[this.linePointer]]) {
                    this.moveCursor();
                    this.preload();
                } else if (!this.isSequenceChar()) {
                    this.feedLine();
                    this.preload();
                }
            }
        } catch (IOException e) {
            this.isEOF = true;
            this.linePointer = this.currentLine.length;
        }
        return this.kmer;
    }

    @Override
    public void close() throws IOException {
        this.reader.close();
    }

    @Override
    public int getK() {
        return this.k;
    }

    @Override
    public byte[] getReverseComplement() {
        return this.kmerReverseComplement;
    }

    @Override
    public KMerCoordinates getCoordinates() {
        //TODO: actually count the values
        return new KMerCoordinates(
            0,
            0,
            0,
            0,
            0,
            this.kmer,
            0
        );
    }
    
}
