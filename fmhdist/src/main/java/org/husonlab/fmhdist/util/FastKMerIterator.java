package org.husonlab.fmhdist.util;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;

import jloda.util.FileUtils;

/**
 * Iterator to extract all valid k-mers from a fasta file (raw, zipped or
 * gzipped), local file system or ftp url.
 *
 * The iterator can skip ambiguous bases (N/n) and calculate the reverse
 * complement for the given k-mer.
 *
 * If the fasta file contains multiple records, there won't be any k-mers in the
 * result of this iterator that span multiple records.
 */
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

    /**
     * Create a new Iterator to extract the kmers from the given file.
     * @param k the size of the k-mer
     * @param fileName the path to the file to fetch
     * @param skipN indicate if k-mers containing the letter "N" or "n" should
     * be skipped. This is typically the case if the base is ambiguous for
     * genomic sequences.
     * @throws IOException
     */
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


    /**
     * Assuming this.nextByte is the start of a new k-mer, this function reads
     * at least the next k-1 letters to preload the next k-mer. If those k-1
     * bytes contain any character that needs special treatment (new-line,
     * beginning of new record, amb. character), this function will read so many
     * bytes until preloaded_kmer contains k-1 bytes (or EOF)
     * @throws IOException
     */
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

    /**
     * Read bytes until the next character is neither whitespace or part of a
     * header-line. After this function returns, this.nextByte points at the
     * first byte of the sequence.
     * @throws IOException
     */
    private void readUntilSequenceStart() throws IOException {  
        if (!this.hasNext())
            return;

        this.nextByte = (byte) this.reader.read();
        while (this.hasNext() && isLineContainingSkippableChar[this.nextByte]) {
            this.skipToNextLine();
        }
    }

    /**
     * Read bytes until the next line. After this function returns,
     * this.nextByte points at the first byte of the new line.
     * @throws IOException
     */
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
     * Returns the reverse complement to the kmer returned by next(). Only
     * correct if called after next()!
     * @return
     */
    public byte[] getReverseComplement() {
        return this.complement;
    }
    
}
