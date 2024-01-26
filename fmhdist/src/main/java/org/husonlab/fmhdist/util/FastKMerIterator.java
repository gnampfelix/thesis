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
 * With each new record in the file, a new k-mer is created, i.e. this iterator
 * won't return k-mers spanning multiple records.
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
    private final byte[] kmer;
    private final byte[] complement;
    private final byte[] preloaded_kmer;
    private final byte[] preloaded_complement;
    private final InputStreamReader reader;
    private final int k;
    
    private byte nextByte;
    private boolean hasPreloaded;

    private int recordIndexInFile = 0;
    private int skippedKmersInFile = 0;
    private int skippedKmersInRecord = 0;
    private int sequenceIndexInRecord = 0;
    private int sequenceIndexInFile = 0;

    private int preloadedRecordIndexInFile = 0;
    private int preloadedSkippedKmersInFile = 0;
    private int preloadedSkippedKmersInRecord = 0;
    private int preloadedSequenceIndexInRecord = 0;
    private int preloadedSequenceIndexInFile = 0;
 
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
        this(k, new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(fileName)), skipN);        
    }

    /**
     * Create a new Iterator to extract the kmers from the given file.
     * @param k the size of the k-mer
     * @param reader the InputStreamReader to read the sequence from.
     * @param skipN indicate if k-mers containing the letter "N" or "n" should
     * be skipped. This is typically the case if the base is ambiguous for
     * genomic sequences.
     * @throws IOException
     */
    public FastKMerIterator(int k, InputStreamReader reader, boolean skipN) throws IOException {
        this.k = k;
        this.kmer = new byte[k];
        this.preloaded_kmer = new byte[k];
        this.complement = new byte[k];
        this.preloaded_complement = new byte[k];

        this.isAmbiguousChar['N'] = skipN;
        this.isAmbiguousChar['n'] = skipN;

        this.reader = reader;
        this.nextByte = (byte) reader.read();
        this.preloadedRecordIndexInFile = -1;
        this.preload();
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

    private boolean isSequenceChar() {
        return !isLineContainingSkippableChar[this.nextByte];
    }

    private boolean isHeaderStart() {
        return this.nextByte == '>';
    }

    private void skipToNextLine() throws IOException {
        while (this.hasNext() && this.nextByte != '\n') {
            this.nextByte = (byte) this.reader.read();
        }
        if (this.hasNext())
            this.nextByte = (byte) this.reader.read();
    }

    private void handleSequenceStart() throws IOException {
        this.preloadedRecordIndexInFile++;
        this.preloadedSequenceIndexInRecord = 0;
        this.preloadedSkippedKmersInRecord = 0;
        skipToNextLine();
    }

    /**
     * Assuming that "this.nextByte" has NOT been handled yet, this function
     * reads as many bytes as needed to fill in the first k-1 bytes of a NEW
     * k-mer. This is needed if
     * 1. A new sequence starts
     * 2. An existing k-mer cannot be extended because the next character is
     *    ambiguous. In this case, the caller MUST already have loaded the byte
     *    AFTER the ambiguous nucleotide.
     *
     * If this function needs to skip k-mers because it encounters ambiguous
     * nucleotides, it updates the according counters.
     * @throws IOException
     */
    private void preload() throws IOException {
        this.hasPreloaded = true;
        int i = 0;
        while(this.hasNext() && i < this.k - 1) {
            if(isSequenceChar()) {
                if (isAmbiguousChar[this.nextByte]) {
                    // we need to discard the previous i k-mers (0-based, thus
                    // +1)
                    this.preloadedSkippedKmersInFile += i + 1;
                    this.preloadedSkippedKmersInRecord += i + 1;
                    i = 0;
                    this.nextByte = (byte) this.reader.read();
                    continue;
                }
                this.preloaded_kmer[i] = toUpperTable[this.nextByte];
                this.preloaded_complement[this.k - i - 1] = toUpperTable[complementTable[this.nextByte]];
                this.nextByte = (byte) this.reader.read();
                i++;
            } else {
                if (isHeaderStart()) {
                    handleSequenceStart();
                    i = 0;
                } else {
                    skipToNextLine();
                }
            }
        }
    }

    private void copyIndices() {
        this.recordIndexInFile = this.preloadedRecordIndexInFile;
        this.sequenceIndexInFile = this.preloadedSequenceIndexInFile;
        this.sequenceIndexInRecord = this.preloadedSequenceIndexInRecord;
        this.skippedKmersInFile = this.preloadedSkippedKmersInFile;
        this.skippedKmersInRecord = this.preloadedSkippedKmersInRecord;
    }

    @Override
    public byte[] next() {
        // First, finalize current k-mer
        if (this.hasPreloaded) {
            System.arraycopy(this.preloaded_kmer, 0, this.kmer, 0, this.k - 1);
            System.arraycopy(this.preloaded_complement, 1, this.complement, 1, this.k-1);
            this.hasPreloaded = false;
        } else {
            System.arraycopy(this.kmer, 1, this.kmer, 0, this.k - 1);
            System.arraycopy(this.complement, 0, this.complement, 1, this.k - 1);
        }

        this.copyIndices();
        this.kmer[this.k-1] = toUpperTable[this.nextByte];
        this.complement[0] = toUpperTable[complementTable[this.nextByte]];
        
        // Now, k-mer is finished. Time to prepare the next one!
        try {
            this.nextByte = (byte) reader.read();
            this.preloadedSequenceIndexInFile++;
            this.preloadedSequenceIndexInRecord++;

            boolean forceNextIteration = true;
            boolean needsPreload = false;

            // There are four possible cases to consider:
            // 1. The next character is a valid sequence character
            // 2. The next character is a ambiguous sequence character
            // 3. The next character is the start of a new header
            // 4. The next character is a \n
            while(hasNext() & forceNextIteration){
                forceNextIteration = false;
                if (isSequenceChar()) {
                    if (isAmbiguousChar[this.nextByte]){
                        // We need to skip the next k k-mers
                        this.preloadedSkippedKmersInFile += this.k;
                        this.preloadedSkippedKmersInRecord += this.k;
                        this.nextByte = (byte) reader.read();
                        needsPreload = true;
                    } 
                } else {
                    if (isHeaderStart()){
                        handleSequenceStart();
                        needsPreload = true;
                    } else {
                        skipToNextLine();
                        forceNextIteration = true;
                    }
                }

                if(needsPreload) {
                    this.preload();
                    // We do not need to break because in above, it is
                    // impossible to needPreload AND forceNextIteration, but to
                    // be explicit here:
                    break;
                }
            }
        } catch (IOException e) {
            this.nextByte = -1;
        }

        // At last: finished!
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

    /**
     * Returns the coordinates of the current k-mer. Only correct if called
     * after next()!
     * @return
     */
    public KMerCoordinates getCoordinates() {
        return new KMerCoordinates(
            this.recordIndexInFile, 
            this.sequenceIndexInFile, 
            this.sequenceIndexInRecord,
            this.sequenceIndexInFile + this.skippedKmersInFile,
            this.sequenceIndexInRecord + this.skippedKmersInRecord, 
            this.kmer);
    }    
}
