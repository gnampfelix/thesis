package org.husonlab.fmhdist.util;

import java.util.Iterator;

public class KMerIterator implements Iterator<byte[]> {

    private Iterator<byte[]> fastaProvider;
    private final int k;
    private byte[] bytes;
    private byte[] kmer;
    private byte[] nextKmer;
    private byte[] remainingBasesInRow;
    private int pos = 0;
    private boolean done;

    // NCBI FASTA Files have a typical line length of 80, so 1000 is more then enough
    private static final int BUFFER_SIZE = 1000;

    public KMerIterator(Iterator<byte[]> fastaProvider, int k) {
        this.k = k;
        this.fastaProvider = fastaProvider;  
        this.kmer = new byte[k];
        this.nextKmer = new byte[k];
        this.remainingBasesInRow = new byte[k-1];
        this.bytes = new byte[BUFFER_SIZE];
        if (this.fastaProvider.hasNext()) {
            this.done = !this.extractKmer();
        }
    }

    // Assumption: call if bytes[pos + k] == \n or \r or \0
    private boolean fetch() {
        if (!this.fastaProvider.hasNext()) {
            return false;
        }
        byte[] nextLine = this.fastaProvider.next();
        
        if (nextLine == null) {
            // this is an error in the fastaProvider - there should not be 
            // a null result! Don't do anything here to fail loudly!
        }

        // new sequence - start new k-mers, remaining bases are discarded
        if (nextLine[0] == '>') {
            if (!this.fastaProvider.hasNext()) {
                // invalid fasta - this shouldn't happen - header without
                // sequence
                return false;
            }
            nextLine = this.fastaProvider.next();
            if (nextLine.length + this.k > bytes.length) {
                // already prepare buffer for remaining k letters for next
                // iteration
                this.bytes = new byte[nextLine.length + this.k]; 
            }         
            System.arraycopy(nextLine, 0, bytes, 0, nextLine.length);          
        } else {
            System.arraycopy(bytes, pos, remainingBasesInRow, 0, this.k - 1);
            if (nextLine.length + this.k > bytes.length) {
                this.bytes = new byte[nextLine.length + this.k];
            }            
            System.arraycopy(remainingBasesInRow, 0, bytes, 0, this.k-1); //copy remaining bases, pos+k == \n or \r or 0!
            System.arraycopy(nextLine, 0, bytes, this.k-1, nextLine.length); //copy new bases
        }
        this.pos = 0;
        return true;
    }

    private boolean extractKmer() {
        if (
            (this.pos + this.k - 1) >= this.bytes.length 
            || this.bytes[this.pos + this.k - 1] == 0
            || this.bytes[this.pos + this.k - 1] == '\n'
            || this.bytes[this.pos + this.k - 1] == '\r') {
            if (!this.fetch()) {
                return false;
            }
        }
        System.arraycopy(this.bytes, this.pos++, this.nextKmer, 0, this.k);
        return true;
    }


    @Override
    public boolean hasNext() {
        return !this.done;
    }

    @Override
    public byte[] next() {
        byte[] tmp = this.kmer;
        this.kmer = this.nextKmer;
        this.nextKmer = tmp;
        this.done = !this.extractKmer();
        return this.kmer;
    }

    public int getK() {
        return this.k;
    }
    
}
