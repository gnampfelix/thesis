package org.husonlab.fmhdist.utils;


import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.Matchers.hasSize;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.husonlab.fmhdist.util.FileProducer;
import org.husonlab.fmhdist.util.InMemoryKMerIterator;
import org.husonlab.fmhdist.util.KMerIterator;
import org.husonlab.fmhdist.util.LineKMerIteratorConsumer;
import org.junit.Ignore;
import org.junit.Test;



public class LineKMerIteratorConsumerTests {
    @Test
    public void shouldIterateKmers() throws IOException {
        FileProducer p = new FileProducer(1);
        final ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(() -> {p.run();});
        try(KMerIterator km = new LineKMerIteratorConsumer(21, "src/test/resources/fastaWith1Seq.fasta", p, 0, true)) {

            List<String>kmers = new ArrayList<>();
            while(km.hasNext()) {
                byte[] next = km.next();
                kmers.add(new String(next));
            }
            assertThat(kmers, contains("ACTCGTATGAACTTTGACTGG",
            "CTCGTATGAACTTTGACTGGT",
            "TCGTATGAACTTTGACTGGTT",
            "CGTATGAACTTTGACTGGTTT",
            "GTATGAACTTTGACTGGTTTT",
            "TATGAACTTTGACTGGTTTTT",
            "ATGAACTTTGACTGGTTTTTG",
            "TGAACTTTGACTGGTTTTTGG",
            "GAACTTTGACTGGTTTTTGGG",
            "AACTTTGACTGGTTTTTGGGG",
            "ACTTTGACTGGTTTTTGGGGC",
            "CTTTGACTGGTTTTTGGGGCG",
            "TTTGACTGGTTTTTGGGGCGC",
            "TTGACTGGTTTTTGGGGCGCG",
            "TGACTGGTTTTTGGGGCGCGA",
            "GACTGGTTTTTGGGGCGCGAG",
            "ACTGGTTTTTGGGGCGCGAGA",
            "CTGGTTTTTGGGGCGCGAGAG",
            "TGGTTTTTGGGGCGCGAGAGT",
            "GGTTTTTGGGGCGCGAGAGTT",
            "GTTTTTGGGGCGCGAGAGTTT",
            "TTTTTGGGGCGCGAGAGTTTG",
            "TTTTGGGGCGCGAGAGTTTGG",
            "TTTGGGGCGCGAGAGTTTGGT",
            "TTGGGGCGCGAGAGTTTGGTT",
            "TGGGGCGCGAGAGTTTGGTTT",
            "GGGGCGCGAGAGTTTGGTTTG",
            "GGGCGCGAGAGTTTGGTTTGG",
            "GGCGCGAGAGTTTGGTTTGGA",
            "GCGCGAGAGTTTGGTTTGGAT",
            "CGCGAGAGTTTGGTTTGGATG",
            "GCGAGAGTTTGGTTTGGATGA",
            "CGAGAGTTTGGTTTGGATGAA",
            "GAGAGTTTGGTTTGGATGAAA",
            "AGAGTTTGGTTTGGATGAAAC",
            "GAGTTTGGTTTGGATGAAACG",
            "AGTTTGGTTTGGATGAAACGC",
            "GTTTGGTTTGGATGAAACGCA",
            "TTTGGTTTGGATGAAACGCAC",
            "TTGGTTTGGATGAAACGCACC",
            "TGGTTTGGATGAAACGCACCC",
            "GGTTTGGATGAAACGCACCCG",
            "GTTTGGATGAAACGCACCCGC",
            "TTTGGATGAAACGCACCCGCT",
            "TTGGATGAAACGCACCCGCTA",
            "TGGATGAAACGCACCCGCTAT",
            "GGATGAAACGCACCCGCTATC",
            "GATGAAACGCACCCGCTATCG",
            "ATGAAACGCACCCGCTATCGA",
            "TGAAACGCACCCGCTATCGAT",
            "GAAACGCACCCGCTATCGATA",
            "AAACGCACCCGCTATCGATAG",
            "AACGCACCCGCTATCGATAGA",
            "ACGCACCCGCTATCGATAGAT",
            "CGCACCCGCTATCGATAGATA",
            "GCACCCGCTATCGATAGATAG",
            "CACCCGCTATCGATAGATAGT",
            "ACCCGCTATCGATAGATAGTT",
            "CCCGCTATCGATAGATAGTTT",
            "CCGCTATCGATAGATAGTTTA",
            "CGCTATCGATAGATAGTTTAC",
            "GCTATCGATAGATAGTTTACC",
            "CTATCGATAGATAGTTTACCT",
            "TATCGATAGATAGTTTACCTT",
            "ATCGATAGATAGTTTACCTTT",
            "TCGATAGATAGTTTACCTTTC",
            "CGATAGATAGTTTACCTTTCT",
            "GATAGATAGTTTACCTTTCTC",
            "ATAGATAGTTTACCTTTCTCT",
            "TAGATAGTTTACCTTTCTCTT",
            "AGATAGTTTACCTTTCTCTTT",
            "GATAGTTTACCTTTCTCTTTT",
            "ATAGTTTACCTTTCTCTTTTT",
            "TAGTTTACCTTTCTCTTTTTT",
            "AGTTTACCTTTCTCTTTTTTT",
            "GTTTACCTTTCTCTTTTTTTG",
            "TTTACCTTTCTCTTTTTTTGA",
            "TTACCTTTCTCTTTTTTTGAA",
            "TACCTTTCTCTTTTTTTGAAA",
            "ACCTTTCTCTTTTTTTGAAAA",
            "CCTTTCTCTTTTTTTGAAAAC",
            "CTTTCTCTTTTTTTGAAAACG",
            "TTTCTCTTTTTTTGAAAACGC",
            "TTCTCTTTTTTTGAAAACGCA",
            "TCTCTTTTTTTGAAAACGCAC",
            "CTCTTTTTTTGAAAACGCACC",
            "TCTTTTTTTGAAAACGCACCC",
            "CTTTTTTTGAAAACGCACCCG",
            "TTTTTTTGAAAACGCACCCGC",
            "TTTTTTGAAAACGCACCCGCT",
            "TTTTTGAAAACGCACCCGCTA",
            "TTTTGAAAACGCACCCGCTAT",
            "TTTGAAAACGCACCCGCTATC",
            "TTGAAAACGCACCCGCTATCG",
            "TGAAAACGCACCCGCTATCGA",
            "GAAAACGCACCCGCTATCGAT",
            "AAAACGCACCCGCTATCGATA",
            "AAACGCACCCGCTATCGATAG",
            "AACGCACCCGCTATCGATAGA",
            "ACGCACCCGCTATCGATAGAT",
            "CGCACCCGCTATCGATAGATA",
            "GCACCCGCTATCGATAGATAG",
            "CACCCGCTATCGATAGATAGG",
            "ACCCGCTATCGATAGATAGGG",
            "CCCGCTATCGATAGATAGGGT",
            "CCGCTATCGATAGATAGGGTT",
            "CGCTATCGATAGATAGGGTTA",
            "GCTATCGATAGATAGGGTTAC",
            "CTATCGATAGATAGGGTTACC",
            "TATCGATAGATAGGGTTACCC",
            "ATCGATAGATAGGGTTACCCT",
            "TCGATAGATAGGGTTACCCTT",
            "CGATAGATAGGGTTACCCTTT",
            "GATAGATAGGGTTACCCTTTC",
            "ATAGATAGGGTTACCCTTTCT",
            "TAGATAGGGTTACCCTTTCTC",
            "AGATAGGGTTACCCTTTCTCT",
            "GATAGGGTTACCCTTTCTCTT",
            "ATAGGGTTACCCTTTCTCTTT",
            "TAGGGTTACCCTTTCTCTTTG",
            "AGGGTTACCCTTTCTCTTTGA",
            "GGGTTACCCTTTCTCTTTGAG",
            "GGTTACCCTTTCTCTTTGAGG",
            "GTTACCCTTTCTCTTTGAGGT",
            "TTACCCTTTCTCTTTGAGGTT",
            "TACCCTTTCTCTTTGAGGTTT",
            "ACCCTTTCTCTTTGAGGTTTT",
            "CCCTTTCTCTTTGAGGTTTTT",
            "CCTTTCTCTTTGAGGTTTTTG",
            "CTTTCTCTTTGAGGTTTTTGG",
            "TTTCTCTTTGAGGTTTTTGGG",
            "TTCTCTTTGAGGTTTTTGGGA",
            "TCTCTTTGAGGTTTTTGGGAA",
            "CTCTTTGAGGTTTTTGGGAAA",
            "TCTTTGAGGTTTTTGGGAAAA",
            "CTTTGAGGTTTTTGGGAAAAC",
            "TTTGAGGTTTTTGGGAAAACG",
            "TTGAGGTTTTTGGGAAAACGC",
            "TGAGGTTTTTGGGAAAACGCA",
            "GAGGTTTTTGGGAAAACGCAC",
            "AGGTTTTTGGGAAAACGCACC",
            "GGTTTTTGGGAAAACGCACCT",
            "GTTTTTGGGAAAACGCACCTG",
            "TTTTTGGGAAAACGCACCTGC",
            "TTTTGGGAAAACGCACCTGCT",
            "TTTGGGAAAACGCACCTGCTA",
            "TTGGGAAAACGCACCTGCTAT",
            "TGGGAAAACGCACCTGCTATC",
            "GGGAAAACGCACCTGCTATCG",
            "GGAAAACGCACCTGCTATCGA",
            "GAAAACGCACCTGCTATCGAT",
            "AAAACGCACCTGCTATCGATA",
            "AAACGCACCTGCTATCGATAG",
            "AACGCACCTGCTATCGATAGT",
            "ACGCACCTGCTATCGATAGTT",
            "CGCACCTGCTATCGATAGTTT",
            "GCACCTGCTATCGATAGTTTA",
            "CACCTGCTATCGATAGTTTAA",
            "ACCTGCTATCGATAGTTTAAC",
            "CCTGCTATCGATAGTTTAACC",
            "CTGCTATCGATAGTTTAACCT",
            "TGCTATCGATAGTTTAACCTT",
            "GCTATCGATAGTTTAACCTTT",
            "CTATCGATAGTTTAACCTTTC",
            "TATCGATAGTTTAACCTTTCT",
            "ATCGATAGTTTAACCTTTCTC",
            "TCGATAGTTTAACCTTTCTCT",
            "CGATAGTTTAACCTTTCTCTT",
            "GATAGTTTAACCTTTCTCTTT",
            "ATAGTTTAACCTTTCTCTTTG",
            "TAGTTTAACCTTTCTCTTTGA",
            "AGTTTAACCTTTCTCTTTGAG",
            "GTTTAACCTTTCTCTTTGAGG",
            "TTTAACCTTTCTCTTTGAGGT",
            "TTAACCTTTCTCTTTGAGGTT",
            "TAACCTTTCTCTTTGAGGTTT",
            "AACCTTTCTCTTTGAGGTTTT",
            "ACCTTTCTCTTTGAGGTTTTT",
            "CCTTTCTCTTTGAGGTTTTTG",
            "CTTTCTCTTTGAGGTTTTTGG",
            "TTTCTCTTTGAGGTTTTTGGA",
            "TTCTCTTTGAGGTTTTTGGAT",
            "TCTCTTTGAGGTTTTTGGATG",
            "CTCTTTGAGGTTTTTGGATGA",
            "TCTTTGAGGTTTTTGGATGAA",
            "CTTTGAGGTTTTTGGATGAAA",
            "TTTGAGGTTTTTGGATGAAAC",
            "TTGAGGTTTTTGGATGAAACG",
            "TGAGGTTTTTGGATGAAACGC",
            "GAGGTTTTTGGATGAAACGCA",
            "AGGTTTTTGGATGAAACGCAC",
            "GGTTTTTGGATGAAACGCACC",
            "GTTTTTGGATGAAACGCACCC",
            "TTTTTGGATGAAACGCACCCG",
            "TTTTGGATGAAACGCACCCGC",
            "TTTGGATGAAACGCACCCGCT",
            "TTGGATGAAACGCACCCGCTA",
            "TGGATGAAACGCACCCGCTAT",
            "GGATGAAACGCACCCGCTATC",
            "GATGAAACGCACCCGCTATCT",
            "ATGAAACGCACCCGCTATCTA",
            "TGAAACGCACCCGCTATCTAT",
            "GAAACGCACCCGCTATCTATC",
            "AAACGCACCCGCTATCTATCG",
            "AACGCACCCGCTATCTATCGT",
            "ACGCACCCGCTATCTATCGTT",
            "CGCACCCGCTATCTATCGTTA",
            "GCACCCGCTATCTATCGTTAC",
            "CACCCGCTATCTATCGTTACA",
            "ACCCGCTATCTATCGTTACAC",
            "CCCGCTATCTATCGTTACACT",
            "CCGCTATCTATCGTTACACTT",
            "CGCTATCTATCGTTACACTTC",
            "GCTATCTATCGTTACACTTCT",
            "CTATCTATCGTTACACTTCTC",
            "TATCTATCGTTACACTTCTCT",
            "ATCTATCGTTACACTTCTCTA",
            "TCTATCGTTACACTTCTCTAT",
            "CTATCGTTACACTTCTCTATC",
            "TATCGTTACACTTCTCTATCG",
            "ATCGTTACACTTCTCTATCGC",
            "TCGTTACACTTCTCTATCGCA",
            "CGTTACACTTCTCTATCGCAT",
            "GTTACACTTCTCTATCGCATT",
            "TTACACTTCTCTATCGCATTA",
            "TACACTTCTCTATCGCATTAA",
            "ACACTTCTCTATCGCATTAAC",
            "CACTTCTCTATCGCATTAACT",
            "ACTTCTCTATCGCATTAACTT",
            "CTTCTCTATCGCATTAACTTG",
            "TTCTCTATCGCATTAACTTGT",
            "TCTCTATCGCATTAACTTGTC",
            "CTCTATCGCATTAACTTGTCT",
            "TCTATCGCATTAACTTGTCTT",
            "CTATCGCATTAACTTGTCTTT",
            "TATCGCATTAACTTGTCTTTT",
            "ATCGCATTAACTTGTCTTTTG",
            "TCGCATTAACTTGTCTTTTGT",
            "CGCATTAACTTGTCTTTTGTA",
            "GCATTAACTTGTCTTTTGTAA",
            "CATTAACTTGTCTTTTGTAAC",
            "ATTAACTTGTCTTTTGTAACT",
            "TTAACTTGTCTTTTGTAACTT",
            "TAACTTGTCTTTTGTAACTTT",
            "AACTTGTCTTTTGTAACTTTG",
            "ACTTGTCTTTTGTAACTTTGC",
            "CTTGTCTTTTGTAACTTTGCT",
            "TTGTCTTTTGTAACTTTGCTC",
            "TGTCTTTTGTAACTTTGCTCG",
            "GTCTTTTGTAACTTTGCTCGT",
            "TCTTTTGTAACTTTGCTCGTT",
            "CTTTTGTAACTTTGCTCGTTT",
            "TTTTGTAACTTTGCTCGTTTT",
            "TTTGTAACTTTGCTCGTTTTG",
            "TTGTAACTTTGCTCGTTTTGT",
            "TGTAACTTTGCTCGTTTTGTT",
            "GTAACTTTGCTCGTTTTGTTA",
            "TAACTTTGCTCGTTTTGTTAG",
            "AACTTTGCTCGTTTTGTTAGG",
            "ACTTTGCTCGTTTTGTTAGGA",
            "CTTTGCTCGTTTTGTTAGGAA",
            "TTTGCTCGTTTTGTTAGGAAA",
            "TTGCTCGTTTTGTTAGGAAAG",
            "TGCTCGTTTTGTTAGGAAAGA",
            "GCTCGTTTTGTTAGGAAAGAT",
            "CTCGTTTTGTTAGGAAAGATA",
            "TCGTTTTGTTAGGAAAGATAG",
            "CGTTTTGTTAGGAAAGATAGA",
            "GTTTTGTTAGGAAAGATAGAT",
            "TTTTGTTAGGAAAGATAGATA",
            "TTTGTTAGGAAAGATAGATAC",
            "TTGTTAGGAAAGATAGATACA",
            "TGTTAGGAAAGATAGATACAC",
            "GTTAGGAAAGATAGATACACG",
            "TTAGGAAAGATAGATACACGA",
            "TAGGAAAGATAGATACACGAA",
            "AGGAAAGATAGATACACGAAG",
            "GGAAAGATAGATACACGAAGG",
            "GAAAGATAGATACACGAAGGA",
            "AAAGATAGATACACGAAGGAA",
            "AAGATAGATACACGAAGGAAG",
            "AGATAGATACACGAAGGAAGC",
            "GATAGATACACGAAGGAAGCA",
            "ATAGATACACGAAGGAAGCAC",
            "TAGATACACGAAGGAAGCACT",
            "AGATACACGAAGGAAGCACTC",
            "GATACACGAAGGAAGCACTCG",
            "ATACACGAAGGAAGCACTCGC",
            "TACACGAAGGAAGCACTCGCT",
            "ACACGAAGGAAGCACTCGCTC",
            "CACGAAGGAAGCACTCGCTCG",
            "ACGAAGGAAGCACTCGCTCGT",
            "CGAAGGAAGCACTCGCTCGTT",
            "GAAGGAAGCACTCGCTCGTTA",
            "AAGGAAGCACTCGCTCGTTAG",
            "AGGAAGCACTCGCTCGTTAGG",
            "GGAAGCACTCGCTCGTTAGGA",
            "GAAGCACTCGCTCGTTAGGAA",
            "AAGCACTCGCTCGTTAGGAAA",
            "AGCACTCGCTCGTTAGGAAAG",
            "GCACTCGCTCGTTAGGAAAGA",
            "CACTCGCTCGTTAGGAAAGAT",
            "ACTCGCTCGTTAGGAAAGATA",
            "CTCGCTCGTTAGGAAAGATAC",
            "TCGCTCGTTAGGAAAGATACA",
            "CGCTCGTTAGGAAAGATACAC",
            "GCTCGTTAGGAAAGATACACA",
            "CTCGTTAGGAAAGATACACAC",
            "TCGTTAGGAAAGATACACACG",
            "CGTTAGGAAAGATACACACGA",
            "GTTAGGAAAGATACACACGAG",
            "TTAGGAAAGATACACACGAGG",
            "TAGGAAAGATACACACGAGGA",
            "AGGAAAGATACACACGAGGAA",
            "GGAAAGATACACACGAGGAAG",
            "GAAAGATACACACGAGGAAGC",
            "AAAGATACACACGAGGAAGCA",
            "AAGATACACACGAGGAAGCAC",
            "AGATACACACGAGGAAGCACT",
            "GATACACACGAGGAAGCACTC",
            "ATACACACGAGGAAGCACTCG",
            "TACACACGAGGAAGCACTCGC",
            "ACACACGAGGAAGCACTCGCT",
            "CACACGAGGAAGCACTCGCTC",
            "ACACGAGGAAGCACTCGCTCG",
            "CACGAGGAAGCACTCGCTCGT",
            "ACGAGGAAGCACTCGCTCGTT",
            "CGAGGAAGCACTCGCTCGTTA",
            "GAGGAAGCACTCGCTCGTTAG",
            "AGGAAGCACTCGCTCGTTAGG",
            "GGAAGCACTCGCTCGTTAGGA",
            "GAAGCACTCGCTCGTTAGGAA",
            "AAGCACTCGCTCGTTAGGAAA",
            "AGCACTCGCTCGTTAGGAAAG",
            "GCACTCGCTCGTTAGGAAAGA",
            "CACTCGCTCGTTAGGAAAGAT",
            "ACTCGCTCGTTAGGAAAGATA",
            "CTCGCTCGTTAGGAAAGATAG",
            "TCGCTCGTTAGGAAAGATAGA",
            "CGCTCGTTAGGAAAGATAGAT",
            "GCTCGTTAGGAAAGATAGATA",
            "CTCGTTAGGAAAGATAGATAC",
            "TCGTTAGGAAAGATAGATACA",
            "CGTTAGGAAAGATAGATACAC",
            "GTTAGGAAAGATAGATACACG",
            "TTAGGAAAGATAGATACACGA",
            "TAGGAAAGATAGATACACGAG",
            "AGGAAAGATAGATACACGAGG",
            "GGAAAGATAGATACACGAGGA",
            "GAAAGATAGATACACGAGGAA",
            "AAAGATAGATACACGAGGAAG",
            "AAGATAGATACACGAGGAAGC",
            "AGATAGATACACGAGGAAGCC",
            "GATAGATACACGAGGAAGCCT",
            "ATAGATACACGAGGAAGCCTC",
            "TAGATACACGAGGAAGCCTCC",
            "AGATACACGAGGAAGCCTCCC",
            "GATACACGAGGAAGCCTCCCT",
            "ATACACGAGGAAGCCTCCCTC",
            "TACACGAGGAAGCCTCCCTCG",
            "ACACGAGGAAGCCTCCCTCGC",
            "CACGAGGAAGCCTCCCTCGCT",
            "ACGAGGAAGCCTCCCTCGCTA",
            "CGAGGAAGCCTCCCTCGCTAG",
            "GAGGAAGCCTCCCTCGCTAGA",
            "AGGAAGCCTCCCTCGCTAGAC",
            "GGAAGCCTCCCTCGCTAGACA",
            "GAAGCCTCCCTCGCTAGACAC",
            "AAGCCTCCCTCGCTAGACACT",
            "AGCCTCCCTCGCTAGACACTA",
            "GCCTCCCTCGCTAGACACTAG",
            "CCTCCCTCGCTAGACACTAGG",
            "CTCCCTCGCTAGACACTAGGA",
            "TCCCTCGCTAGACACTAGGAA",
            "CCCTCGCTAGACACTAGGAAG",
            "CCTCGCTAGACACTAGGAAGT",
            "CTCGCTAGACACTAGGAAGTC",
            "TCGCTAGACACTAGGAAGTCA",
            "CGCTAGACACTAGGAAGTCAT",
            "GCTAGACACTAGGAAGTCATC",
            "CTAGACACTAGGAAGTCATCG",
            "TAGACACTAGGAAGTCATCGC",
            "AGACACTAGGAAGTCATCGCT",
            "GACACTAGGAAGTCATCGCTC",
            "ACACTAGGAAGTCATCGCTCG",
            "CACTAGGAAGTCATCGCTCGT",
            "ACTAGGAAGTCATCGCTCGTT",
            "CTAGGAAGTCATCGCTCGTTA",
            "TAGGAAGTCATCGCTCGTTAG",
            "AGGAAGTCATCGCTCGTTAGG",
            "GGAAGTCATCGCTCGTTAGGA",
            "GAAGTCATCGCTCGTTAGGAA",
            "AAGTCATCGCTCGTTAGGAAA",
            "AGTCATCGCTCGTTAGGAAAG",
            "GTCATCGCTCGTTAGGAAAGA",
            "TCATCGCTCGTTAGGAAAGAC",
            "CATCGCTCGTTAGGAAAGACA",
            "ATCGCTCGTTAGGAAAGACAT",
            "TCGCTCGTTAGGAAAGACATA", 
            "CGCTCGTTAGGAAAGACATAA"));
        }
        p.close();
    }

    @Test
    public void shouldIterateKmersInMultipleSequences() throws IOException {
        FileProducer p = new FileProducer(1);
        final ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(() -> {p.run();});
        try(KMerIterator km = new LineKMerIteratorConsumer(21, "src/test/resources/fastaWith2Seq.fasta", p, 0, true)) {
            List<String>kmers = new ArrayList<>();
            while(km.hasNext()) {
                kmers.add(new String(km.next()));
            }
            assertThat(kmers, hasSize(800)); //both sequences are iterated independently, there is no k-mer spanning seq1 and seq2
        }
        p.close();
    }

    @Test
    public void shouldNotCreateKmers() throws IOException {
        FileProducer p = new FileProducer(1);
        final ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(() -> {p.run();});
        try(KMerIterator km = new LineKMerIteratorConsumer(21, "src/test/resources/fastaWithTwoShortSeq.fasta", p, 0, true)) {
            assertThat(km.hasNext(), equalTo(false));
        }
        p.close();
    }

    @Test
    public void checkComplement() throws IOException {
        FileProducer p = new FileProducer(1);
        final ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(() -> {p.run();});
        try(KMerIterator km = new LineKMerIteratorConsumer(21, "src/test/resources/fastaWith1Seq.fasta", p, 0, true)) {
            assertThat(km.hasNext(), equalTo(true));
            assertThat(new String(km.next()), equalTo("ACTCGTATGAACTTTGACTGG"));
            assertThat(new String(km.getReverseComplement()), equalTo("CCAGTCAAAGTTCATACGAGT"));
        }
        p.close();
    }

    @Test
    public void checkDiscardAmbiguousChars() throws IOException {
        FileProducer p = new FileProducer(1);
        final ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.submit(() -> {p.run();});
        try(KMerIterator km = new LineKMerIteratorConsumer(21, "src/test/resources/fastaWithAmbSeq.fasta", p, 0, true)) {
            int i = 0;
            while(km.hasNext()) {
                km.next();
                i++;
            }
            assertThat(i, equalTo((70*6) - 21 + 1 - 21));
        }
        p.close();
    }

    @Test
    @Ignore
    public void checkCoordinates() throws IOException {
        try(KMerIterator km = new InMemoryKMerIterator(8, "src/test/resources/fastaWithMultipleAmbSeq.fasta", true)) {
            int i = 0;
            int[] ambigPos = new int[]{0, 1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
            while(km.hasNext()) {
                km.next();
                var c = km.getCoordinates();
                if (i < 13) {
                    assertThat(c.getRecordIndexInFile(), equalTo(0));
                    assertThat(c.getSequenceIndexInRecord(), equalTo(i));
                    assertThat(c.getSequenceIndexInRecordIncludingAmbiguous(), equalTo(ambigPos[i]));
                    
                } else {
                    assertThat(c.getRecordIndexInFile(), equalTo(1));
                    assertThat(c.getSequenceIndexInRecord(), equalTo(i-13));
                    assertThat(c.getSequenceIndexInRecordIncludingAmbiguous(), equalTo(ambigPos[i-13]));
                }                
                i++;
            }
            
        }
    }  
}
