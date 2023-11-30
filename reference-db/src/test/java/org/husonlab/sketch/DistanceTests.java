package org.husonlab.sketch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.husonlab.util.KMerIterator;
import org.junit.Test;

import jloda.util.FileLineBytesIterator;
import jloda.util.progress.ProgressSilent;


public class DistanceTests {
    @Test
    public void testIntersectionSize() {
        long[] a = new long[]{1, 2, 3};
        long[] b = new long[]{2};
        assertThat(Distance.getIntersectionSize(a, b), equalTo(1));
    }

    @Test
    public void testJaccardIndex() {
        long[] a = new long[]{1, 2, 3};
        long[] b = new long[]{2};
        assertThat(Distance.calculateJaccardIndex(a, b, 1), equalTo(1.0/3.0));
    }

    @Test
    public void testJaccardIndexWithSelf() {
        long[] a = new long[]{1, 2, 3};
        assertThat(Distance.calculateJaccardIndex(a, a, 1), equalTo(1.0));
    }

    @Test
    public void testContainmentIndex() {
        long[] a = new long[]{1, 2, 3};
        long[] b = new long[]{2};
        assertThat(Distance.calculateContainmentIndex(b, a, 1), equalTo(1.0));
    }

    @Test
    public void testJaccrdDistanceWithSelf() {
        assertThat(Distance.jaccardToDistance(1, 21), equalTo(0.0));
    }

    @Test
    public void shouldCalculateDistancesCorrectly() throws IOException, IncompatibleParameterException {
        final String[] downloadUrls = new String[]{
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/295/175/GCA_012295175.1_ASM1229517v1/GCA_012295175.1_ASM1229517v1_genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/316/315/GCA_011316315.1_UAnd_PInf_RC1-10.1/GCA_011316315.1_UAnd_PInf_RC1-10.1_genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/300/685/GCA_024300685.1_ASM2430068v1/GCA_024300685.1_ASM2430068v1_genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/018/455/GCA_003018455.1_ASM301845v1/GCA_003018455.1_ASM301845v1_genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/722/995/GCA_025722995.1_ASM2572299v1/GCA_025722995.1_ASM2572299v1_genomic.fna.gz",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/333/115/GCA_000333115.2_PhyKer844_4v2/GCA_000333115.2_PhyKer844_4v2_genomic.fna.gz"
        };

        final String[] taxa = new String[]{
            "P. infestans",
            "P. infestans",
            "E. coli",
            "E. coli",
            "P. agathidic",
            "P. kernoviae"
        };

        long start = System.currentTimeMillis();
        final int k = 21;
        final int s = 1000;
        List<FracMinHashSketch> sketches = new ArrayList<>(); 
        for(String url : downloadUrls) {
            try(FileLineBytesIterator it = new FileLineBytesIterator(url)) {
                KMerIterator kmers = new KMerIterator(it, k);
                sketches.add(FracMinHashSketch.compute(
                    url,
                    kmers,
                    200000000,
                    true,
                    s,
                    42,
                    false,
                    false,
                    new ProgressSilent()));
            }
            System.out.println(String.format("%d seconds", (System.currentTimeMillis()-start)/1000));
        }

        // each entry consists of the three different distances
        List<List<List<Double>>> distances = new ArrayList<>();
        for (int i = 0; i < sketches.size(); i++) {
            List<List<Double>> currentRow = new ArrayList<>();
            for (int j = i; j < sketches.size(); j++) {
                List<Double> currentDistances = new ArrayList<>();
                double containment = Distance.calculateContainmentIndex(sketches.get(i).getValues(), sketches.get(j).getValues(), s);
                double jaccard = Distance.calculateJaccardIndex(sketches.get(i).getValues(), sketches.get(j).getValues(), s);
                currentDistances.add(Distance.containmentToDistance(containment, k));
                currentDistances.add(Distance.jaccardToDistance(jaccard, k));
                currentDistances.add(Distance.jaccardToMashDistance(jaccard, k));
                currentRow.add(currentDistances);
            }
            distances.add(currentRow);
        }
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < distances.size(); j++) {
                System.out.print(String.format("%s:\t", taxa[j]));
                for (int l = 0; l < distances.get(j).size(); l++) {
                    System.out.print(String.format("%f\t", distances.get(j).get(l).get(i)));
                }
                System.out.print("\n");
            }
            System.out.print("\n");
        }
        System.out.println(distances);
        System.out.println();
        
    }
}
