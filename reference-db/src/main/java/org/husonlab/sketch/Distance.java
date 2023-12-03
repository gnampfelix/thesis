package org.husonlab.sketch;

/**
 * Distance between to sequences based on jaccard index estimation estimation
 */
public class Distance {
    public static double calculateJaccardIndex(long[] a, long[] b, int s) {
        int intersectionSize = getIntersectionSize(a, b);
        int unionSize = (a.length + b.length - intersectionSize);
        double jHat = (double)intersectionSize/(double)unionSize;
        double j = jHat / (1.0 - Math.pow(1.0 - (1.0 / (double)s), unionSize));
        return Math.min(1.0, j);
    }

    public static double calculateContainmentIndex(long[] a, long[] b, int s) {
        int intersectionSize = getIntersectionSize(a, b);
        double cHat = (double)intersectionSize/(double)a.length;
        double c = cHat/(1.0-Math.pow(1.0-(1.0/(double)s),a.length));
        return Math.min(1.0, c);
    }

    public static double jaccardToDistance(double j, int k) {
        return 1.0-Math.pow((2.0*j/(1.0+j)), 1.0/(double)k);
    }

    public static double jaccardToMashDistance(double j, int k) {
        return -1.0/(double)k * Math.log((2.0*j)/(1.0+j));
    }

    public static double containmentToDistance(double c, int k) {
        return 1.0 - Math.pow(c, 1.0/(double)k);
    }

    /**
     * 
     * @param a sorted array
     * @param b sorted array
     * @return
     */
    public static int getIntersectionSize(long[] a, long[] b) {
        int intersectionSize = 0;
        int i = 0;
        int j = 0;
        while(true) {
            if(a[i] < b[j]) {
                if(++i >= a.length) {
                    break;
                }
            } else if (a[i] > b[j]) {
                if(++j >= b.length) {
                    break;
                }
            } else {
                intersectionSize++;
                if(++i >= a.length || ++j >= b.length) {
                    break;
                }
            }
        }
        return intersectionSize;
    }
}
