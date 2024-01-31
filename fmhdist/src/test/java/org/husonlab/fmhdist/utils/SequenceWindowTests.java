package org.husonlab.fmhdist.utils;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;

import java.util.ArrayList;
import java.util.List;

import org.husonlab.fmhdist.util.SequenceWindow;
import org.junit.Test;

public class SequenceWindowTests {
    @Test
    public void shouldCalculateUniqueKmers() {
        List<byte[]> kmers = new ArrayList<>();
        kmers.add("ATG".getBytes());
        kmers.add("CTT".getBytes());
        kmers.add("CTA".getBytes());
        kmers.add("CTT".getBytes());
        SequenceWindow window = new SequenceWindow(0, 100, kmers);
        assertThat(window.getUniqueKmersInWindow(), equalTo(2));
    }
}
