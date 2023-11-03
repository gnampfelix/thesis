package org.husonlab;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.husonlab.db.ReferenceDatabase;
import org.husonlab.ncbi.Genome;
import org.husonlab.ncbi.NcbiApi;
import org.husonlab.ncbi.TaxonomyTree;
import org.husonlab.sketch.GenomeSketch;

import jloda.fx.util.ArgsOptions;
import jloda.fx.util.ProgramExecutorService;
import jloda.util.FileLineIterator;
import jloda.util.Single;
import jloda.util.UsageException;

public class Application {
    private static NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");

    public static void main(String[] args) throws UsageException {
        final ArgsOptions options = new ArgsOptions(args, Application.class,
                "Generates a reference database based on NCBI viral genomes and FracMinHash");

        final String input = options.getOptionMandatory("-i", "input",
                "New-line delimited list of genome accession codes", "");
        final String output = options.getOption("-o", "output",
                "Path to the output database file. Existing files will be overwritten", "database.db");

        final int kParameter = options.getOption("-k", "kmerSize", "Word size k", 21);
        final int sParameter = options.getOption("-s", "scalingFactor",
                "Scaling factor s. Hash values h are only part of the sketch if h <= H/s", 100);
        final int randomSeed = options.getOption("-rs", "randomSeed", "Hashing random seed", 42);

        ProgramExecutorService.setNumberOfCoresToUse(options.getOption("-t", "threads", "Number of threads", 8));

        options.done();

        try {
            FileLineIterator it = new FileLineIterator(input);
            List<String> accessionCodes = it.stream().map(line -> line.replaceAll("\\s+", ""))
                    .collect(Collectors.toList());
            it.close();

            List<Genome> genomes = api.getGenomes(accessionCodes);
            TaxonomyTree tree = api.getTaxonomyTreeForGenomes(genomes);
            if ((new File(output)).exists()) {
                (new File(output)).delete();
            }

            Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
            final Single<IOException> exception = new Single<>();
            final ExecutorService executor = Executors
                    .newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());
            try {
                genomes.forEach(genome -> executor.submit(() -> {
                    if (exception.isNull()) {
                        try {
                            GenomeSketch sketch = GenomeSketch.sketch(genome, kParameter, sParameter, randomSeed);
                            sketches.add(sketch);

                        } catch (IOException ex) {
                            exception.setIfCurrentValueIsNull(ex);
                        }
                    }
                }));
            } finally {
                executor.shutdown();
                executor.awaitTermination(1000, TimeUnit.DAYS);
            }
            for (Genome g : genomes) {

            }
            ReferenceDatabase db = ReferenceDatabase.create(output);
            db.insertGenomes(genomes);
            db.insertTaxonomy(tree);
            db.insertSketches(sketches);
            db.close();
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
