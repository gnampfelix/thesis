package org.husonlab.fmhdist.cmd;

import java.io.File;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.db.ReferenceDatabase;
import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.ncbi.NcbiApi;
import org.husonlab.fmhdist.ncbi.TaxonomyTree;
import org.husonlab.fmhdist.sketch.GenomeSketch;

import jloda.fx.util.ProgramExecutorService;
import jloda.util.FileLineIterator;
import jloda.util.Single;
import net.openhft.hashing.LongHashFunction;

public class DatabaseCreator {
    private static NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");

    public void run(
            String input,
            String database,
            int kParameter,
            int sParameter,
            LongHashFunction hashFunction,
            int randomSeed) {
        Logger logger = Logger.getLogger(DatabaseCreator.class.getName());

        if ((new File(database)).exists()) {
            (new File(database)).delete();
        }

        try {
            FileLineIterator it = new FileLineIterator(input);
            List<String> accessionCodes = it.stream().map(line -> line.replaceAll("\\s+", ""))
                    .collect(Collectors.toList());
            it.close();

            List<Genome> genomes = api.getGenomes(accessionCodes);
            TaxonomyTree tree = api.getTaxonomyTreeForGenomes(genomes);

            Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
            final Single<Throwable> exception = new Single<>();
            final ExecutorService executor = Executors
                    .newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());

            try {
                genomes.forEach(genome -> executor.submit(() -> {
                    if (exception.isNull()) {
                        int retries = 5;
                        // Sometimes, the connection to NCBI breaks - this is a quick workaround
                        while (retries-- > 0) {
                            try {
                                GenomeSketch sketch = GenomeSketch.sketch(genome, kParameter, sParameter, hashFunction, randomSeed, false);
                                sketches.add(sketch);
                                break;
                            } catch (Exception ex) {
                                logger.warning(ex.getMessage());
                            } catch (Throwable e) {
                                // Somethings wrong here - no way to recover.
                                logger.severe(e.getMessage());
                                exception.setIfCurrentValueIsNull(e);
                                break;
                            }
                        }
                    }
                }));
            } finally {
                executor.shutdown();
                executor.awaitTermination(1000, TimeUnit.DAYS);
            }

            ReferenceDatabase db = ReferenceDatabase.create(database);
            db.insertTaxonomy(tree);
            db.insertSketches(sketches);
            db.insertInfo(kParameter, sParameter, randomSeed);
            db.close();
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
