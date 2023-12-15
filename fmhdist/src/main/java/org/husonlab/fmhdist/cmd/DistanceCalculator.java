package org.husonlab.fmhdist.cmd;

import java.io.FileWriter;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.db.ReferenceDatabase;
import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.sketch.Distance;
import org.husonlab.fmhdist.sketch.GenomeSketch;
import org.husonlab.fmhdist.sketch.IncompatibleParameterException;

import jloda.fx.util.ProgramExecutorService;
import jloda.util.FileLineIterator;
import jloda.util.Single;
import jloda.util.progress.ProgressSilent;
import splitstree6.data.DistancesBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.writers.distances.NexusWriter;

public class DistanceCalculator {

    public void run(
            String input,
            String output,
            String database,
            double maxDistance
            ) {
        Logger logger = Logger.getLogger(DistanceCalculator.class.getName());
        try {
            logger.info("Loading reference DB!");
            ReferenceDatabase db = ReferenceDatabase.open(database);
            Map<String, Integer> info = db.getInfo();
            Collection<GenomeSketch> refSketches = db.getSketches();
            db.close();

            if (!info.containsKey("sketch_k") ||
                    !info.containsKey("sketch_s") ||
                    !info.containsKey("sketch_seed")) {
                throw new IncompatibleParameterException(
                        "reference db does not provide all sketching parameters (s, k, seed)");
            }
            int kParameter = info.get("sketch_k");
            int sParameter = info.get("sketch_s");
            int randomSeed = info.get("sketch_seed");

            logger.info("Reading queries list...");
            FileLineIterator it = new FileLineIterator(input);
            List<Genome> sequencePaths = it.stream().map(line -> line.replaceAll("\\s+", ""))
                    .map(line -> new Genome(Paths.get(line).getFileName().toString(), line))
                    .collect(Collectors.toList());
            it.close();

            Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
            final Single<Throwable> exception = new Single<>();
            final ExecutorService executor = Executors
                    .newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());

            logger.info("Sketching query sequences...");
            logger.info(String.format("Using reference DB parameters: s: %d, k: %d, seed: %d", sParameter, kParameter, randomSeed));
            try {
                sequencePaths.forEach(genome -> executor.submit(() -> {
                    if (exception.isNull()) {
                        try {
                            GenomeSketch sketch = GenomeSketch.sketch(genome, kParameter, sParameter, randomSeed, false,
                                    false, new ProgressSilent());
                            sketches.add(sketch);
                        } catch (Exception ex) {
                            logger.warning(ex.getMessage());
                        } catch (Throwable e) {
                            // Somethings wrong here - no way to recover.
                            logger.severe(e.getMessage());
                            exception.setIfCurrentValueIsNull(e);
                        }
                    }
                }));
            } finally {
                executor.shutdown();
                executor.awaitTermination(1000, TimeUnit.DAYS);
            }



            logger.info("Finding closest reference genomes...");
            Set<GenomeSketch> resultSketchSet = new HashSet<>();
            for (GenomeSketch querySketch : sketches) {
                for (GenomeSketch refSketch : refSketches) {
                    double jaccard = Distance.calculateJaccardIndex(querySketch.getSketch().getValues(),
                            refSketch.getSketch().getValues(), sParameter);
                    double distance = Distance.jaccardToDistance(jaccard, kParameter);
                    if (distance <= maxDistance) {
                        resultSketchSet.add(refSketch);
                    }
                }
                resultSketchSet.add(querySketch);
            }

            logger.info("Calculating pairwise distances...");
            List<GenomeSketch> resultSketchesList = new ArrayList<>(resultSketchSet);

            DistancesBlock distances_jaccard = new DistancesBlock();
            distances_jaccard.setNtax(resultSketchSet.size());

            DistancesBlock distances_mash = new DistancesBlock();
            distances_mash.setNtax(resultSketchSet.size());

            DistancesBlock distances_containment = new DistancesBlock();
            distances_containment.setNtax(resultSketchSet.size());

            TaxaBlock taxa = new TaxaBlock();
            for (int i = 0; i < resultSketchSet.size(); i++) {
                taxa.addTaxonByName(resultSketchesList.get(i).getGenome().getOrganismName());
                for (int j = i; j < resultSketchSet.size(); j++) {
                    double jaccard = Distance.calculateJaccardIndex(
                        resultSketchesList.get(i).getSketch().getValues(),
                        resultSketchesList.get(j).getSketch().getValues(), 
                        sParameter
                    );
                    // Containment is not symmetrical
                    double containment_i = Distance.calculateContainmentIndex(
                        resultSketchesList.get(i).getSketch().getValues(),
                        resultSketchesList.get(j).getSketch().getValues(), 
                        sParameter
                    );
                    double containment_j = Distance.calculateContainmentIndex(
                        resultSketchesList.get(j).getSketch().getValues(),
                        resultSketchesList.get(i).getSketch().getValues(), 
                        sParameter
                    );

                    
                    distances_jaccard.setBoth(i + 1, j + 1, Distance.jaccardToDistance(jaccard, kParameter)); // for some reason, the method is 1-based
                    distances_containment.set(i+1, j+1, Distance.containmentToDistance(containment_i, kParameter));
                    distances_containment.set(j+1, i+1, Distance.containmentToDistance(containment_j, kParameter));
                    distances_mash.setBoth(i+1, j+1, Distance.jaccardToMashDistance(jaccard, kParameter));
                }
            }

            logger.info("Exporting...");
            FileWriter outFile = new FileWriter(output, false);
            outFile.write("#nexus\n");
            NexusWriter writer = new NexusWriter();
            writer.write(outFile, taxa, distances_jaccard);
            outFile.close();

            outFile = new FileWriter(output + ".containment", false);
            outFile.write("#nexus\n");
            writer = new NexusWriter();
            writer.write(outFile, taxa, distances_containment);
            outFile.close();

            outFile = new FileWriter(output + ".mash", false);
            outFile.write("#nexus\n");
            writer = new NexusWriter();
            writer.write(outFile, taxa, distances_containment);
            outFile.close();

        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
