package org.husonlab.cmd;

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

import org.husonlab.db.ReferenceDatabase;
import org.husonlab.ncbi.Genome;
import org.husonlab.ncbi.NcbiApi;
import org.husonlab.sketch.Distance;
import org.husonlab.sketch.GenomeSketch;
import org.husonlab.sketch.IncompatibleParameterException;


import jloda.fx.util.ArgsOptions;
import jloda.fx.util.ProgramExecutorService;
import jloda.util.FileLineIterator;
import jloda.util.Single;
import jloda.util.UsageException;
import jloda.util.progress.ProgressSilent;
import splitstree6.data.DistancesBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.writers.distances.NexusWriter;
import splitstree6.io.writers.distances.PhylipWriter;

public class DistanceCalculator {
    
    private static NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");

    public static void main(String[] args) throws UsageException {
        final ArgsOptions options = new ArgsOptions(args, DistanceCalculator.class,
                "Calculates the distances from a set of query sequences to a reference database");

        final String input = options.getOptionMandatory("-i", "input",
                "New-line delimited list of sequence file paths or URLs to fasta files (gzip ok)", "");
        
        final String database = options.getOptionMandatory("-db", "database",
                "Path to reference database file", "database.db");

        final String output = options.getOption("-o", "output",
                "Path to the output file in which the distances are written. Existing files will be overwritten", "out.nex");

        final int kParameter = options.getOption("-k", "kmerSize", "Word size k", 21);
        final int sParameter = options.getOption("-s", "scalingFactor",
                "Scaling factor s. Hash values h are only part of the sketch if h <= H/s", 2000);
        final int randomSeed = options.getOption("-rs", "randomSeed", "Hashing random seed", 42);
        final double maxDistance = options.getOption("-md", "maxDistance", 
            "The maximum distance that a query sequence should have to have to a reference sequence for the reference sequence to be included in the output", 0.4);

        // TODO: is much slower when multithreading - maybe IO is limiting factor here?
        ProgramExecutorService.setNumberOfCoresToUse(options.getOption("-t", "threads", "Number of threads", 1));

        options.done();
        Logger logger = Logger.getLogger(DBCreation.class.getName());
        try {
            FileLineIterator it = new FileLineIterator(input);
            List<Genome> sequencePaths = it.stream().map(line -> line.replaceAll("\\s+", "")).map(line -> new Genome(Paths.get(line).getFileName().toString(), line))
                    .collect(Collectors.toList());
            it.close();

            Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
            final Single<Throwable> exception = new Single<>();
            final ExecutorService executor = Executors
                    .newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());

            logger.info("Sketching query sequences...");
            try {
                sequencePaths.forEach(genome -> executor.submit(() -> {
                    if (exception.isNull()) {
                        try {
                            GenomeSketch sketch = GenomeSketch.sketch(genome, kParameter, sParameter, randomSeed, false, false, new ProgressSilent());
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
            
            logger.info("Loading reference DB!");
            ReferenceDatabase db = ReferenceDatabase.open(database);
            Map<String, Integer> info = db.getInfo();
            Collection<GenomeSketch> refSketches = db.getSketches();
            db.close();
            if (
                !info.containsKey("sketch_k") || 
                !info.containsKey("sketch_s") || 
                !info.containsKey("sketch_seed") ||
                info.get("sketch_k") != kParameter ||
                info.get("sketch_s") != sParameter ||
                info.get("sketch_seed") != randomSeed
            ) {
                throw new IncompatibleParameterException("passed sketch params (s, k, seed) are not compatible to the sketch params of the passed reference DB!");
            }

            logger.info("Finding closest reference genomes...");
            Set<GenomeSketch> resultSketchSet = new HashSet<>();
            for (GenomeSketch querySketch : sketches) {
                for (GenomeSketch refSketch : refSketches) {
                    double jaccard = Distance.calculateJaccardIndex(querySketch.getSketch().getValues(), refSketch.getSketch().getValues(), sParameter);
                    double distance = Distance.jaccardToDistance(jaccard, kParameter);
                    if (distance <= maxDistance){
                        resultSketchSet.add(refSketch);
                    }
                }
                resultSketchSet.add(querySketch);
            }

            logger.info("Calculating pairwise distances...");
            List<GenomeSketch> resultSketchesList = new ArrayList<>(resultSketchSet);
            DistancesBlock distances = new DistancesBlock();
            distances.setNtax(resultSketchSet.size());
            TaxaBlock taxa = new TaxaBlock();
            for (int i = 0; i < resultSketchSet.size(); i++) {
                taxa.addTaxonByName(resultSketchesList.get(i).getGenome().getAccession().substring(4, resultSketchesList.get(i).getGenome().getAccession().length()));
                for (int j = i; j < resultSketchSet.size(); j++) {
                    double jaccard = Distance.calculateJaccardIndex(resultSketchesList.get(i).getSketch().getValues(), resultSketchesList.get(j).getSketch().getValues(), sParameter);
                    double distance = Distance.jaccardToDistance(jaccard, kParameter);
                    distances.setBoth(i+1, j+1, distance); //for some reason, the method is 1-based
                }
            }

            FileWriter outFile = new FileWriter(output, false);
            PhylipWriter writer = new PhylipWriter();
            writer.write(outFile, taxa, distances);
            outFile.close();
            
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}

