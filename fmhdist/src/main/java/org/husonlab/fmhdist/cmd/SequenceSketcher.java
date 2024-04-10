package org.husonlab.fmhdist.cmd;

import java.io.FileWriter;
import java.nio.file.Paths;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.sketch.GenomeSketch;
import org.husonlab.fmhdist.util.FileProducer;
import org.husonlab.fmhdist.util.KMerCoordinates;

import com.google.common.hash.HashFunction;

import jloda.fx.util.ProgramExecutorService;
import jloda.thirdparty.HexUtils;
import jloda.util.FileLineIterator;
import jloda.util.Single;
import net.openhft.hashing.LongHashFunction;

public class SequenceSketcher {

    public void run(
        String input,
        String output,
        int kParameter,
        int sParameter,
        LongHashFunction hashFunction,
        int randomSeed,
        boolean saveCoordinates
    ) {
        Logger logger = Logger.getLogger(SequenceSketcher.class.getName());
        try {
            logger.info("Parsing input file...");
            FileLineIterator it = new FileLineIterator(input);
            List<Genome> sequencePaths = it.stream()
                .map(line -> {
                    // Use second column as optional label for the genome
                    String[] comp = line.replaceAll("\\s+", "").split(",");
                    if(comp.length > 1) {
                        return new Genome(comp[1], comp[0]);
                    }
                    return new Genome(Paths.get(line).getFileName().toString(), line);
                })
                .collect(Collectors.toList());
            it.close();


            Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
            final Single<Throwable> exception = new Single<>();
            final ExecutorService executor = Executors
                    .newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());

            FileProducer producer = new FileProducer(ProgramExecutorService.getNumberOfCoresToUse());
            final ExecutorService ioExecutor = Executors.newFixedThreadPool(1);
            ioExecutor.submit(() -> {producer.run();});

            logger.info("Sketching sequences...");
            try {
                sequencePaths.forEach(genome -> executor.submit(() -> {
                    int threadIndex = Integer.parseInt(Thread.currentThread().getName().split("-")[3])-1;
                    if (exception.isNull()) {
                        try {
                            GenomeSketch sketch = GenomeSketch.sketch(producer, threadIndex, genome, kParameter, sParameter, hashFunction, randomSeed, saveCoordinates);
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

            producer.close();
            ioExecutor.shutdown();
            ioExecutor.awaitTermination(10, TimeUnit.SECONDS);

            logger.info("Saving sketches...");
            for(GenomeSketch sketch : sketches) {
                logger.fine(String.format("Saving %s...", sketch.getGenome().getOrganismName()));
                FileWriter writer = new FileWriter(Paths.get(output, String.format("%s.sketch", sketch.getGenome().getOrganismName())).toFile());
                writer.write(HexUtils.encodeHexString(sketch.getSketch().getBytes()));
                writer.close();

                if (saveCoordinates) {
                    writer = new FileWriter(Paths.get(output, String.format("%s.sketch.coordinates", sketch.getGenome().getOrganismName())).toFile());
                    List<KMerCoordinates> coordinates = sketch.getSketch().getCoordinates();
                    for (KMerCoordinates coord : coordinates) {
                        writer.write(coord.toString());
                        writer.write("\n");
                    }
                    writer.close();
                }
            } 
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
    
}
