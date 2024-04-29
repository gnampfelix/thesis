package org.husonlab.fmhdist.cmd;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.db.ReferenceDatabase;
import org.husonlab.fmhdist.sketch.Distance;
import org.husonlab.fmhdist.sketch.FracMinHashSketch;
import org.husonlab.fmhdist.sketch.GenomeSketch;
import org.husonlab.fmhdist.sketch.IncompatibleParameterException;
import org.husonlab.fmhdist.util.HashFunctionParser;

import jloda.thirdparty.HexUtils;
import jloda.util.FileLineIterator;
import splitstree6.data.DistancesBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.writers.distances.NexusWriter;

public class ReferenceDistanceCalculator {

    private FracMinHashSketch lineToSketch(String line){        
        try {
            String[] splitLine = line.replaceAll("\\s+", "").split(",");
            byte[] content = Files.readAllBytes(Paths.get(splitLine[0]));
            FracMinHashSketch result = FracMinHashSketch.parse(HexUtils.decodeHexString(new String(content)));
            if(splitLine.length > 1) {
                result.setName(splitLine[1]);
            } else {
                result.setName(splitLine[0]);
            }
            return result; 
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

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
            Map<String, Integer> info = db.getNumericalInfo();
            String hashFunctionName = db.getUsedHashFunction();
            Collection<GenomeSketch> refSketches = db.getSketches();
            db.close();

            if (!info.containsKey("sketch_k") ||
                    !info.containsKey("sketch_s") ||
                    !info.containsKey("sketch_seed") ||
                    hashFunctionName.equals("")) {
                throw new IncompatibleParameterException(
                        "reference db does not provide all sketching parameters (s, k, seed, hash function)");
            }

            if (!HashFunctionParser.getSupportedFunctions().contains(hashFunctionName)) {
                throw new IncompatibleParameterException(String.format("hash function '%s' used in database is not supported", hashFunctionName));
            }
            int kParameter = info.get("sketch_k");
            int sParameter = info.get("sketch_s");
            int randomSeed = info.get("sketch_seed");              

            long hashedMagicNumber = FracMinHashSketch.getHashedMagicNumber(HashFunctionParser.createHashFunction(hashFunctionName, randomSeed));

            logger.info("Reading queries list...");
            FileLineIterator it = new FileLineIterator(input);
            List<FracMinHashSketch> sketches = it
                .stream()
                .map(this::lineToSketch)
                .collect(Collectors.toList());
            it.close();

            logger.info("Finding closest reference genomes...");
            Set<FracMinHashSketch> resultSketchSet = new HashSet<>();
            for (FracMinHashSketch querySketch : sketches) {
                
                if (
                    sParameter != querySketch.getSParam() || 
                    kParameter != querySketch.getKSize() || 
                    randomSeed != querySketch.getSeed() || 
                    hashedMagicNumber != querySketch.getHashedMagicNumber()
                ) {
                    logger.severe("sketches have incompatible sketching parameters");
                    return;
                }

                for (GenomeSketch refSketch : refSketches) {
                    double jaccard = Distance.calculateJaccardIndex(querySketch.getValues(),
                            refSketch.getSketch().getValues(), sParameter);
                    double distance = Distance.jaccardToDistance(jaccard, kParameter);
                    if (distance <= maxDistance) {
                        resultSketchSet.add(refSketch.getSketch());
                    }
                }
                resultSketchSet.add(querySketch);
            }

            logger.info("Calculating pairwise distances...");
            List<FracMinHashSketch> resultSketchesList = new ArrayList<>(resultSketchSet);

            DistancesBlock distances_jaccard = new DistancesBlock();
            distances_jaccard.setNtax(resultSketchSet.size());

            DistancesBlock distances_mash = new DistancesBlock();
            distances_mash.setNtax(resultSketchSet.size());

            DistancesBlock distances_containment = new DistancesBlock();
            distances_containment.setNtax(resultSketchSet.size());

            TaxaBlock taxa = new TaxaBlock();
            for (int i = 0; i < resultSketchSet.size(); i++) {
                taxa.addTaxonByName(resultSketchesList.get(i).getName());
                for (int j = i; j < resultSketchSet.size(); j++) {
                    double jaccard = Distance.calculateJaccardIndex(
                        resultSketchesList.get(i).getValues(),
                        resultSketchesList.get(j).getValues(), 
                        sParameter
                    );
                    // Containment is not symmetrical
                    double containment_i = Distance.calculateContainmentIndex(
                        resultSketchesList.get(i).getValues(),
                        resultSketchesList.get(j).getValues(), 
                        sParameter
                    );
                    double containment_j = Distance.calculateContainmentIndex(
                        resultSketchesList.get(j).getValues(),
                        resultSketchesList.get(i).getValues(), 
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
