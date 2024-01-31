package org.husonlab.fmhdist.cmd;

import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.husonlab.fmhdist.util.KMerCoordinates;
import org.husonlab.fmhdist.util.SequenceWindow;

public class CoordinatesAnalyzer {
    public void run(
        String input,
        String output,
        int windowSize
    ) {
        Logger logger = Logger.getLogger(CoordinatesAnalyzer.class.getName());
        logger.info("Reading coordinates file");
        List<SequenceWindow> windows = new ArrayList<>();

        try (Stream<String> stream = Files.lines(Paths.get(input))) {
            List<KMerCoordinates> coords = stream
                .map(KMerCoordinates::fromString)
                .collect(Collectors.toList());
            
            for(int i = 0; i < coords.size(); i++) {
                List<byte[]> currentKmers = new ArrayList<>();
                KMerCoordinates current = coords.get(i);
                currentKmers.add(current.getKmer());
                int j = i + 1;
                while (j < coords.size() && coords.get(j).getSequenceIndexInFileIncludingAmbiguous() < current.getSequenceIndexInFileIncludingAmbiguous() + windowSize) {
                    currentKmers.add(coords.get(j).getKmer());
                    j++;
                }
                windows.add(new SequenceWindow(current.getSequenceIndexInFileIncludingAmbiguous(), windowSize, currentKmers));
            }
            FileWriter writer = new FileWriter(Paths.get(output).toFile());
            for (SequenceWindow w : windows) {
                writer.write(w.toString());
                writer.write("\n");
            }
            writer.close();
        // TODO: write to file
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
