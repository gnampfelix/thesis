package org.husonlab.fmhdist.cmd;

import jloda.fx.util.ArgsOptions;
import jloda.fx.util.ProgramExecutorService;
import jloda.util.UsageException;

public class Main {
    private final static String CREATE_DB_COMMAND = "db";
    private final static String COMPARE_SKETCH_COMMAND = "dist";
    private final static String COMPARE_REF_SKETCH_COMMAND = "ref_dist";
    private final static String SKETCH_COMMAND = "sketch";

    public static void main(String[] args) throws UsageException {
        final ArgsOptions options = new ArgsOptions(args, Main.class,
                "Tool to calculate and compare FracMinHash sketches for sequence files");

        final String command = options.getCommand(
                new ArgsOptions.Command(
                        CREATE_DB_COMMAND,
                        "Create a new reference database from the given NCBI accession codes."),
                new ArgsOptions.Command(
                        COMPARE_SKETCH_COMMAND,
                        "Calculate the pairwise distances for all query sketches"),
                new ArgsOptions.Command(
                        COMPARE_REF_SKETCH_COMMAND,
                        "Find the closest database entries (in terms of distance) to given query sequences\n"
                                +
                                "and calculate all pairwise distances for those + the query sequences."),
                new ArgsOptions.Command(
                        SKETCH_COMMAND,
                        "Calculate the sketch for all given sequences and store them on the file system"));

        options.comment("Input/Output options");
        final String input = options.getOptionMandatory("-i", "input",
                String.format(
                        "New-line delimited list of\n" +
                                "(a - for %s command) NCBI accession codes\n" +
                                "(b - for other commands) sequence file paths or URLs to fasta files (gzip ok)"
                                +
                                ", a line consists of the mandatory path and (comma-separated) the optional label",
                        CREATE_DB_COMMAND, COMPARE_SKETCH_COMMAND),
                "");

        final String database = options.getOption("-db", "database",
                String.format(
                        "Path to reference database file (input for %s command, output for %s command)",
                        COMPARE_REF_SKETCH_COMMAND, CREATE_DB_COMMAND),
                "database.db");

        final String output = options.getOption("-o", "output",
                "Path to the output file in which the distances are written. Existing files will be overwritten. Format is nexus.",
                "out.nxs");

        options.comment("Algorithm parameters");
        final int kParameter = options.getOption("-k", "kmerSize", "Word size k", 21);
        final int sParameter = options.getOption("-s", "scalingFactor",
                "Scaling factor s. Hash values h are only part of the sketch if h <= H/s", 2000);
        final int randomSeed = options.getOption("-rs", "randomSeed", "Hashing random seed", 42);
        final double maxDistance = options.getOption("-md", "maxDistance",
                "The maximum distance that a query sequence should have to have to a reference sequence for the reference sequence to be included in the output",
                0.4);

        options.comment("Performance options");
        ProgramExecutorService
                .setNumberOfCoresToUse(options.getOption("-t", "threads", "Number of threads", 1));

        options.done();
        switch (command) {
            case CREATE_DB_COMMAND:
                DatabaseCreator dbCreator = new DatabaseCreator();
                dbCreator.run(input, database, kParameter, sParameter, randomSeed);
                break;
            case COMPARE_SKETCH_COMMAND:
                DistanceCalculator distanceCalculator = new DistanceCalculator();
                distanceCalculator.run(input, output);
                break;
            case COMPARE_REF_SKETCH_COMMAND:
                ReferenceDistanceCalculator refDistanceCalculator = new ReferenceDistanceCalculator();
                refDistanceCalculator.run(input, output, database, maxDistance);
                break;
            case SKETCH_COMMAND:
                SequenceSketcher sketcher = new SequenceSketcher();
                sketcher.run(input, output, kParameter, sParameter, randomSeed);
                break;
        }
    }
}
