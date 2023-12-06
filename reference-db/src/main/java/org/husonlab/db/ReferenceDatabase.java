package org.husonlab.db;

import java.io.Closeable;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.logging.Logger;

import org.husonlab.ncbi.Taxon;
import org.husonlab.ncbi.TaxonomyTree;
import org.husonlab.sketch.GenomeSketch;
import org.sqlite.SQLiteConfig;

import jloda.thirdparty.HexUtils;

public class ReferenceDatabase implements Closeable{
    private Connection connection;
    private Logger logger;
    
    public static ReferenceDatabase create(String path) throws SQLException {
        SQLiteConfig config = new SQLiteConfig();
        ReferenceDatabase result = new ReferenceDatabase();
        result.connection = config.createConnection("jdbc:sqlite:"+path);
        result.connection.createStatement().execute("CREATE TABLE bloom_filters (taxon_id INTEGER PRIMARY KEY, bloom_filter TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE tree (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE taxa (taxon_id INTEGER PRIMARY KEY, taxon_name TEXT, taxon_display_name TEXT, parent_id INTEGER REFERENCES taxa(taxon_id)) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE genomes (taxon_id INTEGER PRIMARY KEY, genome_accession TEXT NOT NULL, genome_size INTEGER, fasta_url TEXT) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE frac_min_hash_sketches (taxon_id INTEGER PRIMARY KEY, frac_min_hash_sketch TEXT NOT NULL) WITHOUT ROWID;");
        return result;
    }

    public static ReferenceDatabase open(String path) throws SQLException {
        SQLiteConfig config = new SQLiteConfig();
        ReferenceDatabase result = new ReferenceDatabase();
        result.connection = config.createConnection("jdbc:sqlite:"+path);
        return result;
    }

    private ReferenceDatabase() {
        this.logger = Logger.getLogger(ReferenceDatabase.class.getName());
    }

    public void insertSketches(Collection<GenomeSketch> sketches) throws SQLException {
        this.logger.info("Inserting genomes and sketches in DB...");
        PreparedStatement genome_statement = this.connection.prepareStatement("INSERT INTO genomes (taxon_id, genome_accession, genome_size, fasta_url) VALUES (?, ?, ?, ?);");
        PreparedStatement sketch_statement = this.connection.prepareStatement("INSERT INTO frac_min_hash_sketches (taxon_id, frac_min_hash_sketch) VALUES (?, ?);");
        PreparedStatement check_unique_statement = this.connection.prepareStatement("SELECT 1 FROM genomes WHERE taxon_id = ?");
        for (GenomeSketch g : sketches) {
            // Is Taxon already known?
            check_unique_statement.setInt(1, g.getGenome().getTaxonId());
            ResultSet rs = check_unique_statement.executeQuery();
            if (rs.next()) {
                this.logger.warning(String.format("multiple accession codes for taxon %d found", g.getGenome().getTaxonId()));
                continue;
            }

            genome_statement.setInt(1, g.getGenome().getTaxonId());
            genome_statement.setString(2, g.getGenome().getAccession());
            genome_statement.setLong(3, g.getGenome().getGenomeSize());
            genome_statement.setString(4, g.getGenome().getFastaUrl());
            genome_statement.executeUpdate();

            sketch_statement.setInt(1, g.getGenome().getTaxonId());
            sketch_statement.setString(2, HexUtils.encodeHexString(g.getSketch().getBytes()));
            sketch_statement.executeUpdate();
        }
        this.logger.info("Finished inserting genome sketches in DB!");
    }

    public void insertTaxonomy(TaxonomyTree taxonomy) throws SQLException {
        this.logger.info("Inserting taxa in DB...");
        PreparedStatement s = this.connection.prepareStatement("INSERT INTO taxa (taxon_id, taxon_name, taxon_display_name, parent_id) VALUES (?, ?, ?, ?)");
        for (Taxon t : taxonomy.getTaxa().values()) {
            s.setInt(1, t.getTaxonId());
            s.setString(2, t.getOrganismName());
            s.setString(3, t.getOrganismName());
            int parent = 0;
            if (!taxonomy.getRoot().equals(t)) {
                for (Taxon p : taxonomy.getTree().predecessors(t)) { // there is only one
                    parent = p.getTaxonId();
                }
            }
            s.setInt(4, parent);
            s.executeUpdate();
        }
        this.logger.info("Finished inserting taxa in DB!");
    }

    public void insertInfo(int kSize, int sParam, int seed) throws SQLException{
        this.logger.info("Inserting sketch creation info...");
        PreparedStatement s = this.connection.prepareStatement("INSERT INTO info (key, value) VALUES (?, ?);");
        int[] values = new int[]{kSize, sParam, seed};
        String[] keys = new String[]{"sketch_k", "sketch_s", "sketch_seed"};
        for (int i=0; i < 3; i++) {
            s.setInt(1, values[i]);
            s.setString(2, keys[i]);
            s.executeUpdate();
        }
    }

    @Override
    public void close() throws IOException {
        try {
            this.connection.close();
        } catch (SQLException e) {
            throw new IOException(e);
        }
    }    
}