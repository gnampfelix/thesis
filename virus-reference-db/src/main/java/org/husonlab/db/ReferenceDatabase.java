package org.husonlab.db;

import java.io.Closeable;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.List;

import org.husonlab.ncbi.Genome;
import org.sqlite.SQLiteConfig;

public class ReferenceDatabase implements Closeable{
    private Connection connection;
    
    public static ReferenceDatabase create(String path) throws SQLException {
        SQLiteConfig config = new SQLiteConfig();
        ReferenceDatabase result = new ReferenceDatabase();
        result.connection = config.createConnection("jdbc:sqlite:"+path);
        result.connection.createStatement().execute("CREATE TABLE bloom_filters (taxon_id INTEGER PRIMARY KEY, bloom_filter TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE tree (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE taxa (taxon_id INTEGER PRIMARY KEY, taxon_name TEXT, taxon_display_name TEXT, parent_id INTEGER REFERENCES taxa(taxon_id)) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE genomes (taxon_id INTEGER PRIMARY KEY, genome_accession TEXT NOT NULL, genome_size INTEGER, fasta_url TEXT) WITHOUT ROWID;");
        result.connection.createStatement().execute("CREATE TABLE mash_sketches (taxon_id INTEGER PRIMARY KEY, mash_sketch TEXT NOT NULL) WITHOUT ROWID;");
        return result;
    }

    private ReferenceDatabase() {}

    public void insertGenomes(List<Genome> genomes) throws SQLException {
        PreparedStatement s = this.connection.prepareStatement("INSERT INTO genomes (taxon_id, genome_accession, genome_size, fasta_url) VALUES (?, ?, ?, ?);");
        for (Genome g : genomes) {
            s.setInt(1, g.getTaxonId());
            s.setString(2, g.getAccession());
            s.setInt(3, g.getGenomeSize());
            s.setString(4, g.getFastaUrl());
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