package org.husonlab.fmhdist.cmd;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

import org.jfree.svg.SVGGraphics2D;

import jloda.graph.Edge;
import jloda.graph.Node;
import jloda.graph.NodeArray;
import jloda.phylo.PhyloSplitsGraph;
import jloda.util.parse.NexusStreamParser;
import jloda.util.progress.ProgressSilent;
import splitstree6.algorithms.distances.distances2splits.NeighborNet;
import splitstree6.data.DistancesBlock;
import splitstree6.data.SplitsBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.NexusImporter;
import splitstree6.layout.splits.algorithms.PhylogeneticOutline;
import javafx.geometry.Point2D;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

public class OutlineVisualizer {
    public void run(String input, String output, int width, int height, int scale, int xOffset, int yOffset) {
        Logger logger = Logger.getLogger(OutlineVisualizer.class.getName());
        try {       
            logger.info("Reading distances...");
            FileReader reader = new FileReader(input);
            NexusStreamParser parser = new NexusStreamParser(reader);
            TaxaBlock taxa = new TaxaBlock();
            DistancesBlock distances = new DistancesBlock();
            NexusImporter.parse(parser, taxa, distances);

            logger.info("Calculating splits...");
            SplitsBlock splits = new SplitsBlock();
            NeighborNet nn = new NeighborNet();
            nn.compute(new ProgressSilent(), taxa, distances, splits);

            logger.info("Calculating outline...");
            PhyloSplitsGraph graph = new PhyloSplitsGraph();
            NodeArray<Point2D> nodes = new NodeArray<>(graph);
            BitSet usedSplits = new BitSet();
            ArrayList<ArrayList<Node>> loops = new ArrayList<>();

            PhylogeneticOutline.apply(new ProgressSilent(), true, taxa, splits, graph, nodes, usedSplits, loops, 0, 0);

            logger.info("Creating output...");
            BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            Graphics2D graphics = img.createGraphics();
                        
            graphics.setBackground(Color.red);
            graphics.setPaint(Color.white);
            graphics.fillRect(0, 0, width, height);
            
            graphics.setPaint(Color.black);
            for(Edge e : graph.edges()) {
                Node s = e.getSource();
                Node t = e.getTarget();
                
                Point2D p1 = nodes.get(s);
                Point2D p2 = nodes.get(t);

                int x1 = (int)Math.floor(p1.getX() * scale) + xOffset;
                int y1 = (int)Math.floor(p1.getY() * scale) + yOffset;
                int x2 = (int)Math.floor(p2.getX() * scale) + xOffset;
                int y2 = (int)Math.floor(p2.getY() * scale) + yOffset;

                graphics.drawLine(x1, y1, x2, y2);
            }

            for(Node n : graph.leaves()) {
                String label = n.getLabel();
                Point2D p = nodes.get(n);
                int x = (int)Math.floor(p.getX() * scale) + xOffset;
                int y = (int)Math.floor(p.getY() * scale) + yOffset;
                graphics.drawString(label, x, y);
            }
            
            File imageFile = new File(output);
            if (imageFile.exists()) {
                imageFile.delete();
            }
            imageFile.createNewFile();
            ImageIO.write(img, "jpg", imageFile);
            System.out.println();
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
