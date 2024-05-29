# thesis

This repository contains the main sources for my master thesis. The thesis is
about applications of
[FracMinHash](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2) to
evaluate the phylogenetic context of _Phytophthora_ samples using [Phylogenetic
Outlines](https://academic.oup.com/gbe/article/13/9/evab213/6370152).

## Overview
The repository is structured as follows:
- `analysis`: Computational artifacts generated while performing experiments
  using FracMinHash on different datasets. Each subfolder corresponds to a
  dataset/experiment. Please refer to the individual folder for further details.
- `doc`: The latex source code for the final report.
- `fmhdist`: The source code for a Java implementation of the concepts. Please
  see below for more details.
- `fmhdist-benchmark`: A benchmarking project, mainly used as a utility to
  benchmark the runtime of different hash functions and different approaches to
  insert-and-sort/insert-sorted into sets.
- `misc`: Skripts and utilities used for the analysis.

## `fmhdist` utility
The `fmhdist` source code is a Maven project based on the
[SplitsTree6](https://github.com/husonlab/splitstree6) and
[jloda](https://github.com/husonlab/jloda3) repositories. You need to install
those artifacts in your local Maven repository. After that, navigate to the
`fmhdist` directory and execute `mvn clean install`. You might need to skip the
tests if you don't have internet connection or NCBI is unavailable.

### Example
In `example`, there are five genomic sequences with their corresponding NCBI
accessions. Those are the genomes of some bacteria and some fungi.

Assuming you have the path to the built JAR file stored in `$fhmdist` and are in
the example directory of this repository, you can sketch those sequences using

```bash
for accession in */*.fna; do 
  echo $(pwd)/$accession,$(basename -s .fna $accession) >> sequences.csv; 
done
java -jar $fmhdist sketch --input sequences.csv --output . -t 2
```

and then estimate the distances using
```bash
for sketch in *.sketch; do
  echo $(pwd)/$sketch,$(basename -s .sketch $sketch) >> sketches.csv;
done
java -jar $fmhdist dist --input sketches.csv --output distance.nxs
```

### Functionality
`fmhdist` is a command line utility that is able to

- calculate the sketches for a list of genomic sequences
- estimate pairwise phylogenetic distances for a list of sketches
- estimate pairwise phylogenetic distances for a list of sketches to all closely
  related genomes from a reference database
- create such a reference database based on NCBI accession codes

#### Sketch calculation
```bash
java -jar fmhdist.jar sketch --input <input.csv> --output <directory/to/sketches/> [OPTIONAL PARAMS]
```

The `input.csv` must specify the sequences that should be sketched. URLs and
gzipped FASTA files are accepted, one per line. If needed, you can specify an
optional filename after a comma that will be used as the sketch filename.

The optional parameters are:
- `-k`: The $k$-mer size, default 21
- `-s`: The scaling parameter `s`. Only the hashes of $k$-mers that are below
  $\frac{H}{s}$ will be included in the sketch. $H$ is the maximum hash value
  that the hash function can produce. Default: 2000
- `-rs`: Random seed for the hash function initialization. Default: 42
- `-hf`: The hash function to use. Hash Functions are provided by
  [Zero-Allocation-Hashing](https://github.com/OpenHFT/Zero-Allocation-Hashing). Default: farm
- `-c`: Store coordinates of sketch hashes. This creates an additional file per
  sketch that lists the sketch's hash position in the original genome file. Default: false
- `-t`: The number of threads to use. Default: 1 

#### Distance Estimation
```bash
java -jar fmhdist.jar dist --input <input.csv> --output <directory/to/distances> 
```

The `input.csv` must specify the paths to the sketches that should be used for
distance estimation, one per line. If needed, you can specify an optional taxon
name after a comma that will be used in the generated nexus files.

This command estimates the Jaccard index $J_{frac} = \frac{\hat{J}_{frac}}{1 -
(1 - s)^{|A \cup B|}}$ with $\hat{J}_{frac}(A, B) = \frac{|\mathbf{FRAC}_s(A)
\cap \mathbf{FRAC}_s(B)|}{|\mathbf{FRAC}_s(A) \cup \mathbf{FRAC}_s(B)|}$ and the
containment index $C_{frac}(A, B) = \frac{|\mathbf{FRAC}_S(A) \cap
\mathbf{FRAC}_s(B)|}{|\mathbf{FRAC}_S(A)| (1-(1-s)^{|A|})}$ based on [Deriving
confidence intervals for mutation rates cross a wide range of evolutionary
distances using FracMinHash](https://genome.cshlp.org/content/33/7/1061).

Using those similarities, the command produces three distance matrices:
- one that is based on the estimated jaccard index $J$ and the corresponding
  formula presented in the publication above
- one that is based on the estimated containment index $C$ and the corresponding
- formula presented in the publication above
- one that is based on the estimated jaccard index $J$ and the formula for the
  [Mash-Distance](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)

**Todo**: Directly calculate the outline here.

#### Reference Distance estimation
```bash
java -jar fmhdist.jar ref_dist --input <input.csv> --database <path/to/db.db> --output <directory/to/distances> [-md 0.4] 
```

The `input.csv` must specify the paths to the sketches that should be used for
distance estimation, one per line. If needed, you can specify an optional taxon
name after a comma that will be used in the generated nexus files.

The `path/to/db.db` must either be a path to a valid SQLite reference DB (see
below) or another `input.csv` that will be used as a set of reference sketches.

The advantage of using the SQLite reference database is the automatic
application of the NCBI taxon names to the sequences involved. 

This command will estimate the pairwise distances of all sketches given via
`--input` to the sketches in `--database` and will only include those sketches
from `--database` in the output nexus file to which the distance is below `-md`.
Other than that, the command is identical to `dist`.

#### Reference DB creation
```bash
java -jar fmhdist.jar db --input <input.csv> --output <path/to/db.db> [OPTIONAL PARAMS]
```

This command creates a new reference database by treating every line of
`input.csv` as a NCBI accession code. The corresponding genome is downloaded
using the NCBI API and sketched. If there is already a `path/to/db.db`, it will
be overwritten. 

If multiple accession codes point to the same taxon, only the first is used for
the DB.

Other than that, the `[OPTIONAL PARAMS]` of the `sketch` command apply.

#### Outline Visualization
```bash
java -jar fmhdist.jar outline --input <output.nxs> --output <outline.svg> [OPTIONAL PARAMS]
```

This command visualizes the distances as a phylogenetic outline and stores the result in an SVG file.

The `[OPTIONAL PARAMS]` are special for this one:
- `iw`: the width of the image
- `ih`: the height of the image
- `is`: the scale of the image
- `ix`: the x-offset of the outline in the image
- `iy`: the y-offset of the outline in the image
- `il`: path to a TSV file that specifies alternative labels for the given taxa

## `misc` scripts
The `misc` folder contains scripts that I've used to automate several aspects.
Not all of them are built for simple re-use, I have created them when I needed a
particular thing done. For full transparency, I've included them in this
repository. For the python scripts, several dependencies need to be installed.
The full list can be found in `requirements.txt`.

The following scripts are supposed to be re-used: 

- `compare_splits.py`: This script takes a list of file names, each pointing to
a plaintext splits file as exported by SplitsTree 6. For full details on this
script, see section 3.5 of the thesis.
- `coordinate_complexity_correlation.py`: This script analyses the
  KMerCoordinates as produced by fmhdist. The idea is to supply different
  coordinate files that were created with the same settings but differen hash
  seeds. Also, this script needs complexities as generated by macle for each
  sequence in the FASTA, concatenated into a single text file. Optionally, you
  can provide genome annotations in GTF format. The result is an optional
  scatter plot of window hash counts against window sequence complexity, the
  pearson correlation coefficient R and the p-value of a Mann-Whitney U test.
  You may want to look at `prepare_analysis.sh` to prepare the macle
  complexities, the sequences, the sketches and the coordinates. For more
  details, please see section 3.8 of the thesis.
- `get_empty_intersections.py`: This script takes a list of paths to Mash
  sketches in JSON format and checks for each pairwise comparison if the
  intersection is empty. For more details, please refer to section 3.7 of the
  thesis.
- `prepare_anlaysis.sh`: This shell script prepares the analysis for
  `coordinate_complexity_correlation.py`. This script assumes that the current
  working directory has child directories, each representing a genome (e.g. as
  obtained when downloading from NCBI). This script assumes the `macle`
  executable can be found in the PATH as well as `java` and `splitfasta`. The
  genome will be split into individual sequences and `macle` will be applied on
  each sequence. The results are concatenated such that there is one large list
  of complexities per window size.  Note that the script renames files and
  directories and attempts to prepare GTF annotation files. If no such files are
  found, warnings will be thrown that can be ignored. Supply the following
  parameters:
  - `-w`: comma separated list of window sizes (applied for macle sequence
    complexity calculation)
  - `-s`: comma separated list of scaling parameters (applied for `fmhdist`
    sketch calculation)
  - `-k`: comma separated list of $k$-mer sizes (applied for `fmhdist` sketch
    calculation)
  - `-f`: path to the `fmhdist` jar-file
  - `-r`: comma separated list of random seeds
  - `-t`: number of threads used for sketch calculation
- `visualize_coordinates_complexity.py`: Creates a plot for each sequence in a
  FASTA file that displays the median hash count, the average hash count, the
  min hash count and the max hash counts for each (non-overlapping) window as
  well as the window's sequence complexity as calculated by macle. As an
  example, please see Figure 4.8 of the thesis.

  The other scripts and artifacts are:
- `calculate_ani.sh`: script that I've used to calculate the ANI values for a
    list of hardcoded genomes
- `dict_lookup.awk`: script to replace values in one text file with values
  defined in a dictionary by using `awk`
-  `j_to_d.py`: script that I have used to play around with some jaccard index conversions
-  `phytophthora_genome_sizes`: list of some _Phytophthora_ genome size
-  `phytophthora_sketch_sizes_k21_s2000`: list of some sketch sizes for _Phytophthora_.
-  `plot_genome_sizes.py`/`plot_sketch_sizes.py`: scripts that I've used to get
   a visual idea of sizes in the early stages of this project
- `virus_genome_sizes`: Sizes of some viral genomes
- `visualize_coordinates.py`: script that I've used to get an idea of the hash
  distribution in the original genomes.  
  
