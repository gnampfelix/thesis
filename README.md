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
- `doc` (todo): The latex source code for the final report.
- `fmhdist`: The source code for a Java implementation of the concepts. Please
  see below for more details.
- `fmhdist-benchmark`: A benchmarking project, mainly used as a utility to
  benchmark the runtime of different hash functions and different approaches to
  insert-and-sort/insert-sorted into sets.
- `misc`: Skripts and utilitie used for the analysis.

## `fmhdist` utility
The `fmhdist` source code is a Maven project based on the
[SplitsTree6](https://github.com/husonlab/splitstree6) and
[jloda](https://github.com/husonlab/jloda3) repositories. You need to install
those artifacts in your local Maven repository. After that, navigate to the
`fmhdist` directory and execute `mvn clean install`. You might need to skip the
tests if you don't have internet connection or NCBI is unavailable.

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

## `misc` scripts
**todo**
