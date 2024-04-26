import argparse
import sys
import miniFasta as mf
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import gtfparse as gtf

gene_search_terms = [
    "Crinkler","CRN","PTHR33129","PTHR24362",
    "Cutinase","PF01083","IPR000675",
    "Elicitin","IPR002200","IPR036470","PF00964",
    "NLP","IPR008701","PF05630","Necrosis inducing","NPP1",
    "Phytotoxin","PF09461","IPR018570","PcF",
    "Cystatin","IPR000010","IPR027214","IPR018073","IPR020381","IPR002350","Kazal","PF00050","PF07648","Cathepsin","IPR013201","EPIC","protease inhibitor",
    "RXLR","IPR031825","PF16810",
    "Avirulence","Effector","Av"
]

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-w", "--window-size", help="The window size for density analysis", type=int, default=1000, required=False)
    p.add_argument("-a", "--annotations", help="The path to the annotations file in GTF format", type=str, required=False)
    p.add_argument("-sc", "--scaling", help="The scaling parameter s that was used to generate the sketch on which the coordinates are based", type=int, default=2000)
    p.add_argument("-c", "--coordinates", help="The coordinates file to analyze", required=True, nargs="+")
    p.add_argument("-s", "--sequence", help="The underlying sequence file in FASTA format", required=True)
    p.add_argument("-m", "--macle", help="The path to the macle file", required=True)
    p.add_argument("-o", "--overlap", help="Indicates if the window offset should be 1, otherwise <window-size>", action="store_true")
    p.add_argument("-dr", "--density-range", help="The range in which the density in a window is considered 'expected' is given by w/s +/- w/s * dr.", required=False, type=float, default=0.95)
    return (p.parse_args())

class Coordinates:
    def __init__(
        self, 
        record_index_in_file, 
        sequence_index_in_file, 
        sequence_index_in_record, 
        sequence_index_in_file_including_ambiguous, 
        sequence_index_in_record_including_ambiguous, 
        kmer
    ) -> None:
        self.record_index_in_file = record_index_in_file
        self.sequence_index_in_file = sequence_index_in_file
        self.sequence_index_in_record = sequence_index_in_record
        self.sequence_index_in_file_including_ambiguous = sequence_index_in_file_including_ambiguous
        self.sequence_index_in_record_including_ambiguous = sequence_index_in_record_including_ambiguous
        self.kmer = kmer

def read_coords(filename):
    coords = []
    with open(filename, "r") as f:
        for line in f.readlines():
            parts = line.split(",")
            coords.append(Coordinates(int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), parts[5].strip()))
    return coords

def read_macle_complexity(filename):
    complexities = {}
    first_line_read = False
    midpoint_offset = 0
    with open(filename, "r") as f:
        for line in f.readlines():
            parts = line.split()
            if not first_line_read:
                if parts[0] == "<file>":
                    continue
                midpoint_offset = (int(parts[1]))
                first_line_read = True
            new_complexity = (int(parts[1]) - midpoint_offset, float(parts[2]))
            if parts[0] in complexities:
                complexities[parts[0]].append(new_complexity)
            else:
                complexities[parts[0]] = [new_complexity]
    return (complexities, 2*midpoint_offset)

def read_gtf(filename, gene_search_terms):
    df = gtf.read_gtf(filename).to_pandas()
    relevant_cds = df[
        (df["product"].str.contains("|".join(gene_search_terms), False)) & 
        (df["feature"] == "CDS")
    ]
    relevant_genes = df[
        (df["feature"] == "gene") &
        (df["gene_id"].isin(relevant_cds.gene_id))
    ]
    return dict(
	    relevant_genes.groupby("seqname", observed=True)["start"].agg(list)
    )


def calculate_gene_count_in_windows(window_size, sequence_length, gene_positions, overlap):
    result = np.zeros(calculate_number_of_windows(window_size, sequence_length, overlap), dtype=int)
    for p in gene_positions:
        if overlap:
            index = np.arange(np.max([0, p - window_size + 1]), np.min([len(result), p+1]))
        else:
            index = p // window_size
            if index >= len(result):
                continue
        result[index] += 1
    return result

"""
    If overlap=True, then the each index in the result array is the start position
    of the window.
    If overlap=False, then the start position of the window is index*window_size
"""
def calculate_densities_in_windows(window_size, lengths, relevant_coords, overlap):
    sequence_length = lengths
    densities = np.zeros(calculate_number_of_windows(window_size, sequence_length, overlap), dtype=int)
    for c in relevant_coords:
        if overlap:
            pos = c.sequence_index_in_record_including_ambiguous
            index = np.arange(np.max([0, pos - window_size + 1]), np.min([len(densities), pos+1]))
        else:
            index = c.sequence_index_in_record_including_ambiguous // window_size
            if index >= len(densities):
                continue
        densities[index] += 1

    return densities

"""
    Returns index and the info if the position is identical to the start
    position of a window.
"""
def calculate_index(position, window_size, overlap):
    if overlap:
        return (position, True)
    return (position // window_size, position % window_size == 0)

"""
    overlap only separates between offset = 1 and offset = window_size  
"""
def calculate_number_of_windows(window_size, sequence_length, overlap):
    if overlap:
        return sequence_length - window_size + 1
    return sequence_length // window_size

""""
    Create the contingency table. Assume data is a list of (density,
    complexity). Output 2x2 table: high complexity (> 0.5) vs low complexity (<
    0.5) and usual density (w/s +/- epsilon) and unusual density (remaining
    density).
    Assume all negative complexity windows are already filtered out.
    
    w is the window size
    s is the scaling parameter
"""
def create_hash2comp_observations(data, w, s, dr):
    result = np.zeros(shape=(2, 2), dtype=int)
    
    complexity_threshold = 0.5
    e = w/s
    expected_range = ((e - e * dr), (e + e * dr))

    for d, c in data:
        if c < complexity_threshold:
            if d < expected_range[0] or d > expected_range[1]:
                #unexpected density, low complexity
                result[0][0] += 1
            else:
                result[0][1] += 1
        else:
            if d < expected_range[0] or d > expected_range[1]:
                result[1][0] += 1
            else:
                result[1][1] += 1
    
    # row 0 is low complexity, row 1 is high complexity
    # col 0 is unexpected density, col 1 is expected density
    return result

"""
    Create the contingency table. Assume data is a list of (density, effecotr
    count). Output 2x2 table: effector absent vs effector present
    and usual density (w/s +/- epsilon) and unusual density (remaining density).
    Assume all negative complexity windows are already filtered out.
    
    w is the window size 
    s is the scaling parameter
"""
def create_hash2effector_observations(data, w, s, dr):
    result = np.zeros(shape=(2, 2), dtype=int)
    
    expected_density = w/s
    expected_range = ((expected_density - expected_density * dr), (expected_density + expected_density * dr))

    for d, e in data:
        if e == 0:
            if d < expected_range[0] or d > expected_range[1]:
                #unexpected density, absent effectors
                result[0][0] += 1
            else:
                result[0][1] += 1
        else:
            if d < expected_range[0] or d > expected_range[1]:
                result[1][0] += 1
            else:
                result[1][1] += 1
    
    # row 0 is absent effectors, row 1 is present effectors
    # col 0 is unexpected density, col 1 is expected density
    return result

"""
    Splits a list of tuples of (hash_density, other_value) into a tuple of lists
    where the first list is other_value for all hash_density in the unexpected
    range and the second list is other_value for all hash_density in the
    expected range.
    Expected range is calculated by (w/s) +/- (w/s*dr)
"""
def split(data, w, s, dr):
    e = w/s
    expected_range = ((e - e * dr), (e + e * dr))
    unexpeceted_density_complexities = np.array([c for d, c in data if d < expected_range[0] or d > expected_range[1]])
    expected_density_complexities = np.array([c for d, c in data if not (d < expected_range[0] or d > expected_range[1])])
    return (unexpeceted_density_complexities, expected_density_complexities)



def main():
    args = create_parser()
    seq = mf.read(args.sequence)
    complexities, macle_window_size = read_macle_complexity(args.macle)
    if macle_window_size != args.window_size:
        print("error: window size of macle computation does not equal the input window size")
        return
    

    if (args.annotations != None):
        annotations = read_gtf(args.annotations, gene_search_terms)

    coords = [(c, read_coords(c)) for c in args.coordinates]

    record_index = 0
    current_coords_pointer = [0 for _ in coords]

    hash_density_complexity_list = []
    hash_density_effector_density_list = []
    
    for s in seq:

        macle_header = s.getHead().split()[0].strip().replace(">", "")

        if len(macle_header) > 32:
            macle_header = macle_header[:32]
        
        length = len(s.body)
        length_without_ambig = length - (s.body.count("N") + s.body.count("n"))
        
        if calculate_number_of_windows(args.window_size, length_without_ambig, args.overlap) <= 0:
            print(f"warning: {macle_header} is too short for density calculations with window size {args.window_size}", file=sys.stderr)
            record_index += 1
            continue      
        
        current_densities = []
        for i in range(len(coords)):
            path, c = coords[i]
            coords_start = 0
            coords_end = 0
            while current_coords_pointer[i] < len(c) and c[current_coords_pointer[i]].record_index_in_file < record_index:
                current_coords_pointer[i] += 1

            if current_coords_pointer[i] >= len(c) or c[current_coords_pointer[i]].record_index_in_file > record_index:
                print(f"warning: no coordinates for {macle_header} found in {path}", file=sys.stderr)
            else:
                coords_start = current_coords_pointer[i]
                while current_coords_pointer[i] < len(c) and c[current_coords_pointer[i]].record_index_in_file == record_index:
                    current_coords_pointer[i] += 1
                coords_end = current_coords_pointer[i] #end is always exclusive, so this works
            
            window_density = calculate_densities_in_windows(args.window_size, length, c[coords_start: coords_end], args.overlap)  
            current_densities.append(window_density)

        record_index += 1
        median_density = np.median(current_densities, axis=0)

        # likely that a sequence does not contain an annotation, so no warning
        # here
        effector_densities = []
        if args.annotations != None and macle_header in annotations:
            effector_densities = calculate_gene_count_in_windows(args.window_size, length, annotations[macle_header], args.overlap)            
        
        if (macle_header not in complexities):
                print(f"warning: no complexities for {macle_header} found", file=sys.stderr)
                continue
        else:
            for pos, c_m in complexities[macle_header]:
                # only include windows with non-negative macle complexity and
                # only those macle complexities that can be mapped to a density
                # window (with smaller macle window sizes, some windows at the
                # end cannot be mapped)
                index, is_at_window_start = calculate_index(pos, args.window_size, args.overlap)
                if c_m >= 0 and len(median_density) > index and is_at_window_start:    
                    hash_density_complexity_list.append((median_density[index], c_m))
                    if len(effector_densities) > index:
                        hash_density_effector_density_list.append((median_density[index], effector_densities[index]))

    # Hash Density to Complexity analysis
    x = [x_i for x_i, _ in hash_density_complexity_list]
    y = [y_i for _, y_i in hash_density_complexity_list]
    r=np.corrcoef(x, y, rowvar=True)

    observations = create_hash2comp_observations(hash_density_complexity_list, args.window_size, args.scaling, args.density_range)
    if(len(np.flatnonzero(observations)) == len(observations.flatten())):
        chi2_result = sc.stats.chi2_contingency(observations)
    else:
        chi2_result = None
    unexpected_density_complexities, expected_density_complexities = split(hash_density_complexity_list, args.window_size, args.scaling, args.density_range)
    if len(unexpected_density_complexities) > 0 and len(expected_density_complexities) > 0:
        mannwhitneyu_results = sc.stats.mannwhitneyu(unexpected_density_complexities, expected_density_complexities)
    else:
        mannwhitneyu_results = None
    
    print()
    print("Statistics for hash count vs complexity of window")
    print(f"Number of windows analyzed: {len(hash_density_complexity_list)}")
    def overlapping_string():
        if not args.overlap:
            return "not "
        return ""
    
    print(f"windows are {overlapping_string()}overlapping")
    print(f"correlation of window density vs complexity r={r[0][1]}")
    if chi2_result == None:
        print("at least one observation was empty, no χ2 performed")
    else:
        print(f"χ2({chi2_result.dof}, N={np.sum(observations)})={chi2_result.statistic}, p={chi2_result.pvalue}")

    if mannwhitneyu_results == None:
        print("at least one observation was empty, no Mann-Whitney-Test performed")
    else:
        print(f"Mann-Whitney-Test: z={mannwhitneyu_results.statistic}, p={mannwhitneyu_results.pvalue}")

    plt.scatter(x, y, s=0.1)
    plt.xlabel(f"hashes in window with size $w={args.window_size}$")    
    plt.ylabel(f"$C_m$ in window with size $w={macle_window_size}$")          
    plt.show()

    # Hash Density to Effector density analysis
    if args.annotations:
        x = [x_i for x_i, _ in hash_density_effector_density_list]
        y = [y_i for _, y_i in hash_density_effector_density_list]
        r=np.corrcoef(x, y, rowvar=True)

        observations = create_hash2effector_observations(hash_density_effector_density_list, args.window_size, args.scaling, args.density_range)
        if(len(np.flatnonzero(observations)) == len(observations.flatten())):
            chi2_result = sc.stats.chi2_contingency(observations)
        else:
            chi2_result = None
        
        chi2_result = sc.stats.chi2_contingency(observations)

        unexpected_density_effectors, expected_density_effectors = split(hash_density_effector_density_list, args.window_size, args.scaling, args.density_range)
        if len(unexpected_density_effectors) > 0 and len(expected_density_effectors) > 0: 
            mannwhitneyu_results = sc.stats.mannwhitneyu(unexpected_density_effectors, expected_density_effectors)
        else:
            mannwhitneyu_results = None
        
        print()
        print("Statistics for hash count vs effector gene count")
        print(f"Number of windows analyzed: {len(hash_density_complexity_list)}")
        def overlapping_string():
            if not args.overlap:
                return "not "
            return ""
        
        print(f"windows are {overlapping_string()}overlapping")
        print(f"correlation of window density vs effector count r={r[0][1]}")
        if chi2_result == None:
            print("at least one observation was empty, no χ2 performed")
        else:
            print(f"χ2({chi2_result.dof}, N={np.sum(observations)})={chi2_result.statistic}, p={chi2_result.pvalue}")
        if mannwhitneyu_results == None:
            print("at least one observation was empty, no Mann-Whitney-Test performed")
        else:
            print(f"Mann-Whitney-Test: z={mannwhitneyu_results.statistic}, p={mannwhitneyu_results.pvalue}")

        plt.scatter(x, y, s=0.1)
        plt.xlabel(f"hashes in window with size $w={args.window_size}$")    
        plt.ylabel(f"effector genes in window with size $w={macle_window_size}$")          
        plt.show()
main()
