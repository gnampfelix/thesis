import argparse
import sys
import miniFasta as mf
import matplotlib.pyplot as plt
import numpy as np

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-w", "--window-size", help="The window size for density analysis", type=int, default=1000, required=False)
    p.add_argument("-c", "--coordinates", help="The coordinates file to analyze", required=True, nargs="+")
    p.add_argument("-s", "--sequence", help="The underlying sequence file in FASTA format", required=True)
    p.add_argument("-m", "--macle", help="The path to the macle file", required=True)
    p.add_argument("-o", "--output", help="The path where the plots should be saved", required=False, default=".")
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


def calculate_densities_in_windows(window_size, lengths, relevant_coords):
    sequence_length = lengths
    densities = np.zeros(calculate_number_of_windows(window_size, sequence_length), dtype=int)
    for c in relevant_coords:
        pos = c.sequence_index_in_record_including_ambiguous
        densities[np.arange(np.max([0, pos - window_size + 1]), np.min([len(densities), pos+1]))] += 1

    return densities

def calculate_number_of_windows(window_size, sequence_length):
    return sequence_length - window_size + 1

def main():
    args = create_parser()
    seq = mf.read(args.sequence)
    complexities, macle_window_size = read_macle_complexity(args.macle)
    if macle_window_size != args.window_size:
        print("error: window size of macle computation does not equal the input window size")
        return
    
    coords = [(c, read_coords(c)) for c in args.coordinates]

    record_index = 0
    current_coords_pointer = [0 for _ in coords]

    densities = []
    
    for s in seq:

        macle_header = s.getHead().split()[0].strip().replace(">", "")

        if len(macle_header) > 32:
            macle_header = macle_header[:32]
        
        length = len(s.body)
        length_without_ambig = length - (s.body.count("N") + s.body.count("n"))
        
        if calculate_number_of_windows(args.window_size, length_without_ambig) <= 0:
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
            
            window_density = calculate_densities_in_windows(args.window_size, length, c[coords_start: coords_end])  
            current_densities.append(window_density)

        record_index += 1
        median_density = np.median(current_densities, axis=0)

        if (macle_header not in complexities):
                print(f"warning: no complexities for {macle_header} found", file=sys.stderr)
                continue
        else:
            for pos, c_m in complexities[macle_header]:
                # only include windows with non-negative macle complexity and
                # only those macle complexities that can be mapped to a density
                # window (with smaller macle window sizes, some windows at the
                # end cannot be mapped)
                if c_m >= 0 and len(median_density) > pos:    
                    densities.append((median_density[pos], c_m))

    x = [x_i for x_i, _ in densities]
    y = [y_i for _, y_i in densities]
    r=np.corrcoef(x, y, rowvar=True)
    print(f"correlation of window density vs complexity r={r[0][1]}")

    plt.scatter(x, y, s=0.1)
    plt.xlabel(f"hashes in window with size $w={args.window_size}$")    
    plt.ylabel(f"$C_m$ in window with size $w={macle_window_size}$")          
    plt.show()
main()
