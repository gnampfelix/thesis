import argparse
import sys
import miniFasta as mf
import matplotlib.pyplot as plt
import numpy as np

"""
Script to plot hash counts of a given window and the sequence complexity of the
same window. Plots are generated per sequence in the FASTA file, but only if
there are enough windows in that sequence (i.e. at least 1 window / at least -wt
windows). For each window, the lowest hash count, the highest hash count, the
median hash count and the average hash count are depicted (assumption: the user
passes different coordinate files, e.g. generated by using different hash
seeds). A complexity below 0 indicates that the window contains ambiguous bases.
However, those windows are still plotted to paint a complete picture. 
"""


def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-w", "--window-size", help="The window size for density analysis", type=int, default=1000, required=False)
    p.add_argument("-sc", "--scaling", help="The scaling parameter s that was used to generate the sketch on which the coordinates are based", type=int, default=2000)
    p.add_argument("-c", "--coordinates", help="The coordinates file to analyze", required=True, nargs="+")
    p.add_argument("-s", "--sequence", help="The underlying sequence file in FASTA format", required=True)
    p.add_argument("-m", "--macle", help="The path to the macle file", required=True)
    p.add_argument("-o", "--overlap", help="Indicates if the window offset should be 1, otherwise <window-size>", action="store_true")
    p.add_argument("-dr", "--density-range", help="The range in which the density in a window is considered 'expected' is given by w/s +/- w/s * dr.", required=False, type=float, default=0.95)
    p.add_argument("-i", "--interactive", help="Show plots", action="store_true")
    p.add_argument("-wt", "--window-threshold", help="minimum number of windows to be analyzed for a sequence befor plot is generated", type=int, default=1)
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
def is_window_start(position, window_size, overlap):
    return overlap or position % window_size == 0

"""
    overlap only separates between offset = 1 and offset = window_size  
"""
def calculate_number_of_windows(window_size, sequence_length, overlap):
    if overlap:
        return sequence_length - window_size + 1
    return sequence_length // window_size

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
    new_rc_params = {'text.usetex': False,
        "svg.fonttype": 'none',
        "font.size": 20
    }
    plt.rcParams.update(new_rc_params)
    args = create_parser()
    seq = mf.read(args.sequence)
    complexities, macle_window_size = read_macle_complexity(args.macle)
    if macle_window_size != args.window_size:
        print("error: window size of macle computation does not equal the input window size")
        return
    

    coords = [(c, read_coords(c)) for c in args.coordinates]

    record_index = 0
    current_coords_pointer = [0 for _ in coords]

    
    for s in seq:
        hash_density_complexity_list = []
        hash_density_effector_density_list = []

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
        median_densities = np.median(current_densities, axis=0)
        mean_densities = np.mean(current_densities, axis=0)
        max_densities = np.max(current_densities, axis=0)
        min_densities = np.min(current_densities, axis=0)

        window_index = [i for i in range(len(mean_densities))]
        if (macle_header not in complexities):
            print(f"warning: no complexities for {macle_header} found")
        elif (len(median_densities) < args.window_threshold):
            print(f"warning: only {len(median_densities)} window analyzed for {macle_header}, no plot shown")
        else:
            print(f"showing plot for {macle_header}")
            fig, ax1 = plt.subplots()
            #ax1.set_title(f"hash counts and sequence complexity for sequence")

            p1 = ax1.plot(window_index, mean_densities, label="mean hash count", color="darkgray")
            p2 = ax1.plot(window_index, median_densities, label="median hash count", color="black")
            p3 = ax1.fill_between(window_index, min_densities, max_densities, facecolor="lightgray", label="hash count range")

            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            relevant_complexities = [comp for index, comp in complexities[macle_header] if is_window_start(index, args.window_size, args.overlap)]
            p4 = ax2.plot(window_index, relevant_complexities, color="forestgreen", label="\\$C_m\\$ complexity")

            ax1.set_ylabel("hashes in window")
            ax1.set_xlabel(f"window index with \\$w={args.window_size}\\$")
            
            ax1.tick_params(axis="y")
            ax1.set_ylim([0, 20])
            ax2.tick_params(axis="y")

            ax2.set_ylim([0, 1.1])
            ax2.set_ylabel("\\$C_m\\$") 

            all_plots = p1 + p2 # Those two are lists, concatenate
            all_plots.append(p3) # this is not a list, append
            if p4:
                all_plots += p4
            labels = [l.get_label() for l in all_plots]            

            ax2.legend(handles=all_plots, labels=labels, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol= 2)
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            #fig.savefig(f"{args.output}/{macle_header}.png")
            if args.interactive:
                plt.show()
            plt.close(fig)
               
main()
