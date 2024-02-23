import argparse
import miniFasta as mf
import matplotlib.pyplot as plt
import numpy as np

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-w", "--window-size", help="The window size for density analysis", type=int, default=1000, required=False)
    p.add_argument("-k", "--window-interval", help="The window interval for density analysis", type=int, default=100, required=False)
    p.add_argument("-c", "--coordinates", help="The coordinates file to analyze", required=True, nargs="+")
    p.add_argument("-s", "--sequence", help="The underlying sequence file in FASTA format", required=True)
    p.add_argument("-m", "--macle", help="The path to the macle file", required=True)
    p.add_argument("-rm", "--repeat-masker", help="The path to the RepeatMasker .out file", required=True)
    p.add_argument("-o", "--output", help="The path where the plots should be saved", required=False, default=".")
    p.add_argument("-mw", "--meta-window-size", help="Number of consecutive windows with high or 0 density to trigger the image save for a sequence", type=int, required=False, default=200)
    p.add_argument("-i", "--ignore", help="List of header names to ignore in the analysis (e.g. because those sequences contain amb. nucleotides)", required=False, default="")
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

def read_macle_complexity(filename, midpoint_offset):
    complexities = {}
    with open(filename, "r") as f:
        for line in f.readlines():
            parts = line.split()
            new_complexity = (int(parts[1]) - midpoint_offset, float(parts[2]))
            if parts[0] in complexities:
                complexities[parts[0]].append(new_complexity)
            else:
                complexities[parts[0]] = [new_complexity]
    return complexities

def read_repeat_masker_output(filename):
    intervals = {}
    with open(filename, "r") as f:
        # skip the first three lines as headers
        for line in f.readlines()[3:]:
            parts = line.split()
            new_interval = (int(parts[5]), int(parts[6]))
            if parts[4] in intervals:
                intervals[parts[4]].append(new_interval)
            else:
                intervals[parts[4]] = [new_interval]
    return intervals

"""
Calculate the number of masked positions in a given window
"""
def calculate_masked_nucleotide_ratio(window_size, sequence_length, relevant_rm_data):
    windows = np.zeros(calculate_number_of_windows(window_size, sequence_length), dtype=int)
    for start, end in relevant_rm_data:
        # We _could_ calculate things here using closed form. However, I am not
        # able to enumerate all possible cases here, so I will just iterate the
        # nucleotides in each window.
        masked_nucleotides_in_window = np.zeros(len(windows), dtype=int)
        for i in range(start, end + 1):
            masked_nucleotides_in_window[np.arange(np.max([0, i-window_size]), np.min([i+1, len(masked_nucleotides_in_window)]))] += 1
        windows += masked_nucleotides_in_window
    return windows/window_size
"""
meta_window_size: The _number_ of consecutive windows that must be "interesting".
assumption: window start positions are offset by 1
"""
def has_interesting_density(densities, meta_window_size):
    np_densities = densities
    zero_windows = np.where(np_densities == 0)[0]
    for i in range(len(zero_windows)-meta_window_size):
        if zero_windows[i] + meta_window_size == zero_windows[i+meta_window_size]:
            return True 
    return False

def calculate_densities_in_windows(window_size, sequence_length, relevant_coords, use_ambig):
    densities = np.zeros(calculate_number_of_windows(window_size, sequence_length), dtype=int)
    for c in relevant_coords:
        if use_ambig:
            pos = c.sequence_index_in_record_including_ambiguous
        else:
            pos = c.sequence_index_in_record
        current_offset = np.zeros(len(densities), dtype=int)
        # All window_size windows starting before pos include the k-mer, and the window starting with pos itself!
        current_offset[np.arange(np.max([0, pos - window_size + 1]), np.min([len(densities), pos+1]))] = 1
        densities += current_offset
    return densities

def calculate_number_of_windows(window_size, sequence_length):
    return sequence_length - window_size + 1

def read_ignore_list(filename):
    if filename == "":
        return []
    print(f"using ignore list {filename}")
    with open(filename, "r") as f:
        return f.readlines()

def main():
    args = create_parser()
    midpoint_offset = args.window_size // 2
    seq = mf.read(args.sequence)
    complexities = read_macle_complexity(args.macle, midpoint_offset)
    repeat_masker_intervals = read_repeat_masker_output(args.repeat_masker)
    coords = [(c, read_coords(c)) for c in args.coordinates]

    ignore_list = [i.strip() for i in read_ignore_list(args.ignore)]

    record_index = 0
    current_coords_pointer = [0 for _ in coords]
    for s in seq:

        macle_header = s.getHead().split()[0].strip().replace(">", "")

        if len(macle_header) > 32:
            macle_header = macle_header[:32]

        if macle_header in ignore_list:
            print(f"ignoring {macle_header}...")
            continue
        
        length = len(s.body)
        length_without_ambig = length - (s.body.count("N") + s.body.count("n"))
        
        if calculate_number_of_windows(args.window_size, length_without_ambig) <= 0:
            print(f"warning: {macle_header} is too short for density calculations with window size {args.window_size}")
            record_index += 1
            continue      
        
        print(f"analyzing {macle_header}...")
        densities = []
        plot_densities = []
        for i in range(len(coords)):
            path, c = coords[i]
            coords_start = 0
            coords_end = 0
            while current_coords_pointer[i] < len(c) and c[current_coords_pointer[i]].record_index_in_file < record_index:
                current_coords_pointer[i] += 1

            if current_coords_pointer[i] >= len(c) or c[current_coords_pointer[i]].record_index_in_file > record_index:
                print(f"warning: no coordinates for {macle_header} found in {path}")
            else:
                coords_start = current_coords_pointer[i]
                while current_coords_pointer[i] < len(c) and c[current_coords_pointer[i]].record_index_in_file == record_index:
                    current_coords_pointer[i] += 1
                coords_end = current_coords_pointer[i] #end is always exclusive, so this works
            
            # Ignore all amb. k-mers to determine "what is interesting", but
            # include them for plotting
            densities.append(calculate_densities_in_windows(args.window_size, length_without_ambig, c[coords_start: coords_end], False))
            plot_densities.append(calculate_densities_in_windows(args.window_size, length, c[coords_start: coords_end], True))
     
        record_index += 1
        is_interesting = False
        for d in densities:
            if has_interesting_density(d, args.meta_window_size):
                is_interesting = True
                break
      
        if is_interesting:
            max_densities = np.max(plot_densities, axis=0)
            min_densities = np.min(plot_densities, axis=0)
            mean_densities = np.mean(plot_densities, axis=0)
            median_densities = np.median(plot_densities, axis=0)

            fig, ax1 = plt.subplots()
            ax1.set_title(macle_header)

            p1 = ax1.plot(mean_densities, label="mean #hashes", color="darkgray")
            p2 = ax1.plot(median_densities, label="median #hashes", color="black")
            p3 = ax1.fill_between(np.arange(0, len(max_densities)), min_densities, max_densities, facecolor="lightgray", label="#hash range")

            
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            if (macle_header not in complexities):
                print(f"warning: no complexities for {macle_header} found")
            else:
                p4 = ax2.plot([pos for pos, _ in complexities[macle_header]], [c for _, c in complexities[macle_header]], color="forestgreen", label="$C_m$ complexity")

            if macle_header in repeat_masker_intervals:
                rm_windows = calculate_masked_nucleotide_ratio(args.window_size, len(s), repeat_masker_intervals[macle_header])
                p5 = ax2.plot(np.arange(0, len(rm_windows)), rm_windows, color="salmon", label="masked nucleotide ratio")
            

            ax1.set_ylabel("#hashes in window")
            ax1.set_xlabel(f"window start position in genome with $w={args.window_size}$")
            
            ax1.tick_params(axis="y")
            ax1.set_ylim([0, 20])
            ax2.tick_params(axis="y")

            ax2.set_ylim([0, 1])
            ax2.set_ylabel("$C_m$/Masked nucleotides") 

            all_plots = p1 + p2 # Those two are lists, concatenate
            all_plots.append(p3) # this is not a list, append
            if p4:
                all_plots += p4
            if p5:
                all_plots += p5
            labels = [l.get_label() for l in all_plots]            

            ax2.legend(handles=all_plots, labels=labels, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol= 2)
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            fig.savefig(f"{args.output}/{macle_header}.png")
            plt.close(fig)
main()
