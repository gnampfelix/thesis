import argparse
import miniFasta as mf
import matplotlib.pyplot as plt

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("-w", "--window-size", help="The window size for density analysis", type=int, default=1000, required=False)
    p.add_argument("-k", "--window-interval", help="The window interval for density analysis", type=int, default=100, required=False)
    p.add_argument("-c", "--coordinates", help="The coordinates file to analyze", required=True)
    p.add_argument("-s", "--sequence", help="The underlying sequence file in FASTA format", required=True)
    p.add_argument("-m", "--macle", help="The path to the macle file", required=True)
    p.add_argument("-o", "--output", help="The path where the plots should be saved", required=False, default=".")
    p.add_argument("-mw", "--meta-window-size", help="Number of consecutive windows with high or 0 density to trigger the image save for a sequence", type=int, required=False, default=200)
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
        self.sequence_index_in_refcord_including_ambiguous = sequence_index_in_record_including_ambiguous
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

"""
meta_window_size: The _number_ of consecutive windows that must be "interesting".
"""
def has_interesting_density(densities, meta_window_size):
    for i in range(len(densities)):
        if densities[i][1] == 0:
            aborted = False
            for j in range(meta_window_size):
                if i+j+1 >= len(densities) or densities[j+i+1][1] != 0:
                    aborted = True
                    break
            if not aborted:
                return True
    return False



def calculate_densities_in_windows(window_size, window_interval, sequence_length, relevant_coords):
    densities = []
    smallest_coords_index = 0
    for window_start in range(0, sequence_length-window_size):
        window_density = 0
        for current_coords_index in range(smallest_coords_index, len(relevant_coords)):
            if relevant_coords[current_coords_index].sequence_index_in_record < window_start:
                smallest_coords_index += 1
                continue
            if relevant_coords[current_coords_index].sequence_index_in_record > window_start + window_size:
                break
            window_density += 1
        densities.append((window_start, window_density))
    return densities

def main():
    args = create_parser()
    midpoint_offset = args.window_size // 2
    seq = mf.read(args.sequence)
    complexities = read_macle_complexity(args.macle, midpoint_offset)
    coords = read_coords(args.coordinates)

    record_index = 0
    current_coords_pointer = 0
    for s in seq:
        no_coords = False
        coords_start = 0
        coords_end = 0
        while current_coords_pointer < len(coords) and coords[current_coords_pointer].record_index_in_file < record_index:
            current_coords_pointer += 1
        
        macle_header = s.getHead().split()[0].strip().replace(">", "")
        if len(macle_header) > 32:
            macle_header = macle_header[:32]

        if coords[current_coords_pointer].record_index_in_file > record_index:
            print(f"warning: no coordinates for {macle_header} found")
            no_coords = True
        else:
            coords_start = current_coords_pointer
            while current_coords_pointer < len(coords) and coords[current_coords_pointer].record_index_in_file == record_index:
                current_coords_pointer += 1
            coords_end = current_coords_pointer #end is always exclusive, so this works
     
        record_index += 1
      
        densities = calculate_densities_in_windows(args.window_size, args.window_interval, len(s), coords[coords_start: coords_end])
        if has_interesting_density(densities, args.meta_window_size):
            fig, ax1 = plt.subplots()
            ax1.set_title(macle_header)

            ax1.set_xlabel(f"window start in genome with $w={args.window_size}$")
            ax1.set_ylabel("number of hashes in the sketch in the window", color="forestgreen")
            ax1.plot([pos for pos, _ in densities], [dens for _, dens in densities], color="forestgreen")
            ax1.tick_params(axis="y", color="forestgreen")
            ax1.set_ylim([0, 20])
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            ax2.set_ylabel("$C_m$ value of window", color="salmon")  # we already handled the x-label with ax1
            if (macle_header not in complexities):
                print(f"warning: no complexities for {macle_header} found")
            else:
                ax2.plot([pos for pos, _ in complexities[macle_header]], [c for _, c in complexities[macle_header]], color="salmon")
            
            ax2.tick_params(axis="y", color="salmon")
            ax2.set_ylim([0, 1])

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            fig.savefig(f"{args.output}/{macle_header}.png")
            plt.close(fig)
main()
