import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

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


coords = []
with open("fmhdist/PhyInf1.fna.gz.sketch.coordinates", "r") as f:
    for line in f.readlines():
        parts = line.split(",")
        coords.append(Coordinates(int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), parts[5].strip()))
    
def findGapsLargerThan(c, threshold):
    result = []
    for i in range(1, len(c)):
        delta = c[i].sequence_index_in_file_including_ambiguous - c[i-1].sequence_index_in_file_including_ambiguous
        if delta > threshold:
            result.append((delta, i-1, c[i-1].sequence_index_in_file_including_ambiguous))
    return result

def collect_windows(c, window_size):
    windows = []
    for i in range(len(c)):
        current_kmers = [c[i].kmer]
        j = i + 1
        while (
            j < len(c) and 
            c[j].sequence_index_in_file_including_ambiguous - c[i].sequence_index_in_file_including_ambiguous < window_size
        ):
            current_kmers.append(c[j].kmer)
            j += 1
        unique_kmers = list(Counter(current_kmers).keys())
        windows.append((i, c[i].sequence_index_in_file_including_ambiguous, current_kmers, unique_kmers))
    return windows


window_size = 10000
density_threshold = 15
gap_threshold = 40000

pos = [c.sequence_index_in_file_including_ambiguous for c in coords]
windows = collect_windows(coords, window_size)
densest_windows_unique = [(i, p, k, u) for i, p, k, u in windows if len(u) > density_threshold]
gaps = findGapsLargerThan(coords, gap_threshold)
min_gap_size = np.min([d for d, _, _ in gaps])
max_gap_size = np.max([d for d, _, _ in gaps])

# print("densest windows including repeats")
# for w in densest_windows:
#     print(f"window starting at {w[1]}")
#     for kmer in w[2]:
#         print(kmer)
#     print()
#plt.scatter([p for _, p, _, _ in densest_windows], len(densest_windows) * [0.1], s=[len(k) for _, _, k, _ in densest_windows], c="black")

print("densest windows excluding repeats")
for w in densest_windows_unique:
    print(f"window starting at {w[1]}")
    for kmer in w[3]:
        print(kmer)
    print()

print()
print("results")
print(f"there are {len(gaps)} gaps that are longer than {gap_threshold}:")
print([size for size, _, _ in gaps])
print()
print(f"there are {len(densest_windows_unique)} windows of size {window_size} with more than {density_threshold} unique k-mers in it")
print("-see window content above-")
print()


plt.scatter(
    [p for _, p, _, _ in densest_windows_unique], 
    len(densest_windows_unique) * [0.1], 
    s=[len(u) for _, _, _, u in densest_windows_unique], 
    c="black", 
    label=f"position of windows of size {window_size} with more than {density_threshold} unique $k$-mers in it"
)

plt.scatter(pos, len(pos) * [0.00], s=1, c="orange", label="hash origin position including ambigouos $k$-mers")

plt.scatter(
    [pos for _, _, pos in gaps], 
    len(gaps)* [-0.1], 
    s=[(size-min_gap_size)/(max_gap_size-min_gap_size) * 100 for size, _, _ in gaps], 
    c="green", 
    label=f"start position of gaps larger than $t={gap_threshold}$"
)

plt.xlim(0, np.max(pos))
plt.ylim(-1, 1)
plt.gca().get_yaxis().set_visible(False)
plt.legend()
plt.xlabel("genome position of $k$-mers in $S_{FRAC}$")
plt.show()