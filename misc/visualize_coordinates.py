import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def findGapsLargerThan(sortedList, threshold):
    result = []
    for i in range(1, len(sortedList)):
        delta = sortedList[i] - sortedList[i-1]
        if delta > threshold:
            result.append((delta, i, sortedList[i]))
    return result

def collectWindowForRow(row, df, window_size):
    start_pos = row["sequenceIndexInFile"]
    window = df[(df["sequenceIndexInFile"] < start_pos + window_size) & (df["sequenceIndexInFile"] >= start_pos)]
    size = len(window)
    counts = window["kmer"].value_counts()
    unique_ratio = len(counts[counts == 1]) / size
    return (start_pos, size, unique_ratio)

def collectWindowsForDf(sorted_df, window_size):
    rows = []
    for i in range(len(sorted_df)):
        current_row = sorted_df.iloc[i, :]
        window_start = current_row["sequenceIndexInFile"]
        kmers = [current_row["kmer"]]
        j = i + 1
        while j < len(sorted_df) and sorted_df.iloc[j, :]["sequenceIndexInFile"] - window_start < window_size:
            kmers.append(sorted_df.iloc[j, :]["kmer"])
            j += 1
        unique_kmers = [kmer for kmer in kmers if kmers.count(kmer) == 1]
        current_window = [window_start, len(kmers), len(unique_kmers)/len(kmers)]
        rows.append(current_window)
    return pd.DataFrame(rows, columns=["window_start", "kmer_count", "unique_kmer_ratio"])

def collectWindows(sorted_list, window_size):
    result = []
    for i in range(len(sorted_list)):
        current_window = [sorted_list[i]]
        j = i + 1
        while j < len(sorted_list) and sorted_list[j] - sorted_list[i] < window_size:
            current_window.append(sorted_list[j])
            j += 1
        result.append(current_window)
    return result

def findDensestWindows(window_list, threshold):
    result = []
    for window in window_list:
        if len(window) > threshold:
            result.append((window[0], len(window)))
    return result

df = pd.read_csv("fmhdist/PhyInf1.fna.gz.sketch.coordinates", names=[
    "recordIndex",
    "sequenceIndexInFile",
    "sequenceIndexInRecord",
    "sequenceIndexInFileIncludingAmbiguous",
    "sequenceIndexInRecordIncludingAmbiguous",
    "kmer"
    ])

df["kmer"] = df["kmer"].astype("string")
counts = df["kmer"].value_counts()
df["counts"] = df.apply(lambda row: counts[row["kmer"]],axis=1)

pos = df["sequenceIndexInFile"].to_numpy()
pos_ambig = df["sequenceIndexInFileIncludingAmbiguous"].to_numpy()
counts = df["counts"].to_numpy()

window_size = 10000
density_threshold = 20
gap_threshold = 20000
pos = np.array(pos, dtype=np.uint)
pos_ambig = np.array(pos_ambig, dtype=np.uint)

windows = collectWindowsForDf(df, window_size)

largerThanT = findGapsLargerThan(pos_ambig, gap_threshold)
max_gap_size = np.max([size for size, _, _ in largerThanT])
min_gap_size = np.min([size for size, _, _ in largerThanT])

densestWindows = windows[windows["kmer_count"] > density_threshold]
densestWindowsMaxDensity = np.max(densestWindows["kmer_count"].to_numpy())
densestWindowsMinDensity = np.min(densestWindows["kmer_count"].to_numpy())

plt.scatter(pos, len(pos) * [0], s=1, c="blue", label="hash origin position excluding ambiguous $k$-mers")
plt.scatter(pos_ambig, len(pos_ambig) * [0.01], s=1, c="orange", label="hash origin position including ambigouos $k$-mers")

plt.scatter([pos for _, _, pos in largerThanT], len(largerThanT)* [0.11], s=[(size-min_gap_size)/(max_gap_size-min_gap_size) * 100 for size, _, _ in largerThanT], c="green", label=f"gaps larger than $t={gap_threshold}$")
plt.scatter(densestWindows["window_start"], len(densestWindows)*[-0.11], s=(densestWindows["kmer_count"].to_numpy()-densestWindowsMinDensity)/(densestWindowsMaxDensity-densestWindowsMinDensity) * 100, c=densestWindows["unique_kmer_ratio"], label=f"windows with size $s={window_size}$ that are denser than $t={density_threshold}$")

plt.xlim(0, np.max(pos_ambig))
plt.ylim(-1, 1)
plt.gca().get_yaxis().set_visible(False)
plt.legend()
plt.xlabel("genome position of $k$-mers in $S_{FRAC}$")
plt.show()