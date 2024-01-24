import numpy as np
import matplotlib.pyplot as plt

def findGapsLargerThan(sortedList, threshold):
    result = []
    for i in range(1, len(sortedList)):
        delta = sortedList[i] - sortedList[i-1]
        if delta > threshold:
            result.append((delta, i, sortedList[i]))
    return result

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

pos = []
pos_ambig = []
with open("fmhdist/PhyInf1.fna.gz.sketch.coordinates", "r") as f:
    for line in f.readlines():
        parts = line.split(":")
        current = []
        for part in parts:
            current.append(part.split(",")[1].replace(")", "").strip())
        pos.append(current[0])
        pos_ambig.append(current[1])

window_size = 10000
density_threshold = 20
gap_threshold = 20000
pos = np.array(pos, dtype=np.uint)
pos_ambig = np.array(pos_ambig, dtype=np.uint)

windows = collectWindows(pos_ambig, window_size)


largerThanT = findGapsLargerThan(pos_ambig, gap_threshold)
max_gap_size = np.max([size for size, _, _ in largerThanT])
min_gap_size = np.min([size for size, _, _ in largerThanT])

densestWindows = findDensestWindows(windows, density_threshold)
densestWindowsMaxDensity = np.max([size for _, size in densestWindows])
densestWindowsMinDensity = np.min([size for _, size in densestWindows])

plt.scatter(pos, len(pos) * [0], s=1, c="blue", label="hash origin position excluding ambiguous $k$-mers")
plt.scatter(pos_ambig, len(pos_ambig) * [0.01], s=1, c="orange", label="hash origin position including ambigouos $k$-mers")

plt.scatter([pos for _, _, pos in largerThanT], len(largerThanT)* [0.11], s=[(size-min_gap_size)/(max_gap_size-min_gap_size) * 100 for size, _, _ in largerThanT], c="green", label=f"gaps larger than $t={gap_threshold}$")
plt.scatter([pos for pos, _ in densestWindows], len(densestWindows)*[-0.11], s=[(density-densestWindowsMinDensity)/(densestWindowsMaxDensity-densestWindowsMinDensity) * 100 for _, density in densestWindows], c="black", label=f"windows with size $s={window_size}$ that are denser than $t={density_threshold}$")

plt.xlim(0, np.max(pos_ambig))
plt.ylim(-1, 1)
plt.gca().get_yaxis().set_visible(False)
plt.legend()
plt.xlabel("genome position of $k$-mers in $S_{FRAC}$")
plt.show()