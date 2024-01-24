import numpy as np
import matplotlib.pyplot as plt


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


pos = np.array(pos, dtype=np.uint)
pos_ambig = np.array(pos_ambig, dtype=np.uint)

plt.scatter(pos, len(pos) * [0], s=1, c="blue", label="excluding ambiguous $k$-mers")
plt.scatter(pos_ambig, len(pos_ambig) * [0.01], s=1, c="orange", label="including ambigouos $k$-mers")
plt.xlim(0, np.max(pos_ambig))
plt.ylim(-1, 1)
plt.gca().get_yaxis().set_visible(False)
plt.legend()
plt.xlabel("genome position of $k$-mers in $S_{FRAC}$")
plt.show()