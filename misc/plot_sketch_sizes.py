import matplotlib.pyplot as plt
import numpy as np

arr = np.loadtxt("misc/phytophthera_sketch_sizes_k21_s2000");
plt.hist(arr, bins=40, cumulative=True);
plt.show()



plt.scatter(arr, np.random.rand(arr.size), s=2)
plt.vlines([arr.mean(), np.median(arr)], 0, 1, colors=["red", "green"])
plt.show();

print(arr.mean(), np.median(arr), arr.min(), arr.max())
