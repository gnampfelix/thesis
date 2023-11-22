import matplotlib.pyplot as plt
import numpy as np

arr = np.loadtxt("misc/phytophthera_genome_sizes");
plt.scatter(arr, np.random.rand(arr.size), s=2)
plt.vlines([arr.mean(), np.median(arr)], 0, 1, colors=["red", "green"])
plt.show();

print(arr.mean(), np.median(arr), arr.min(), arr.max())

# filtered_arr = arr[arr < 500000];
# plt.scatter(filtered_arr, np.random.rand(filtered_arr.size), s=0.2)
# plt.show();

# filtered_arr = arr[arr < 100000];
# plt.scatter(filtered_arr, np.random.rand(filtered_arr.size), s=0.2)
# plt.vlines([arr.mean(), np.median(arr)], 0, 1, colors=["red", "green"])
# plt.show();