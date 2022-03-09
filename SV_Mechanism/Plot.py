import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import BuildSVAndTestData
from matplotlib.patches import Rectangle


"""

RAG

"""
ax1 = plt.subplot(221)

sv_scores_single = []
random_scores_single = []

sv_rag_data = BuildSVAndTestData.read_FIMO_output("..\\Resources\\SV.breakpoints.padding.25.FIMO.tsv")
for sv_id in sv_rag_data:
    for rag_score in sv_rag_data[sv_id]:
        sv_scores_single.append(rag_score)

random_rag_data = BuildSVAndTestData.read_FIMO_output("..\\Resources\\Random.breakpoints.padding.25.FIMO.tsv")
for sv_id in random_rag_data:
    for rag_score in random_rag_data[sv_id]:
        random_scores_single.append(rag_score)

sns.boxplot(data=[sv_scores_single, random_scores_single], ax=ax1)
sns.stripplot(data=[sv_scores_single, random_scores_single], ax=ax1, s=2.5, color=".25")

ax1.plot([0, 1], [17, 17], c="k")
ax1.plot([0, 0], [15, 17], c="k")
ax1.plot([1, 1], [15, 17], c="k")
ax1.text(0.5, 17.5, "n.s.", ha="center")
ax1.set_ylim(-15, 20)
ax1.set_xticks([0, 1])
ax1.set_xticklabels(["SV\n(n=156)", "Random\n(n=200)"])
ax1.set_ylabel("FIMO RAG Heptamer Score")
ax1.set_title("RAG Heptamer Enrichment")

"""

Microhomology

"""
ax2 = plt.subplot(222)

sv_microhomology_data = []
random_microhomology_data = []

for line in open("Data\\SV_microhomology.txt"):
    sv_microhomology_data.append(int(line))
for line in open("Data\\Random_microhomology.txt"):
    random_microhomology_data.append(int(line))

max = np.amax([np.amax(sv_microhomology_data), np.amax(random_microhomology_data)])

sv_hist, sv_bin_edges = np.histogram(sv_microhomology_data, bins=max, range=(0, max))
random_hist, random_bin_edges = np.histogram(random_microhomology_data, bins=max, range=(0, max))

for i in range(0, len(sv_hist)):
    if sv_hist[i] > 0:
        ax2.add_patch(Rectangle((1 - sv_hist[i] / np.amax(sv_hist) / 2, sv_bin_edges[i]), sv_hist[i] / np.amax(sv_hist), 1, facecolor=sns.color_palette()[0], edgecolor=(63 / 255, 63 / 255, 63 / 255)))
    if random_hist[i] > 0:
        ax2.add_patch(Rectangle((3 - random_hist[i] / np.amax(random_hist) / 2, random_bin_edges[i]), random_hist[i] / np.amax(random_hist), 1, facecolor=sns.color_palette()[1], edgecolor=(63 / 255, 63 / 255, 63 / 255)))
print(sns.color_palette()[0])
ax2.set_xlim(0, 4)
ax2.set_ylim(0, 14)
ax2.plot([1, 3], [12, 12], c="k")
ax2.plot([1, 1], [11, 12], c="k")
ax2.plot([3, 3], [11, 12], c="k")
ax2.text(2, 12.5, "n.s.", ha="center")

ax2.set_xticks([1, 3])
ax2.set_xticklabels(["SV\n(n=" + str(len(sv_microhomology_data)) + ")", "Random\n(n=" + str(len(random_microhomology_data)) + ")"])
ax2.set_ylabel("Number of microhomology bases\naround breakpoint")
ax2.set_title("Microhomology")


"""

LCR

"""

ax3 = plt.subplot(223)
lcr_sv_hist_data = []

lcr_sv_scatter_data = {}
lcr_sv_type_data = {}
for line in open("Data\\SV_LCR_distance.txt"):
    lcr_sv_hist_data.append(int(line.split("\t")[4]))

    sv_id = line.split("\t")[3]
    dist = int(line.split("\t")[4])
    type = line.rstrip().split("\t")[5]

    if dist == 0:
        dist = 1

    if sv_id in lcr_sv_scatter_data:
        lcr_sv_scatter_data[sv_id].append(dist)
        lcr_sv_type_data[sv_id].append(type)
    else:
        lcr_sv_scatter_data[sv_id] = [dist]
        lcr_sv_type_data[sv_id] = [type]

x = []
y = []
for key in lcr_sv_scatter_data:
    if len(lcr_sv_scatter_data[key]) == 2:
        x.append(lcr_sv_scatter_data[key][0])
        y.append(lcr_sv_scatter_data[key][1])

        if lcr_sv_scatter_data[key][0] == 1 and lcr_sv_scatter_data[key][1] == 1:
            print(key)
            print(lcr_sv_type_data[key])


#a, b, c = ax3.hist(lcr_sv_hist_data, range=(0, 1000), bins=2000, cumulative=True)


size = 5

data_for_2d_hist = np.zeros((size, size))
for i_x in range(0, size):
    for i_y in range(0, size):
        for j in range(0, len(x)):
            if x[j] >= np.power(10, i_x) and x[j] < np.power(10, 1 + i_x) and y[j] >= np.power(10, i_y) and y[j] < np.power(10, 1 + i_y):
                data_for_2d_hist[i_x, i_y] += 1

for i_x in range(0, size):
    for i_y in range(0, size):
        color = 1 - data_for_2d_hist[i_x, i_y] / np.amax(data_for_2d_hist)
        ax3.add_patch(Rectangle((i_x, i_y), 1, 1, facecolor=(color, color, color)))

ax3.set_xlim(0, 5)
ax3.set_ylim(0, 5)
ax3.set_xticks([0, 1, 2, 3, 4, 5])
ax3.set_xticklabels(["0", "10", "100", "1000", "10000", "100000"])
ax3.set_yticks([0, 1, 2, 3, 4, 5])
ax3.set_yticklabels(["0", "10", "100", "1000", "10000", "100000"])



"""
lcr_random_hist_data = []
for line in open("Data\\Random_LCR_distance.txt"):
    lcr_random_hist_data.append(int(line.split("\t")[4]))
d, e, f = ax4.hist(lcr_random_hist_data, range=(0, 1000), bins=2000, cumulative=True)
"""


plt.show()
"""
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)

g = np.subtract(np.divide(a, len(lcr_sv_hist_data)), np.divide(d, len(lcr_random_hist_data)))
ax1.plot(range(0, len(g)), np.divide(a, len(lcr_sv_hist_data)))
ax2.plot(range(0, len(g)), np.divide(d, len(lcr_random_hist_data)))
ax3.plot(range(0, len(g)), g)
ax1.set_ylim(0, 1.1)
ax2.set_ylim(0, 1.1)
ax3.set_ylim(-0.505, 0.505)
plt.show()
"""