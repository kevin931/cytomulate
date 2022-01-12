# This is only a temporary file
# used for developing the whole workflow

# Read file and cell type
from utilities import FileIO
from cell_type import CellType

# data visual
import numpy as np
import matplotlib.pyplot as plot

# data analysis
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import networkx as nx
from networkx.algorithms import tree
from networkx.algorithms.community import greedy_modularity_communities

temp = FileIO()

expression_matrix = temp.load_data("./test_data/expression_matrix.txt")
channels = expression_matrix[0]
expression_matrix = expression_matrix[1]

# I think we should expect the users to have clustered
# the dataset
# Or we can write another method in maybe
# utilities that can cluster
labels = temp.load_data("./test_data/labels.txt", dtype = "str")
labels = labels[1]

unique_labels = np.unique(labels)

# We first initialize the CellType objects
# we also associate with eacy cell type
# their corresponding data
n_markers = np.shape(expression_matrix)[1]
cell_types = {}
cell_types_data = {}
id = 0
for c_type in unique_labels:
    cell_types[c_type] = CellType(id = id, name = c_type,
                                  n_markers = n_markers)
    ind = np.where(labels == c_type)[0]
    cell_types_data[c_type] = expression_matrix[ind, :]
    cell_types[c_type].overall_mean = np.mean(expression_matrix[ind, :], axis=0)
    cell_types[c_type].overall_cov = np.cov(expression_matrix[ind, :])
    id += 1

# If cell differentiation is sought after
# We will use the overall means of the each cell type
# to find pairwise distances
# Then we will construct a weighted complete graph
# Then construct an MST
# Then use a community detection to partition the trees
cell_differentiation = True
cell_means = np.zeros((len(unique_labels), n_markers))
for i in range(len(unique_labels)):
    cell_means[i,:] = cell_types[unique_labels[i]].overall_mean

complete_G = nx.Graph()
for i in range(len(unique_labels)):
    for j in range(i, len(unique_labels)):
        complete_G.add_edge(i,j, weight = np.linalg.norm(cell_means[i,:] - cell_means[j,:]))

mst = tree.minimum_spanning_edges(complete_G, algorithm="kruskal", data=False)

mst_G = nx.Graph()
mst_G.add_edges_from(list(mst))
forest = list(greedy_modularity_communities(mst_G))

root_list = []
for t in range(len(forest)):
    nodes_list = list(forest[t])
    root_list.append(np.random.choice(nodes_list))
    doing_list = [root_list[t]]
    done_list = []
    while len(doing_list) > 0:
        parent_cell = doing_list.pop(0)
        done_list.append(parent_cell)
        for e in list(mst_G.edges):
            if parent_cell in e:
                child_cell = (set(e) - {parent_cell}).pop()
                if (child_cell not in done_list) and (child_cell in nodes_list):
                    doing_list.append(child_cell)
                    cell_types[unique_labels[child_cell]].parent = parent_cell
                    cell_types[unique_labels[parent_cell]].children.append(child_cell)







background_noise_variance = np.Inf
max_components = 9
min_components = 1
cv_types = ["spherical", "tied", "diag", "full"]


c_type = "Basophils"
df = cell_types_data[c_type]
col_median = np.median(df, axis=0).reshape(-1,1)

# We will use K-means with 2 groups to
# group markers into two groups
# expressed and unexpressed
kmeans = KMeans(n_clusters=2).fit(col_median)

group_0_mean = np.mean(col_median[kmeans.labels_ == 0])
group_1_mean = np.mean(col_median[kmeans.labels_ == 1])
if group_1_mean > group_0_mean:
    expressed_group = 1
    unexpressed_group = 0
else:
    expressed_group = 0
    unexpressed_group = 1

expressed_markers = np.where(kmeans.labels_ == expressed_group)[0]
unexpressed_markers = np.where(kmeans.labels_ == unexpressed_group)[0]

# For unexpressed markers
# we fit gaussian mixture model with 2 components to each marginal
# as it seems reasonable to assume unexpressed markers are independent
# one with background noise
# one with lowly expressed protein
# The main purpose is to speed up the algorithm
# since the majority of the markers are un/lowly expressed
unexpressed_data = df[:, unexpressed_markers]
expressed_data = df[:, expressed_markers]

model_for_expressed_markers = {}
model_for_unexpressed_markers = {}

for m in range(len(unexpressed_markers)):
    gm = GaussianMixture(n_components=2).fit(unexpressed_data[:, m].reshape(-1, 1))
    est_background_noise_variance = gm.covariances_[np.argmin(gm.means_)][0][0]
    if est_background_noise_variance < background_noise_variance:
        background_noise_variance = est_background_noise_variance
    model_for_unexpressed_markers[unexpressed_markers[m]] = gm

# Since we expect only a few markers to be highly expressed
# we can probably afford fitting the entire data
# with multivariate GMM as well as some model selection
smallest_bic = np.Inf
current_bic = 0
best_gm = None
for n_components in range(min_components, max_components + 1):
    for cv_type in cv_types:
        gm = GaussianMixture(n_components=n_components,
                             covariance_type=cv_type).fit(expressed_data)
        current_bic = gm.bic(expressed_data)
        if current_bic < smallest_bic:
            smallest_bic = current_bic
            best_gm = gm

model_for_expressed_markers["all"] = best_gm








plot.scatter(unexpressed_data[:,0], expressed_data[:,0])

gm_labels = gm.predict(unexpressed_data[:,1].reshape(-1,1))

a = gm.sample(1000)[0]
b = gm1.sample(1000)[0]
plot.scatter(a,b)

plot.scatter(np.arange(946), unexpressed_data[:,1].reshape(-1,1),
             c = gm_labels, s= 1)












ind = np.where(labels == "Monocytes")[0]
y = expression_matrix[ind, 0]

plot.hist(y, bins = 100)
plot.show()

density = gaussian_kde(y)
x = np.linspace(-1,5,300)
density.covariance_factor = lambda : .25
density._compute_covariance()
plot.plot(x, density(x))
plot.show()

fig, ax = plot.subplots()
ax.scatter(np.arange(len(ind[0])), expression_matrix[ind, 0],
           s = 0.1, c = "black")
plot.show()


