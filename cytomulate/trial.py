# This is only a temporary file
# used for developing the whole workflow

# Read file and cell type
from utilities import FileIO
from cell_type import CellType

# data visual
import numpy as np
import matplotlib.pyplot as plot

# data analysis
from collections import Counter
from scipy.stats import gaussian_kde
from scipy.optimize import minimize
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import networkx as nx
from networkx.algorithms import tree
from networkx.algorithms.community import greedy_modularity_communities

from typing import Union, Optional, Any, List, Callable

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
frequency = Counter(labels)
observed_p = np.zeros(len(unique_labels))
counter = 0
for l in unique_labels:
    observed_p[counter] = frequency[l]/len(labels)
    counter += 1



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
    cell_types[c_type].overall_cov = np.cov(expression_matrix[ind, :], rowvar=False)
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
tree_list = []
root_list = []
for t in range(len(forest)):
    tree_list.append(nx.DiGraph())
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
                    tree_list[t].add_edge(parent_cell, child_cell)




list(nx.topological_sort(tree_list[0]))
# If we have bead information, we can use those to estimate background noise variance
# If not, we will use the cells
background_noise_variance = np.Inf
max_components = 9
min_components = 1
cv_types = ["full"]

for c_type in cell_types:
    #c_type = "Mature_B_cells"
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

    cell_types[c_type].expressed_markers = np.where(kmeans.labels_ == expressed_group)[0]
    cell_types[c_type].unexpressed_markers = np.where(kmeans.labels_ == unexpressed_group)[0]

    # For unexpressed markers
    # we fit gaussian mixture model with 2 components to each marginal
    # as it seems reasonable to assume unexpressed markers are independent
    # one with background noise
    # one with lowly expressed protein
    # The main purpose is to speed up the algorithm
    # since the majority of the markers are un/lowly expressed
    unexpressed_data = df[:, cell_types[c_type].unexpressed_markers]
    expressed_data = df[:, cell_types[c_type].expressed_markers]

    cell_types[c_type].model_for_expressed_markers = {}
    cell_types[c_type].model_for_unexpressed_markers = {}

    for m in range(len(cell_types[c_type].unexpressed_markers)):
        gm = GaussianMixture(n_components=2).fit(unexpressed_data[:, m].reshape(-1, 1))
        est_background_noise_variance = gm.covariances_[np.argmin(gm.means_)][0][0]
        if est_background_noise_variance < background_noise_variance:
            background_noise_variance = est_background_noise_variance
        cell_types[c_type].model_for_unexpressed_markers[cell_types[c_type].unexpressed_markers[m]] = gm

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

    cell_types[c_type].model_for_expressed_markers["all"] = best_gm




def moment_of_beta(alpha, beta, moment_type):
    if moment_type == "mean":
        moment = alpha / (alpha + beta)
    elif moment_type == "variance":
        moment = (alpha * beta) / (np.power(alpha + beta, 2) * (alpha + beta + 1))
    elif moment_type == "second":
        moment = (alpha * beta) / (np.power(alpha + beta, 2) * (alpha + beta + 1)) + \
                 np.power(alpha / (alpha + beta) , 2)
    else:
        raise ValueError('Unknown moment type')

    return moment

def moment_of_multinomial(p, n, moment_type):
    d = len(p)
    p = p.reshape(-1)
    if moment_type == "mean":
        moment = n * p
    elif moment_type == "variance":
        moment = np.zeros((d,d))
        for i in range(d):
            for j in range(d):
                if i == j:
                    moment[i,j] = p[i] * (1 - p[i])
                else:
                    moment[i,j] = -p[i] * p[j]
        moment = n * moment
    elif moment_type == "second":
        v = np.zeros((d,d))
        for i in range(d):
            for j in range(d):
                if i == j:
                    v[i,j] = p[i] * (1 - p[i])
                else:
                    v[i,j] = -p[i] * p[j]
        m2 = np.reshape(n * p, (-1, 1)) @ np.reshape(n * p, (1, -1))

        moment = v + m2
    else:
        raise ValueError('Unknown moment type')

    return moment



def rearrange_mean(expressed_mean_parameters,
                   unexpressed_mean_parameters,
                   expressed_markers,
                   unexpressed_markers):
    n_expressed = len(expressed_markers)
    n_unexpressed = len(unexpressed_markers)
    n_markers = n_expressed + n_unexpressed

    mean_parameters = np.zeros((n_markers))
    mean_parameters[expressed_markers] = expressed_mean_parameters
    mean_parameters[unexpressed_markers] = unexpressed_mean_parameters
    return mean_parameters


def rearrange_covariance(expressed_covariance_parameters,
                         unexpressed_covariance_parameters,
                         expressed_markers,
                         unexpressed_markers):
    n_expressed = len(expressed_markers)
    n_unexpressed = len(unexpressed_markers)
    n_markers = n_expressed + n_unexpressed

    covariance_parameters = np.zeros((n_markers, n_markers))
    covariance_parameters[np.ix_(expressed_markers, expressed_markers)] = expressed_covariance_parameters
    covariance_parameters[np.ix_(unexpressed_markers, unexpressed_markers)] = unexpressed_covariance_parameters

    return covariance_parameters


def parameters_to_mean_and_covariance(parameters,
                                      expressed_markers,
                                      unexpressed_markers,
                                      children_cell_types):
    n_expressed = len(expressed_markers)
    n_unexpressed = len(unexpressed_markers)
    n_markers = n_expressed + n_unexpressed


    mean_parameters = parameters[0:n_markers]
    covariance_parameters = parameters[n_markers : (len(parameters) - len(children_cell_types) * n_markers)]
    pseudo_time_parameters = parameters[(len(parameters) - len(children_cell_types) * n_markers) : ]

    expressed_mean = mean_parameters[0:len(expressed_markers)]
    unexpressed_mean = mean_parameters[(len(expressed_markers)) : ]
    expressed_mean = np.exp(expressed_mean)
    unexpressed_mean = np.exp(unexpressed_mean)

    mean = rearrange_mean(expressed_mean,
                          unexpressed_mean,
                          expressed_markers,
                          unexpressed_markers)

    temp_basis = covariance_parameters[0:np.power(n_expressed,2)]
    temp_diag = covariance_parameters[np.power(n_expressed,2) : (np.power(n_expressed,2) + n_expressed)]
    temp_basis = temp_basis.reshape((n_expressed, n_expressed))
    Q, R = np.linalg.qr(temp_basis)
    temp_diag = np.exp(temp_diag)
    temp_diag = np.diag(temp_diag)

    expressed_covariance = np.matmul(np.matmul(Q, temp_diag), Q.transpose())

    unexpressed_covariance = covariance_parameters[(np.power(n_expressed,2) + n_expressed) : ]
    unexpressed_covariance = np.exp(unexpressed_covariance)
    unexpressed_covariance = np.diag(unexpressed_covariance)

    covariance = rearrange_covariance(expressed_covariance,
                                      unexpressed_covariance,
                                      expressed_markers,
                                      unexpressed_markers)

    pseudo_time = 1/(1 + np.exp(-pseudo_time_parameters))

    pseudo_time = np.reshape(pseudo_time, (-1, len(children_cell_types)), order='F')

    return mean, covariance, pseudo_time



def unconditional_mean_and_covariance(mean_parameters,
                                      covariance_parameters,
                                      pseudo_time_parameters,
                                      background_noise_variance,
                                      children_cell_types):
    #
    slopes_of_paths_to_children = np.zeros((len(mean_parameters), len(children_cell_types)))
    children_cell_types_ids = np.zeros(len(children_cell_types))
    p = np.ones((len(children_cell_types))) / len(children_cell_types)
    counter = 0
    for c_type in children_cell_types:
        children_cell_types_ids[counter] = c_type.id
        slopes_of_paths_to_children[:, counter] = c_type.overall_mean - mean_parameters
        counter += 1

    pseudo_time_mean = np.zeros((len(mean_parameters), len(children_cell_types)))
    pseudo_time_covariance = np.zeros((len(mean_parameters), len(children_cell_types)))

    for i in range(len(mean_parameters)):
        for j in range(len(children_cell_types)):
            pseudo_time_mean[i, j] = moment_of_beta(pseudo_time_parameters[i, j], 1, "mean")
            pseudo_time_covariance[i, j] = moment_of_beta(pseudo_time_parameters[i, j], 1, "variance")

    differentiation_path_conditional_mean = np.zeros((len(mean_parameters), len(children_cell_types)))
    differentiation_path_conditional_covariance = np.zeros((len(mean_parameters), len(children_cell_types)))
    for j in range(len(children_cell_types)):
        differentiation_path_conditional_mean[:, j] = slopes_of_paths_to_children[:, j] * \
                                                     pseudo_time_mean[:, j]
        differentiation_path_conditional_covariance[:, j] = np.power(slopes_of_paths_to_children[:, j], 2) * \
                                                           pseudo_time_covariance[:, j]

    Omega_mean = moment_of_multinomial(p, 1, "mean").reshape((-1,1))
    Omega_covariance = moment_of_multinomial(p, 1, "variance")
    Omega_second = moment_of_multinomial(p, 1, "second")

    differentiation_path_unconditional_mean = differentiation_path_conditional_mean @ Omega_mean
    differentiation_path_unconditional_covariance = differentiation_path_conditional_mean @ Omega_covariance @ \
                                                    differentiation_path_conditional_mean.transpose()
    for j in range(len(children_cell_types)):
        differentiation_path_unconditional_covariance += Omega_second[j, j] * np.diag(differentiation_path_conditional_covariance[:, j])

    mean = mean_parameters + differentiation_path_unconditional_mean.reshape((-1))

    covariance = covariance_parameters + differentiation_path_unconditional_covariance + background_noise_variance * np.eye(len(mean_parameters))

    return mean, covariance


def distance_between_observed_and_adjusted(parameters,
                                           expressed_markers,
                                           unexpressed_markers,
                                           background_noise_variance,
                                           observed_mean,
                                           observed_covariance,
                                           children_cell_types):
    mean, covariance, pseudo_time = parameters_to_mean_and_covariance(parameters,
                                                                      expressed_markers,
                                                                      unexpressed_markers,
                                                                      children_cell_types)
    mean, covariance = unconditional_mean_and_covariance(mean,
                                                         covariance,
                                                         pseudo_time,
                                                         background_noise_variance,
                                                         children_cell_types)

    dist = np.linalg.norm(mean - observed_mean) + np.linalg.norm(covariance - observed_covariance)

    return dist

def adjust_mean_and_covariance(parameters0,
                               expressed_markers,
                               unexpressed_markers,
                               background_noise_variance,
                               observed_mean,
                               observed_covariance,
                               children_cell_types):
    res = minimize(distance_between_observed_and_adjusted,
                   x0 = parameters0,
                   args = (expressed_markers,
                           unexpressed_markers,
                           background_noise_variance,
                           observed_mean,
                           observed_covariance,
                           children_cell_types),
                   method='nelder-mead',
                   options={'xatol': 1e-8, 'disp': True})

    return res



def adjust_gaussian_mixture_means(cell_type,
                                  target_mean):
    target_mean = target_mean.reshape((-1))
    expressed_markers = cell_type.expressed_markers
    unexpressed_markers = cell_type.unexpressed_markers
    expressed_mean = target_mean[expressed_markers]

    # We first adjust means for the expressed markers
    for c in range(cell_type.model_for_expressed_markers["all"].n_components):
        cell_type.model_for_expressed_markers["all"].means_[c,:] += \
            cell_type.model_for_expressed_markers["all"].weights_[c] * \
            (expressed_mean - cell_type.model_for_expressed_markers["all"].means_[c,:])

    # Then we adjust means for the unexpressed markers
    for m in unexpressed_markers:
        max_ind = np.argmax(cell_type.model_for_unexpressed_markers[m].means_)
        min_ind = -1 * max_ind + 1
        cell_type.model_for_unexpressed_markers[m].means_[min_ind, :] = 0
        cell_type.model_for_unexpressed_markers[m].means_[max_ind, :] = \
            target_mean[m] / cell_type.model_for_unexpressed_markers[m].weights[max_ind]

    return cell_type



def adjust_gaussian_mixture_covariance(cell_type,
                                       target_covariance):

    expressed_markers = cell_type.expressed_markers
    unexpressed_markers = cell_type.unexpressed_markers

    temp_covariance = target_covariance + cell_type.overall_mean.reshape((-1, 1)) @ cell_type.overall_mean.reshape((1, -1))
    expressed_cov = temp_covariance[np.ix_(expressed_markers, expressed_markers)]
    unexpressed_cov = np.diag(temp_covariance[np.ix_(unexpressed_markers, unexpressed_markers)])

    # For expressed markers
    n_components = cell_type.model_for_expressed_markers["all"].n_components
    for i in range(n_components):
        expressed_cov -= cell_type.model_for_expressed_markers["all"].weights_[i] * \
                         (cell_type.model_for_expressed_markers["all"].means_[i,:].reshape((-1, 1)) @
                          cell_type.model_for_expressed_markers["all"].means_[i,:].reshape((1, -1)))

    for i in range(n_components):
        cell_type.model_for_expressed_markers["all"].covariances_[i, :, :] += cell_type.model_for_expressed_markers["all"].weights_[i] * \
                                                                             (expressed_cov - cell_type.model_for_expressed_markers["all"].covariances_[i, :, :])
    # For unexpressed markers
    counter = 0
    for m in unexpressed_markers:
        max_ind = np.argmax(cell_type.model_for_unexpressed_markers[m].means_)
        min_ind = -1 * max_ind + 1
        unexpressed_cov[counter] -= cell_type.model_for_unexpressed_markers[m].weights_[max_ind] * \
                                    cell_type.model_for_unexpressed_markers[m].means_[max_ind]
        cell_type.model_for_unexpressed_markers[0].covariances_[min_ind, :] = 0
        cell_type.model_for_unexpressed_markers[0].covariances_[max_ind, :] = unexpressed_cov[counter] / \
                                                                              cell_type.model_for_unexpressed_markers[m].weights[max_ind]
        counter += 1

    return cell_type





tree0 = tree_list[0]
topological_sorting_0 = list(reversed(list(nx.topological_sort(tree0))))

for c_id in topological_sorting_0:
    c_type = unique_labels[c_id]
    successor_ids = list(tree0.successors(c_id))
    if len(successor_ids) == 0:
        # This is a leaf node
        continue
    else:
        children_cell_types = []
        for s in successor_ids:
            children_cell_types.append(cell_types[unique_labels[s]])

            parameters = np.array([])

            expressed_markers = cell_types[c_type].expressed_markers
            unexpressed_markers = cell_types[c_type].unexpressed_markers

            observed_mean = cell_types[c_type].overall_mean
            observed_covariance = cell_types[c_type].overall_cov

            log_mean = np.log(observed_mean)
            parameters = np.concatenate((parameters,
                                         log_mean[expressed_markers],
                                         log_mean[unexpressed_markers]))

            cov = observed_covariance
            expressed_cov = cov[np.ix_(expressed_markers, expressed_markers)]
            w, v = np.linalg.eig(expressed_cov)
            parameters = np.concatenate((parameters,
                                         v.reshape((-1)),
                                         np.log(w)))

            unexpressed_cov = cov[np.ix_(unexpressed_markers, unexpressed_markers)]
            parameters = np.concatenate((parameters,
                                         np.log(np.diag(unexpressed_cov))))

            parameters = np.concatenate((parameters,
                                         np.ones(len(successor_ids) * n_markers) * (0.4)))

            res = adjust_mean_and_covariance(parameters,
                                             expressed_markers,
                                             unexpressed_markers,
                                             background_noise_variance,
                                             observed_mean,
                                             observed_covariance,
                                             children_cell_types)

            mean, covariance, pseudo_time = parameters_to_mean_and_covariance(res.x,
                                                                              expressed_markers,
                                                                              unexpressed_markers,
                                                                              children_cell_types)
            cell_types[c_type].overall_mean = mean
            cell_types[c_type].overall_cov = covariance

            cell_types[c_type] = adjust_gaussian_mixture_means(cell_types[c_type],
                                                               mean)
            cell_types[c_type] = adjust_gaussian_mixture_covariance(cell_types[c_type],
                                                                    covariance)


def linear_function(start_value: Union[int, float], end_value: Union[int, float]) -> Callable[
    [Union[int, float]], Union[int, float]]:
    """ Generate a linear function

    :param start_value: the starting value
    :param end_value: the ending value
    :return: a linear function (line segment) that interpolate the two points
    """

    def line_segment(t: Union[int, float]) -> Union[int, float]:
        return start_value + t * (end_value - start_value)

    return line_segment


n_cells = 5000

simulation_data = np.zeros((n_cells, n_markers))

cell_type_indices = np.random.choice(unique_labels, size=n_cells,
                                     replace=True, p = observed_p)

for n in range(n_cells):
    c_type = cell_type_indices[n]
    cell_type = cell_types[c_type]

    x = np.zeros(n_markers)
    x[cell_type.expressed_markers] = cell_type.model_for_expressed_markers["all"].sample(1)[0][0]
    for m in cell_type.unexpressed_markers:
        x[m] = cell_type.model_for_unexpressed_markers[m].sample(1)[0][0]

    g = np.zeros(n_markers)
    children_cell_types = cell_type.children
    if len(children_cell_types) > 0:
        child_id = np.random.choice(children_cell_types, size = 1,
                                    p = np.ones(len(children_cell_types))/len(children_cell_types))
        child_cell = cell_types[unique_labels[child_id]]

        for m in range(n_markers):
            g[m] = linear_function(cell_type.overall_mean[m], child_cell.overall_mean[m])(np.random.beta(0.4,1,1))

    E = np.random.multivariate_normal(np.zeros(n_markers), np.eye(n_markers), size = 1) * np.sqrt(background_noise_variance)

    y = x + g + E

    simulation_data[n,:] = y








