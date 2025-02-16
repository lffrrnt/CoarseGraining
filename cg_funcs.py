import numpy as np
import scipy
from scipy import stats

def normalize(data):
    """
    This function finds the nonzero indices of the signal and 
    sets the average of the nonzero signal to 1 by preserving the 
    zeros. 
    --------------------------------------
    Inputs: Data matrix in shape (N, T) where T is the length of the 
    time series, and N is the number of cells/sensors.
    --------------------------------------
    Return: Normalized data matrix in shape (N, T).
    """
    indices = {}
    for i in range(data.shape[0]):
        indices[i] = np.nonzero(data[i,:])[0]
    norm_data = data.copy()
    for i in range(data.shape[0]):
        norm_data[i,indices[i]] = data[i,indices[i]] / np.mean(data[i,indices[i]]) 
    return norm_data


def sort_correlations(matrix):
    """
    This function sorts the pairwise correlations by excluding
    the upper triangle and the diagonal of the correlation matrix.
    -------------------------------------
    Inputs: Correlation matrix (N, N)
    -------------------------------------
    Return: Sorted indices array.
    """
    sorted_corr_indices = []
    sorted_indices = np.unravel_index(np.argsort(np.abs(matrix), axis=None), matrix.shape)
    for i, j in zip(sorted_indices[0], sorted_indices[1]):
        if (i < j) & (i != j): 
            sorted_corr_indices.append((i,j))
    return np.array(sorted_corr_indices)

def find_NK_k(number):
    """
    Takes the number of cells/sensors as input
    and return the number of clusters at each RG step and
    k is the last step that we can not coarse anymore.
    """
    NK = []
    steps = np.floor(np.log(number)/np.log(2)) 
    for i in range(int(steps)):
        clusters = np.floor(number / 2)
        NK.append(clusters)
        number = clusters
    return NK, int(steps)

def coarse_graining_step(data, corr):
    """
    Performs one coarse-graining iteration. First, the number of new variables is determined.
    We set Nan values in the correlation matrix to 0, those zeros are because of the silent cells.
    Then we sort the correlations and take the maximum pair as a new variable sum their activity.
    We do it until we pair everything. 
    -----------------------------------------------------------------------------------------
    Return: The coarse-grained variables and the indices we grouped in one step.
    """
    x = data.copy()
    NK = int(x.shape[0]/2)
    #x_ = np.zeros((NK, x.shape[1]))
    x_ = []
    original_indices = np.zeros((NK, 2))
    corr[np.where(np.isnan(corr)==True)] = 0
    sorted_pairs = sort_correlations(corr)
    for i in range(NK):
        max_corr = sorted_pairs[-1]
        original_indices[i] = max_corr
        #x_[i,:] = x[max_corr[0],:] + x[max_corr[1],:]
        x_.append(x[max_corr[0],:] + x[max_corr[1],:])
        sorted_pairs = np.delete(sorted_pairs, (np.where(sorted_pairs == max_corr[0])[0]), axis=0)
        sorted_pairs = np.delete(sorted_pairs, (np.where(sorted_pairs == max_corr[1])[0]), axis=0)

    return np.array(x_), original_indices

def find_coarsed_variables(new_variables, k, normalization=True):

    cg_variables = []
    cluster_indices = []

    cg_variables.append(new_variables)
    corr_matrix = np.corrcoef(new_variables)
    for i in range(k):
        new_variables, indices = coarse_graining_step(new_variables, corr_matrix)
        if normalization == True:
            new_variables = normalize(new_variables)
        print("Performed coarse graining:", i,"/",k)
        corr_matrix = np.corrcoef(new_variables)
        cg_variables.append(new_variables)
        cluster_indices.append(indices)
        
    return cg_variables, cluster_indices 

#Functions to find the in-cluster-indices through coarse-graining and eigenspectrum of their covariance matrices

def find_in_cluster_indices(indices, k):
    cluster_list = {}
    cluster_list[0] = dict(enumerate(indices[0]))
    for c in range(k-1):
        clusters = {}
        for idx, i in enumerate(indices[c+1]):
            clusters[idx] = np.concatenate((cluster_list[c][int(i[0])], cluster_list[c][int(i[1])]))
        cluster_list[c+1] = clusters
    return cluster_list

def covariance_in_clusters(in_cluster_indices, data):
    cluster_data = {}
    for i in in_cluster_indices.keys():
        temp = []
        for j in in_cluster_indices[i]:
            temp.append(data[:, int(j)])
            cluster_data[i] = temp
    cov_clusters = {}
    for i in cluster_data.keys():
        cov_clusters[i] = np.cov(cluster_data[i])

    return cov_clusters

def find_mean_eigenvalues(cov_matrices):
    mean_eigvals = {}
    mean_eigvals_sem = {}
    for i in cov_matrices.keys():
        temp = []
        for j in cov_matrices[i].keys():
            #temp.append(np.linalg.eig(cov_matrices[i][j])[0])
            temp.append(scipy.linalg.eigh(cov_matrices[i][j], eigvals_only=True))
        mean_eigvals[i] = np.array(temp).mean(axis=0)
        mean_eigvals_sem[i] = stats.sem(np.array(temp), axis=0)
    return mean_eigvals, mean_eigvals_sem
