# Probably a bad idea...
import warnings
from typing import Union

import numpy as np
import pandas as pd
import tqdm
from scipy.sparse import csgraph

warnings.filterwarnings("ignore", category=RuntimeWarning) 

gen = np.random.default_rng(1337)

def small_world_propensity(
    W: Union[np.ndarray, list], bin: Union[bool, list] = False
) -> pd.DataFrame:
    """Find the small-world propensity and related measures of a network W.

    Args:
        W (Union[np.ndarray, list]): Either an adjacency matrix or a list of adjacency matrices.
        bin (Union[bool, list], optional): Indicated whether the matrix (list of matrices) is (are) binary.
            Defaults to False.

    Returns:
        pd.DataFrame: Dataframe containing the small-world propensity and related measures.
    """
    # Check if W is a list of matrices
    if isinstance(W, list):
        df = pd.DataFrame()
        for i in tqdm.tqdm(range(len(W))):
            df_temp = _small_world_propensity(W[i], bin=bin[i])
            df = pd.concat([df, df_temp], ignore_index=True)
        return df
    else:
        df = _small_world_propensity(W, bin=bin)
        return df


def get_avg_rad_eff(W: np.ndarray) -> int:
    n = len(W)
    num_con = len(np.where(W > 0)[0])
    avg_deg_unw = num_con / n
    avg_rad_unw = avg_deg_unw / 2
    avg_rad_eff = np.ceil(avg_rad_unw)

    return int(avg_rad_eff)


def get_average_paths(W: np.ndarray) -> float:
    """Get the average path length of a network.

    Args:
        W (np.ndarray): Adjacency matrix of the network.

    Returns:
        float: Average path length of the network.
    """
    path_matrix = csgraph.shortest_path(1/W, directed=False, unweighted=False)
    L_W = np.triu(path_matrix).sum() / (len(W) * (len(W) - 1) / 2)

    return L_W


def get_clustering_coefficient(W: np.ndarray) -> float:
    """Get the clustering coefficient of a network.

    Args:
        W (np.ndarray): Adjacency matrix of the network.

    Returns:
        float: Clustering coefficient of the network.
    """
    K = np.where(W > 0, 1, 0).sum(axis=1)
    W2 = W / W.max()
    cyc3 = np.diagonal(np.linalg.matrix_power(W2 ** (1/3), 3))
    K = np.where(cyc3 == 0, np.inf, K)
    C = cyc3 / (K * K-1)

    return C.mean()


def _small_world_propensity(W: np.ndarray, bin: bool = False) -> pd.DataFrame:
    """Finds the small-world propensity and related measures of a single network W.

    Args:
        W (np.ndarray): Adjacency matrix of the network.
        bin (bool, optional): Is the matrix binary or not. Defaults to False.

    Returns:
        pd.DataFrame: Dataframe containing the small-world propensity and related measures.
    """

    # Check if the matrix is symmetric
    if not np.allclose(W, W.T, rtol=1e-05, atol=1e-08):
        W = make_symmetric(W, bin=bin)

    W = W / np.max(W)

    avg_rad_eff = get_avg_rad_eff(W)

    # Make regular and random networks
    W_reg = regular_matrix_generator(W, int(avg_rad_eff))
    W_rand = randomize_matrix(W)

    # Clustering
    C_W = get_clustering_coefficient(W)
    C_reg = get_clustering_coefficient(W_reg)
    C_rand = get_clustering_coefficient(W_rand)

    # Path lengths
    L_W = get_average_paths(W)
    L_reg = get_average_paths(W_reg)
    L_rand = get_average_paths(W_rand)

    # Delta L
    A = L_W - L_rand
    if A < 0:
        A = 0

    diff_path = A / (L_reg - L_rand)

    if np.isinf(L_W) or np.isinf(L_reg) or np.isinf(L_rand):
        diff_path = 1

    if diff_path > 1:
        diff_path = 1

    delta_L = diff_path

    # Delta C
    B = C_reg - C_W
    if B < 0:
        B = 0

    diff_clus = B / (C_reg - C_rand)
    if np.isnan(C_reg) or np.isnan(C_W) or np.isnan(C_rand):
        diff_clus = 1

    if diff_clus > 1:
        diff_clus = 1

    delta_C = diff_clus
    
    # Small-world propensity
    SWP = 1 - (np.sqrt((delta_C) ** 2 + (delta_L) ** 2) / np.sqrt(2))

    alpha = np.arctan(delta_L / delta_C)
    delta = (4 * alpha / np.pi) - 1

    df = pd.DataFrame(
        {
            "Network C": C_W,
            "Network L": L_W,
            "ΔC": delta_C,
            "ΔL": delta_L,
            "SWP": SWP,
            "α": alpha,
            "δ": delta,
            "Regular C": C_reg,
            "Random C": C_rand,
            "Regular L": L_reg,
            "Random L": L_rand,
        },
        index=[0],
    )

    return df


def randomize_matrix(A: np.ndarray) -> np.ndarray:
    """Randomly rewire the edges of a network.

    Args:
        A (np.ndarray): Adjacency matrix of the network.

    Returns:
        np.ndarray: Adjacency matrix of the randomized network.
    """
    num_nodes = A.shape[0]
    A_rand = np.zeros((num_nodes, num_nodes))
    mask = np.triu(np.ones((num_nodes, num_nodes)), 1)

    # Find the indices where mask > 0 in column-major order
    grab_indices = np.column_stack(np.nonzero(mask.T))

    # Access A with the indices
    orig_edges = A[grab_indices[:, 0], grab_indices[:, 1]]
    num_edges = len(orig_edges)
    rand_index = np.random.choice(num_edges, num_edges, replace=False)
    randomized_edges = orig_edges[rand_index]
    edge = 0
    for i in range(num_nodes - 1):
        for j in range(i + 1, num_nodes):
            A_rand[i, j] = randomized_edges[edge]
            A_rand[j, i] = randomized_edges[edge]
            edge += 1  # Move to next edge
    return A_rand


def regular_matrix_generator(G: np.ndarray, r: int) -> np.ndarray:
    """Generate a regular matrix from a given matrix.

    Args:
        G (np.ndarray): Adjacency matrix of the network.
        r (int): The average effective radius of the network.

    Returns:
        np.ndarray: Adjacency matrix of the regularized network.
    """

    n = len(G)
    G = np.triu(G)
    B = G.flatten(order="F")
    B = np.sort(B)[::-1]
    num_els = np.ceil(len(B) / (2 * n))
    num_zeros = 2 * n * num_els - n * n
    B = np.concatenate((B, np.zeros(int(num_zeros))))
    B = B.reshape((n, -1), order="F")

    M = np.zeros((n, n))

    for i in range(n):
        for z in range(r):
            a = gen.integers(0, n)
            while (B[a, z] == 0 and z != r - 1) or (
                B[a, z] == 0 and z == r - 1 and len(B[:, r - 1].nonzero()[0]) != 0
            ):
                a = gen.integers(0, n)

            y_coord = (i + z + 1) % n
            M[i, y_coord] = B[a, z]
            M[y_coord, i] = B[a, z]

            B[a, z] = 0

    return M


def make_symmetric(A: np.ndarray, bin: bool = False) -> np.ndarray:
    """Take an adjacency matrix and make it symmetric.

    Args:
        A (np.ndarray): Adjacency matrix of the network.
        bin (bool, optional): Is the network binary not. Defaults to False.

    Returns:
        np.ndarray: Symmetric adjacency matrix of the network.
    """
    W = np.zeros(A.shape)
    if bin:
        for i in range(A.shape[0]):
            for j in range(i, A.shape[1]):
                if A[i, j] or A[j, i]:
                    W[i, j] = 1
                    W[j, i] = 1
    else:
        for i in range(A.shape[0]):
            for j in range(i, A.shape[1]):
                if A[i, j] or A[j, i]:
                    val = (A[i, j] + A[j, i]) / 2
                    W[i, j] = val
                    W[j, i] = val

    return W
