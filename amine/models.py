#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
+-------------------------------------------------------------------------------------+
| This file is part of AMINE                                                          |
|                                                                                     |
| AMINE is free software: you can redistribute it and/or modify it under the terms of |
| the GNU General Public License as published by the Free Software Foundation, either |
| version 3 of the License, or (at your option) any later version.                    |
| You should have received a copy of the GNU General Public License along with AMINE. |
| If not, see <http://www.gnu.org/licenses/>.                                         |
|                                                                                     |
| Author: Claude Pasquier (I3S Laboratory, CNRS, Université Côte d'Azur)              |
| Contact: claude.pasquier@univ-cotedazur.fr                                          |
| Created on decembre 20, 2022                                                        |
+-------------------------------------------------------------------------------------+

Various graph models.

Each class defined here is instanciated with a networkx graph
After instanciation, it is possible to get, according to the model,
the closest nodes to a specific node with the method 'get_most_similar'
"""

import pathlib
from abc import ABC, abstractmethod
from typing import Union, Iterable

import networkx as nx
import numpy as np
from gensim.models import FastText, Word2Vec
from scipy.spatial import distance

from .dimension_reduction import node2vec
from .dimension_reduction.pecanpy import node2vec as n2v


class Model(ABC):
    """
    abstract class.

    Three implementation are proposed:
        - Node2vec,
        - RandomWalk
        - SVD.
    """

    @abstractmethod
    def get_most_similar(self, elt: str, number: int):
        """
        Collect similar nodes ; method implemented by each subclass.

        Parameters
        ----------
        elt    : str
                 node identifier
        number : int
                 number of most similar nodes to collect

        """


class Node2vec(Model):
    """Node2vec model."""

    def __init__(self):
        """Declare variables."""
        self.model = None
        self.num_walks = 20  # 10
        self.walk_length = 100  # 80
        self.directed = False
        self.param_p = 1  # 4  # 0.15
        self.param_q = 1  # 2
        self.dimensions = 64  # 128
        self.window_size = 5  # 10
        self.workers = 4
        self.epoch = 10  # 10

    def init(self, G: nx.Graph, list_nodes: Iterable = None, precomputed: str = None):
        """
        Initialize the model with a weighted graph.

        Parameters
        ----------
        G           : nx.Graph
                      the graph used to initialize the model
        list_nodes  : Iterable, optional
                      specify an order of the nodes to be used, default is None
        precomputed : str or file-like object, optional
                      None or path to the precomputed model that must be used, default is None

        """
        if list_nodes is None:
            list_nodes = list(G.nodes)
        if precomputed and pathlib.Path(precomputed).is_file():
            self.model = Word2Vec.load(precomputed)
        else:
            for node in G.nodes():
                for nbr in sorted(G.neighbors(node)):
                    G[node][nbr]["weight"] = 1 - abs(
                        G.nodes[node]["weight"] - G.nodes[nbr]["weight"]
                    )
            self.compute_embedding(G, list_nodes)
            if precomputed:
                self.save(precomputed)

    def compute_embedding(self, G: nx.Graph, list_nodes: list):
        """
        Compute embedding.

        Parameters
        ----------
        G           : nx.Graph
                      the processed graph
        list_nodes  : list of nodes
                      the list of start nodes from the randomwalk

        """
        use_pecanpy = False
        if use_pecanpy:
            # from pecanpy import node2vec as n2v
            graph = n2v.SparseOTF(
                p=self.param_p,
                q=self.param_q,
                workers=self.workers,
                verbose=False,
                extend=True,
            )
            A = np.array(
                nx.adjacency_matrix(
                    G, nodelist=sorted(G.nodes), weight="weight"
                ).todense(),
                dtype=np.float_,
            )
            # isolated_nodes = np.where(~A.any(axis=1))[0]
            # print(np.where(~A.any(axis=0))[0])
            # print(nx.is_connected(G))
            # A = np.delete(A, isolated_nodes, axis=0)
            # A = np.delete(A, isolated_nodes, axis=1)
            graph.from_mat(A, sorted(G.nodes))
            walks = graph.simulate_walks(
                num_walks=self.num_walks,
                walk_length=self.walk_length,
                list_nodes=list_nodes,
            )
        else:
            graph = node2vec.Graph(G, self.directed, self.param_p, self.param_q)
            graph.preprocess_transition_probs()
            walks = graph.simulate_walks(
                self.num_walks, self.walk_length, nodes=list_nodes
            )

        # Learn embeddings by optimizing the Skipgram objective using SGD.
        walks = [list(map(str, walk)) for walk in walks]
        # import pickle
        # with open("/home/cpasquie/Téléchargements/test.txt", "wb") as fp:   #Pickling
        #     pickle.dump(walks, fp)
        # dd
        # with open("/home/cpasquie/Téléchargements/test.txt", "rb") as fp:   # Unpickling
        #     walks = pickle.load(fp)
        use_fasttext = False
        if use_fasttext:
            self.model = FastText(
                vector_size=self.dimensions,
                window=self.window_size,
                min_count=1,
                sentences=walks,
                epochs=self.epoch,
                max_n=0,
                sg=1,
            )
        else:
            self.model = Word2Vec(
                walks,
                vector_size=self.dimensions,  # size=self.dimensions,
                window=self.window_size,
                min_count=5,
                negative=5,
                sg=1,
                workers=self.workers,
                epochs=self.epoch,
            )  # iter=self.epoch)

    def save(self, fname_or_handle: str):
        """
            Save the model to file.

        Parameters
        ----------
        fname_or_handle : str or file-like object
                          path or handle to file where the model will be persisted

        """
        self.model.save(fname_or_handle)

    def load(self, fname_or_handle: str):
        """
        Load a previously saved model from a file.

        Parameters
        ----------
        fname_or_handle : str or file-like object
                          path or handle to file that contains the model

        """
        self.model = Word2Vec.load(fname_or_handle)

    def get_most_similar(self, elt: str, number: int):
        """
        Collect similar nodes.

        Parameters
        ----------
        elt    : str
                 node identifier
        number : int
                 number of most similar nodes to collect

        """
        return [int(x[0]) for x in self.model.wv.similar_by_word(str(elt), topn=number)]

    def get_distance(self, elt1: str, elt2: str):
        """
        Return the distance between two elements.

        Parameters
        ----------
        elt1 : str
            first element
        elt2 : str
            second element

        """
        return self.model.wv.distance(str(elt1), str(elt2))

    def get_vector(self, elt: Union[str, int]):
        """
        Get the vector encoding the element

        Parameters
        ----------
        elt : Union[str, int]
            the element

        Returns
        -------
        vector
            the vector encoding the element
        """
        return self.model.wv.get_vector(str(elt))


class RandomWalk(Model):
    """
    RandomWalk model.
    """

    # convergence criterion - when vector L1 norm drops below 10^(-6)
    # (this is the same as the original RWR paper)
    conv_threshold = 0.000001

    def __init__(self):
        """Declare variables."""
        self.nodelist = None
        self.walk_length = 0
        self.restart_prob = 0
        self.T = None

    def init(self, G: nx.Graph):
        """
        Initialize the model with a weighted graph.

        Parameters
        ----------
        G   : nx.Graph
              the graph used to initialize the model

        """
        self.nodelist = sorted(G.nodes)
        self.walk_length = 200
        self.restart_prob = 0.7

        for node in G.nodes():
            for nbr in sorted(G.neighbors(node)):
                G[node][nbr]["weight"] = 1 - abs(
                    G.nodes[node]["weight"] - G.nodes[nbr]["weight"]
                )

        # Create the adjacency matrix of G
        A = np.array(
            nx.adjacency_matrix(G, nodelist=self.nodelist, weight="weight").todense(),
            dtype=np.float_,
        )
        # Create the degree matrix
        D = np.diag(np.sum(A, axis=0))

        # The Laplacian matrix L, not used here is equal to D - A

        # Compute the inverse of D
        # Several solutions are possible
        #     - first solution: numpy.inverse
        #       inverse_of_d = numpy.linalg.inv(D)
        #     - second solution: numpy.solve
        #       inverse_of_d = numpy.linalg.solve(D, numpy.identity(len(nodes_list))
        #     - third solution, as the matrix is diagonal, one can use
        #       the inverse of the diagonal values
        #       inverse_of_d = np.diag(1 / np.diag(D))

        inverse_of_d = np.diag(1 / np.diag(D))

        # compute the transition matrix
        self.T = np.dot(inverse_of_d, A)

    def get_most_similar(self, elt: str, number: int):
        """
        Collect similar nodes.

        Parameters
        ----------
        elt    : str
                 node identifier
        number : int
                 number of most similar nodes to collect

        """
        arr = [0] * len(self.nodelist)
        arr[elt] = 1
        p_0 = np.array(arr)
        state_matrix = p_0
        for _ in range(self.walk_length):

            # evaluate the next state vector
            p_1 = (1 - self.restart_prob) * np.dot(
                state_matrix, self.T
            ) + self.restart_prob * p_0

            # calculate L1 norm of difference between p^(t + 1) and p^(t),
            # for checking the convergence condition
            diff_norm = np.linalg.norm(np.subtract(p_1, state_matrix), 1)
            if diff_norm < RandomWalk.conv_threshold:
                break
        state_matrix = p_1
        result = sorted(
            enumerate(state_matrix.tolist()), key=lambda res: res[1], reverse=True
        )
        return [int(x[0]) for x in result][1 : number + 1]


class SVD(Model):
    """SVD model."""

    def __init__(self):
        """Declare variables."""
        self.nodelist = None
        self.most_similar = []

    def init(self, G: nx.Graph):
        """
        Initialize the model with a weighted graph.

        Parameters
        ----------
        G   : nx.Graph
              the graph used to initialize the model

        """
        self.nodelist = sorted(G.nodes)
        A = np.array(
            nx.adjacency_matrix(G, sorted(G.nodes), weight=None).todense(),
            dtype=np.float_,
        )
        U, S, _ = np.linalg.svd(A, full_matrices=False)
        reduced_dimension = A.shape[0] // 5
        reduced_matrix = U * S
        reduced_matrix = reduced_matrix[:, 0:reduced_dimension]
        self.most_similar = []
        for ctr in range(reduced_matrix.shape[0]):
            dist = distance.cdist(
                reduced_matrix[ctr : ctr + 1], reduced_matrix[0:], "cosine"
            )
            self.most_similar.append(
                [
                    x[0]
                    for x in sorted(
                        list(enumerate(dist[0].tolist())), key=lambda x: x[1]
                    )
                ][1:]
            )

    def get_most_similar(self, elt: str, number: int):
        """
        Collect similar nodes.

        Parameters
        ----------
        elt    : str
                 node identifier
        number : int
                 number of most similar nodes to collect

        """
        return self.most_similar[elt][:number]
