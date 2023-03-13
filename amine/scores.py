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

Various scoring methods

A first group of functions computes the score of a set of nodes on a graph,
optionnaly specifying the attribute that is used to store the nodes' weights

A second group of functions measures the accuracy of a prediction given
a predicted set of nodes and a set of nodes considered as ground truth
"""

import math
import statistics
from typing import Set

import networkx as nx
from scipy.stats import hypergeom, norm
from sklearn import metrics


class Scores:
    """The class grouping all scoring methods"""

    @staticmethod
    def modularity(G: nx.Graph, cluster: Set) -> float:
        """
        Compute the modularity of a set of nodes

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes

        Returns
        -------
        float
            the modularity of the set of nodes
        """

        edge_inside = 0
        total_degree = 0
        modularity = 0
        for n in cluster:
            neigh = set(list(G.neighbors(n)))
            total_degree += len(neigh)
            inside = len(neigh & cluster)
            edge_inside += inside
        modularity += (edge_inside / 2 / G.graph["nb_edges"]) - (
            total_degree / 2 / G.graph["nb_edges"]
        ) ** 2
        return modularity

    @staticmethod
    def clustering(G: nx.Graph, cluster: Set) -> float:
        """
        Compute the clustering coefficient of a set of nodes that is equal to
        the sum of the clustering coefficient of the nodes of the cluster

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes

        Returns
        -------
        float
            the clustering coefficient
        """
        cl = nx.clustering(G, nodes=cluster)
        return sum(list(cl.values())) / len(cluster)

    @staticmethod
    def density(G: nx.Graph, cluster: Set) -> float:
        """
        Return the density of the subgraph composed of nodes belonging to the cluster

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes

        Returns
        -------
        float
            the density of the cluster
        """
        subgraph = nx.subgraph(G, cluster)
        return nx.density(subgraph)

    @staticmethod
    def normalized_cut(G: nx.Graph, cluster: Set) -> float:
        """
        Compute the normalized cut of the cluster
        the formulae is described in:
        David Gleich. Hierarchical Directed Spectral Graph Partitioning.
        <https://www.cs.purdue.edu/homes/dgleich/publications/Gleich%202005%20-%20hierarchical%20directed%20spectral.pdf>

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes

        Returns
        -------
        float
            the normalized cut
        """
        return nx.normalized_cut_size(G, cluster, G.nodes - cluster)

    @staticmethod
    def average_shortest_path(G: nx.Graph, cluster: Set) -> float:
        """
        compute average shortest path in a module composed of vertices belonging to cluster

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes

        Returns
        -------
        float
            the average shortest path
        """
        sum_shortest_path_length = 0
        for n1 in cluster:
            for n2 in cluster:
                if n1 == n2:
                    continue
                sum_shortest_path_length += nx.shortest_path_length(
                    G, source=n1, target=n2
                )
        nb_nodes = len(cluster)
        return sum_shortest_path_length / (nb_nodes * (nb_nodes - 1))

    @staticmethod
    def aggregation_from_pvalue(G: nx.Graph, cluster: Set, attr: str) -> float:
        """
        compute the Z score of a module with Stouffer's Z method used in:
        Li, D., Pan, Z., Hu, G., Zhu, Z., & He, S. (2017).
        Active module identification in intracellular networks using a memetic algorithm
        with a new binary decoding scheme.
        BMC Genomics, 18(S2), 209. https://doi.org/10.1186/s12864-017-3495-y

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes
        attr : str
            the name of the attribute storing the weigh of nodes

        Returns
        -------
        float
            the Z score of a cluster based on p-values
        """
        if "zscores" not in G.graph:
            scores = {}
            for n in G.nodes:
                if attr in G.nodes[n]:
                    scores[n] = min(10, norm.ppf(1 - G.nodes[n][attr]))
                    # if norm.ppf returns inf, 10 is used
            G.graph["zscores"] = scores
        score = 0
        len_cluster = 0
        for n in cluster:
            score += G.graph["zscores"][n]
            len_cluster += 1
        return score / math.sqrt(len_cluster)

    @staticmethod
    def aggregation_from_normdist(G: nx.Graph, cluster: Set, attr: str) -> float:
        """
        compute the Z score of a module with Stouffer's Z method used in:
        Li, D., Pan, Z., Hu, G., Zhu, Z., & He, S. (2017).
        Active module identification in intracellular networks using a memetic algorithm
        with a new binary decoding scheme.
        BMC Genomics, 18(S2), 209. https://doi.org/10.1186/s12864-017-3495-y

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes
        attr : str
            the name of the attribute storing the weigh of nodes

        Returns
        -------
        float
            the Z score of a cluster based on values with normal distribution
        """
        if "zscores" not in G.graph:
            values = []
            for n in G.nodes:
                if attr in G.nodes[n]:
                    values.append(G.nodes[n][attr])
            mean = statistics.mean(values)
            dev = statistics.pstdev(values)

            zscores = {}
            for n in G.nodes:
                if attr in G.nodes[n]:
                    zscores[n] = (G.nodes[n][attr] - mean) / dev
            G.graph["zscores"] = zscores
        score = 0
        len_cluster = 0
        for n in cluster:
            if n in G.graph["zscores"]:
                score += G.graph["zscores"][n]
                len_cluster += 1
        return score / math.sqrt(len_cluster)

    @staticmethod
    def hypergeometric_score(G: nx.Graph, cluster: Set, attr: str) -> float:
        """
        compute the hypergeometric score of a module with a method inspired by:
        Breitling, Rainer, Anna Amtmann, and Pawel Herzyk.
        Graph-based iterative Group Analysis enhances microarray interpretation.
        BMC bioinformatics 5.1 (2004): 100.

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        cluster : Set
            a set of nodes
        attr : str
            the name of the attribute storing the weigh of nodes

        Returns
        -------
        float
            the hypergeometric score of a cluster
        """
        if "ranking_scores" not in G.graph:
            scores = []
            for n in G.nodes:
                if attr in G.nodes[n]:
                    scores.append((n, G.nodes[n][attr]))
            scores.sort(key=lambda x: x[1], reverse=True)
            ranks = {y[0]: x + 1 for x, y in enumerate(scores)}
            G.graph["ranking_scores"] = ranks

        len_cluster = 0
        worse_rank = 0
        for n in cluster:
            if n in G.graph["ranking_scores"]:
                worse_rank = max(worse_rank, G.graph["ranking_scores"][n])
                len_cluster += 1
        return 1 - hypergeom.pmf(
            len_cluster, G.graph["nb_nodes"], worse_rank, len_cluster
        )

    @staticmethod
    def measure_f1(G: nx.Graph, real_cluster: Set, pred_cluster: Set) -> float:
        """
        Compute the F1 score of a predicted cluster

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        real_cluster : Set
            the set of nodes composing the real cluster
        pred_cluster : Set
            the set of nodes composing the predicted cluster

        Returns
        -------
        float
            the F1 score
        """
        prediction = [0] * G.graph["nb_nodes"]
        for x in pred_cluster:
            prediction[x] = 1
        real = [0] * G.graph["nb_nodes"]
        for x in real_cluster:
            real[x] = 1
        return metrics.f1_score(real, prediction)

    # faster way to compute recall
    # @staticmethod
    # def measure_recall(G :nx.Graph, real_cluster :Set, pred_cluster :Set) -> float:
    #     return len(pred_cluster & real_cluster) / len(real_cluster)

    @staticmethod
    def measure_recall(G: nx.Graph, real_cluster: Set, pred_cluster: Set) -> float:
        """
        Compute the recall of a predicted cluster

        Parameters
        ----------
        G : nx.Graph
            a networkx graph
        real_cluster : Set
            the set of nodes composing the real cluster
        pred_cluster : Set
            the set of nodes composing the predicted cluster

        Returns
        -------
        float
            the recall score
        """
        prediction = [0] * G.graph["nb_nodes"]
        for x in pred_cluster:
            prediction[x] = 1
        real = [0] * G.graph["nb_nodes"]
        for x in real_cluster:
            real[x] = 1
        return metrics.recall_score(real, prediction)
