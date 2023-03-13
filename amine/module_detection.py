#!/usr/bin/env python3
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

Module containing the ModuleDetection class.
"""

from random import Random
from typing import Callable, Iterable
import statistics
import itertools
import networkx as nx
from scipy.stats import norm
from .models import Model


class ModuleDetection:
    """A class representing the module detection algorithm."""

    def __init__(
        self,
        G: nx.Graph,
        model: Model,
        fitness_function: Callable,
        background_correction: bool = True,
        precomputed_model: Model = None,
    ):
        """
        Class initialization.

        Parameters
        ----------
        G                     : nx.Graph
                                the graph
        model                 : Model
                                the model
        fitness_function      : Callable
                                the fitness function must accept two parameters
                                (a graph and a set of nodes) and return a float
        background_correction : bool, optional
                                specify if a background correction must be performed,
                                default is True
        precomputed_model     : Model, optional
                                specify a precomputed model to be used, default is None

        """
        self.rng = Random()
        self.G = G
        self.mod = model
        self.precomputed_model = precomputed_model
        self.fitness_function = fitness_function
        self.background_correction = background_correction
        self.bkg_scores_by_cluster_size = {}
        self.connected_components = None

    def get_bkg_score(self, size: int, ensure_connected: bool = True):
        """
        Calculate the mean and standard deviation of the score of random nodes.

        Parameters
        ----------
        size : int
               the number of node in the set
        ensure_connected : bool, optional
               indicates if the set of nodes must constitutes a connected module, default is True

        """
        try:
            return self.bkg_scores_by_cluster_size[size]
        except KeyError:
            scores_for_random_cluster = []
            # print(f"extractring background for clusters of size {size}")
            if ensure_connected:
                if not self.connected_components:
                    self.connected_components = [
                        self.G.subgraph(c).copy()
                        for c in nx.connected_components(self.G)
                    ]
                possible_seeds = []
                for component in self.connected_components:
                    if len(component.nodes) > size:
                        # print(f"component size={len(component.nodes)}")
                        # if len(component.nodes) > size:
                        possible_seeds.extend(component.nodes)

                # by uncommenting the code below, we restrict the seeds to be in the 10%
                # more expressed nodes, resulting in an increase of the mean background zscore
                # list_nodes = sorted([x for x in self.G.nodes],
                #                   key=lambda x: self.G.nodes[x]['weight'], reverse=True)
                # list_nodes = list_nodes[:len(list_nodes)//10]
                # possible_seeds = list(set(possible_seeds).intersection(list_nodes))

            for _ in range(10000):
                if ensure_connected:
                    seed = self.rng.choice(possible_seeds)
                    random_cluster = set([seed])
                    # seeds = set([seed])
                    neighbors = set([n for n in self.G[seed]])
                    while len(random_cluster) < size:
                        node_added = self.rng.choice(list(neighbors))
                        random_cluster.add(node_added)
                        if len(random_cluster) == size:
                            break
                        new_neighbors = {n for n in self.G[node_added]} - random_cluster
                        neighbors.update(new_neighbors)
                else:
                    random_cluster = self.rng.sample(self.G.nodes, size)
                scores_for_random_cluster.append(
                    self.fitness_function(self.G, random_cluster)
                )
            mean = statistics.mean(scores_for_random_cluster)
            dev = statistics.pstdev(scores_for_random_cluster)
            self.bkg_scores_by_cluster_size[size] = (mean, dev)
        return (mean, dev)

    def get_score(
        self,
        cluster: Iterable,
        background_correction: bool,
        ensure_connected_bgmodule=True,
    ):
        """
        Get the score of a cluster.

        Parameters
        ----------
        cluster                   : Iterable
                                    set of nodes
        background_correction     : bool
                                    specify if a background correction must be applied
        ensure_connected_bgmodule : bool, optional
                                    indicates if the set of background nodes must constitutes
                                    a connected module, default is True

        """
        score = self.fitness_function(self.G, cluster)
        pvalue = 0
        if background_correction:
            mean, dev = self.get_bkg_score(len(cluster), ensure_connected_bgmodule)
            score = (score - mean) / dev
            pvalue = 1 - norm.cdf(score)
        return score, pvalue

    def predict_modules(
        self,
        relevant_nodes: Iterable = None,
        max_nb_modules: int = None,
        cutoff: float = None,
    ) -> Iterable:
        """
        Predict a list of active module.

        Parameters
        ----------
        relevant_nodes : Iterable, optional
                         nodes on which research will be focused, default is None
        max_nb_modules : int, optional
                         maximum number of modules to return.
                         If None, returns all modules, default is None
        cutoff         : float, optional
                         cutoff used for pvalue, default is None

        """

        # memorization of the best clusters
        #
        # each solution is a tuple composed of:
        #     - the center node
        #     - the diameter (number of closest nodes)
        #     - the list of nodes in the cluster
        #     - the fitness
        #     - the pvalue of the fitness is background correction is used
        best = []

        # max and min size of a module
        max_size = self.G.graph["nb_nodes"] - 1
        min_size = 3

        relevant_nodes_id = set([])
        if relevant_nodes:
            relevant_nodes_id = set(
                [x for x in self.G.nodes if self.G.nodes[x]["label"] in relevant_nodes]
            )
        prioritary_nodes = sorted(
            [x for x in self.G.nodes if x in relevant_nodes_id],
            key=lambda x: self.G.nodes[x]["weight"],
            reverse=False,
        )
        other_nodes = sorted(
            [x for x in self.G.nodes if x not in relevant_nodes_id],
            key=lambda x: self.G.nodes[x]["weight"],
            reverse=False,
        )

        list_nodes = prioritary_nodes + other_nodes

        for node in list_nodes:
            similar = self.mod.get_most_similar(node, max_size - 1)
            clus = [node]
            best_for_this_node = [0, 0, None, 0, 1]

            for x in range(max_size - 1):
                try:
                    clus = clus + [similar[x]]
                except IndexError:
                    break
                if len(clus) < min_size:
                    continue
                score, pvalue = self.get_score(clus, background_correction=False)

                if score > best_for_this_node[3]:
                    best_for_this_node = [node, x + 1, clus, score, pvalue]
                else:
                    break

            if best_for_this_node == [0, 0, None, 0, 1]:
                continue

            best.append(best_for_this_node)

        # Sort the best modules by their scores
        best = sorted([x for x in best], key=lambda x: x[3], reverse=True)

        # enlarge the clusters and put the results in predicted
        predicted = []
        while best:
            predicted.append(best[0])
            best = best[1:]

            strategy = 1
            if strategy == 0:
                continue  # no merging

            # collect solutions centered on nodes belonging to the processed module
            while "the score of the module is increasing":
                # Collect in other_modules, the set of nodes belonging to other modules
                # whose center is in predicted
                other_modules = set(
                    [frozenset(x[2]) for x in best if x[0] in set(predicted[-1][2])]
                )

                # collect in 'candidates' the list of potential additions to the current module
                candidates = set()
                if strategy == 1:
                    # Union with other modules
                    candidates = other_modules
                elif strategy == 2:
                    # Union with a combination of other modules
                    for size in range(3):
                        for combination in itertools.combinations(
                            other_modules, size + 1
                        ):
                            new_clus = set(predicted[-1][2])
                            for solution in combination:
                                new_clus |= solution
                            candidates.add(frozenset(new_clus))
                elif strategy == 3:
                    # Union with neighborood
                    for x, _ in self.G.adj.items():
                        if x not in predicted[-1][2]:
                            continue
                elif strategy == 4:
                    # Union with the genes belonging to other modules
                    for x in other_modules:
                        candidates.add(frozenset(x))
                elif strategy == 5:
                    # Union with combination of at most 3 genes from other modules
                    for x in other_modules:
                        for size in range(3):
                            for combination in itertools.combinations(
                                x, min(len(x), size + 1)
                            ):
                                new_clus = set(predicted[-1][2]) | set(combination)
                                candidates.add(frozenset(new_clus))
                elif strategy == 6:
                    # Union with combination of the genes in other modules
                    for x in other_modules:
                        new_genes = set(predicted[-1][2]) - set(x)
                        for size in range(len(new_genes)):
                            for combination in itertools.combinations(
                                new_genes, size + 1
                            ):
                                new_clus = set(predicted[-1][2]) | set(combination)
                                candidates.add(frozenset(new_clus))

                best_combination = predicted[-1]

                for candidate in candidates:
                    if set(predicted[-1][2]).issuperset(candidate):
                        continue
                    new_clus = set(predicted[-1][2]) | candidate
                    score, pvalue = self.get_score(
                        new_clus, background_correction=False
                    )
                    if score > best_combination[3]:
                        best_combination = [0, 0, sorted(new_clus), score, pvalue]

                # put in new best all modules that doesn't intersect with the processed one
                try_another_combination = False
                if best_combination[3] > predicted[-1][3]:
                    predicted[-1] = best_combination
                    try_another_combination = True
                # recompute new best
                if relevant_nodes_id:
                    best = [
                        x
                        for x in best
                        if not (set(x[2]) & set(predicted[-1][2]))
                        or (set(x[2]) & relevant_nodes_id)
                    ]
                else:
                    best = [
                        x
                        for x in best
                        if not len(set(x[2]) & set(predicted[-1][2])) > 0
                    ]

                if not try_another_combination:
                    break

        sorted_modules = [
            x[2:] for x in sorted(predicted, key=lambda x: x[3], reverse=True)
        ]
        for x in sorted_modules:
            x[2] = self.get_score(x[0], background_correction=True)[1]

        # remove duplicates
        predicted_modules = []
        last = set()
        for x in sorted_modules:
            if set(x[0]) != last:
                predicted_modules.append([sorted(x[0]), x[1], x[2]])
                last = set(x[0])

        predicted_modules = [
            x for x in predicted_modules if set(x[0]) & relevant_nodes_id
        ] + [x for x in predicted_modules if not (set(x[0]) & relevant_nodes_id)]

        if max_nb_modules:
            predicted_modules = predicted_modules[:max_nb_modules]
        if cutoff:
            predicted_modules = [x for x in predicted_modules if x[2] <= cutoff]

        # output format for each predicated module
        #     1 - list of nodes (space separated)
        #     2 - zscore
        #     3 - pvalue
        return predicted_modules
