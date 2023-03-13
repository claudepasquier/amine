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

Module creating and initializing various datasets.

The created graphs are networkx graphs with the following additional characteristics:
    - each node can have a weight which is specified with the attribute 'weight'
    - each graph can be associated with ground truth information that can be
      a clustering or a specific module. The attribute storing the ground truth is
      specified with the graph attribute 'groups'
"""

import os
import re
import random
from random import Random
from typing import Any, Dict, Set, Iterable, Union, List
from zipfile import ZipFile
import math
import statistics
import xlrd
from scipy.stats import norm

import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import truncnorm
from numpy.random import RandomState

from .parameters import Param
from .scores import Scores
from .graph_generation import gencat
from .amine_exception import AmineException


class Datasets:
    """Class for accessing datasets."""

    @staticmethod
    def get_ppi_string_graph(
        specie: str,
        min_neighborhood: int = 0,
        min_fusion: int = 0,
        min_cooccurence: int = 0,
        min_coexpression: int = 0,
        min_experimental: int = 0,
        min_database: int = 0,
        min_textmining: int = 0,
        min_combined_score: int = 0,
        weight: Union[float, str] = 1.0,
    ) -> nx.Graph:
        """
        Read a String PPI graph, applying various filters on the interaction scores.


        Parameters
        ----------
        specie             : str
                             the specie modeled by the interaction graph
        min_neighborhood   : int, optional
                             minimum neighborhood score, by default 0
        min_fusion         : int, optional
                             minimum fusion score, by default 0
        min_cooccurence    : int, optional
                             minimum cooccurence score, by default 0
        min_coexpression   : int, optional
                             minimum coexpression score, by default 0
        min_experimental   : int, optional
                             minimum experimental score, by default 0
        min_database       : int, optional
                             minimum database score, by default 0
        min_textmining     : int, optional
                             minimum textmining score, by default 0
        min_combined_score : int, optional
                             minimum combined score, by default 0
        weight             : Union[float, str], optional
                             the value of weight between nodes, specified as
                             a float or an attribute name, by default 1.0
        Returns
        -------
        nx.Graph
            the created graph as an nx.Graph instance

        """

        df = pd.read_csv(Param.get_string_ppi_file(specie)[1], sep=" ", header=0)
        df = df[df.neighborhood >= min_neighborhood]
        df = df[df.fusion >= min_fusion]
        df = df[df.cooccurence >= min_cooccurence]
        df = df[df.coexpression >= min_coexpression]
        df = df[df.experimental >= min_experimental]
        df = df[df.database >= min_database]
        df = df[df.textmining >= min_textmining]
        df = df[df.combined_score >= min_combined_score]
        columns = list(df.columns)[2:]
        if weight in columns:
            columns = list(set(columns) - set([weight]))
        df = df.drop(columns=columns)

        df2 = pd.read_csv(Param.get_string_info_file(specie)[1], sep="\t", header=0)
        df2 = df2.drop(columns=["protein_size", "annotation"])
        mapping = {
            str(x[0]): str(x[1]) for x in list(df2.itertuples(index=False, name=None))
        }
        df["protein1"] = df["protein1"].map(mapping)
        df["protein2"] = df["protein2"].map(mapping)
        if weight in list(df.columns):
            df["weight"] = df[weight].map(lambda x: x / 1000)
            del df[weight]
        else:
            s_length = len(df["protein1"])
            df["weight"] = pd.Series([weight] * s_length, index=df.index)
        print(df.head())
        elist = [
            (x[0], x[1], {"weight": x[2]})
            for x in list(df.itertuples(index=False, name=None))
        ]
        G = Datasets.create_graph(elist)
        return G

    @staticmethod
    def get_ppi_biogrid_graph(
        specie: str, inter_type: str = None, weight: float = 1.0
    ) -> nx.Graph:
        """
        Read a BIOGRID graph

        Parameters
        ----------
        specie      : str
                      the specie modeled by the interaction graph
        inter_type  : str, optional
                      the type of interaction that can be 'physical' or 'genetic'
                      None means both interaction types, by default None
        weight      : float, optional
                      the value of weight between nodes, by default 1.0
        Returns
        -------
        nx.Graph
            the created graph as an nx.Graph instance

        """

        df = pd.read_csv(
            Param.get_biogrid_file()[1],
            sep="\t",
            header=0,
            usecols=[
                "Official Symbol Interactor A",
                "Official Symbol Interactor B",
                "Organism ID Interactor A",
                "Organism ID Interactor B",
                "Experimental System Type",
            ],
        )

        # filter by specie
        df = df[df["Organism ID Interactor A"] == int(Param.specie2id[specie])]
        df = df[df["Organism ID Interactor B"] == int(Param.specie2id[specie])]

        # filter by interaction type
        if inter_type:
            df = df[df["Experimental System Type"] == inter_type]

        elist = [
            (str(x[0]), str(x[1]), {"weight": weight})
            for x in list(df.itertuples(index=False, name=None))
        ]
        return Datasets.create_graph(elist)

    @staticmethod
    def get_ppi_intact_graph(
        specie: str, min_confidence: float = 0.0, weight: float = 1.0
    ) -> nx.Graph:
        """
        Read a BIOGRID graph

        Parameters
        ----------
        specie             : str
                             the specie modeled by the interaction graph
        min_confidence     : float, optional
                             the minimum confidence value, by default 0.0
        weight             : float, optional
                             the value of weight between nodes, by default 1.0
        Returns
        -------
        nx.Graph
            the created graph as an nx.Graph instance

        """

        # regular expression for parsing confidence string
        reg_confidence = re.compile(r"intact-miscore:([\d\.]+)")

        def get_confidence(conf: str) -> float:
            try:
                return float(reg_confidence.search(conf).group(1))
            except ValueError:
                return -1.0

        # regular expression for parsing gene name
        reg_name = re.compile(r"uniprotkb:([a-zA-Z0-9\-_]+)\(gene name\)")

        def get_gene_name(aliases: str) -> float:
            try:
                return reg_name.search(aliases).group(1)
            except AttributeError:
                return ""

        zf = ZipFile(Param.get_intact_file()[1])
        df = pd.read_csv(
            zf.open("intact.txt"),
            sep="\t",
            header=0,
            usecols=[
                "Alias(es) interactor A",
                "Alias(es) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Confidence value(s)",
            ],
        )

        # filter by specie
        df = df[df["Taxid interactor A"].str.contains(Param.specie2id[specie])]
        df = df[df["Taxid interactor B"].str.contains(Param.specie2id[specie])]

        # extract uniprot IDs and confidence scores as float
        df["protein1"] = df["Alias(es) interactor A"].map(get_gene_name)
        df = df[df["protein1"] != ""]
        df["protein2"] = df["Alias(es) interactor B"].map(get_gene_name)
        df = df[df["protein2"] != ""]
        df["score"] = df["Confidence value(s)"].map(get_confidence)
        df = df[df["score"] >= min_confidence]
        df = df.drop(
            columns=[
                "Alias(es) interactor A",
                "Alias(es) interactor B",
                "Taxid interactor A",
                "Taxid interactor B",
                "Confidence value(s)",
            ]
        )

        elist = [
            (str(x[0]), str(x[1]), {"weight": weight})
            for x in list(df.itertuples(index=False, name=None))
        ]
        return Datasets.create_graph(elist)

    @staticmethod
    def get_custom_graph(
        path: str,
        source_col: int = 0,
        target_col: int = 1,
        weight_col: int = None,
        header: Union[int, List[int], None] = None,
    ) -> nx.Graph:
        """
        Read a custom graph

        Parameters
        ----------
        path       : str
                     the path of the interaction graph
        source_col : int, optional
                     the column index of the source node, by default 0
        target_col : int, optional
                     the column index of the target node, by default 1
        weight_col : int, optional
                     the index of the column specifying the weight, by default None
        header     : Union[int, List[int], None], optional
                     Row number(s) to use as the column names, by default None
        Returns
        -------
        nx.Graph
            the created graph as an nx.Graph instance

        Raises
        ------
        Exception
            Ff there is a problem reading the edgelist file
        """

        col_indexes = [source_col, target_col]
        col_names = ["source", "target"]
        if weight_col:
            col_indexes.append(weight_col)
            col_names.append("weight")

        # sort columns
        z = sorted(zip(col_indexes, col_names))
        col_indexes, col_names = zip(*z)

        try:
            df = pd.read_excel(
                path, usecols=col_indexes, header=header, names=col_names
            )
        except (xlrd.biffh.XLRDError, ValueError):
            try:
                df = pd.read_csv(
                    path,
                    usecols=col_indexes,
                    header=header,
                    sep=None,
                    engine="python",
                    names=col_names,
                )
            except pd.errors.ParserError as e:
                raise AmineException(f"Unable to open '{path}'") from e

        if weight_col:
            elist = [
                (x[0], x[1], {"weight": x[2]})
                for x in list(df.itertuples(index=False, name=None))
            ]
        else:
            elist = [
                (str(x[0]), str(x[1]), {"weight": 1})
                for x in list(df.itertuples(index=False, name=None))
            ]
        return Datasets.create_graph(elist)

    @staticmethod
    def set_nodes_value(
        G: nx.Graph,
        attr: str,
        node2value: dict,
        default_value: float = None,
        case_sensitive: bool = True,
    ) -> nx.Graph:
        """
        Set values to nodes.

        Parameters
        ----------
        G              : nx.Graph
                         the graph to weight
        attr           : str
                         the name of the attribute storing the value
        node2value     : dict
                         a dict mapping node ids to values
        default_value  : float, optional
                         a default value for nodes not in node2value, by default None
        case_sensitive : bool, optional
                         indicate if the mapping of node ids is case sensitive, by default True

        Returns
        -------
        nx.Graph
            the updated graph as an nx.Graph instance
        """
        nb_mapped = 0
        print(len(G.nodes))
        if not case_sensitive:
            node2value = {x.upper(): y for x, y in node2value.items()}
            node2id = {x: G.nodes[x]["label"].upper() for x in G.nodes}
        else:
            node2id = {x: G.nodes[x]["label"] for x in G.nodes}
        for node, label in node2id.items():
            if label in node2value:
                nb_mapped += 1
                G.nodes[node][attr] = node2value[label]
            else:
                if default_value:
                    G.nodes[node][attr] = default_value
                else:
                    pass
        if Param.verbose:
            print(f"{nb_mapped} values mapped to node network")

        return G

    @staticmethod
    def map_query_to_gene_name(query: Iterable, specie: str) -> Dict:
        """
        Replace gene node with node name.

        Parameters
        ----------
        query          : Iterable
                         query terms
        specie         : str
                         the specie modeled by the interaction graph

        Returns
        -------
        Dict
            the mapping between query terms and gene names

        """
        df = pd.read_csv(Param.get_string_alias_file(specie)[1], sep="\t", header=1)
        lst = list(df.itertuples(index=False, name=None))

        # Identification of source.
        id2source = {x[1]: x[2].split() for x in lst}
        sources = {}
        for node in query:
            try:
                for src in id2source[node]:
                    try:
                        sources[src] += 1
                    except KeyError:
                        sources[src] = 1
            except KeyError:
                pass
        # Selected source is the one that allows to map the largest number of ids
        source = sorted(
            [(x, y) for x, y in sources.items()], key=lambda x: x[1], reverse=True
        )[0]
        query2stringid = {x[1]: x[0] for x in lst if source[0] in x[2]}

        df = pd.read_csv(
            Param.get_string_info_file(specie)[1],
            sep="\t",
            header=0,
            names=["id", "name", "size", "annotation"],
        )
        df = df.astype({"name": "str"})
        lst = list(df.itertuples(index=False, name=None))

        stringid2name = {x[0]: x[1] for x in lst}
        mapping = {}
        for node in query:
            try:
                mapping[node] = stringid2name[query2stringid[node]]
            except KeyError:
                pass
        return mapping

    @staticmethod
    def get_guyon_graph(graph_nb: int) -> nx.Graph:
        """
        Get a Guyon's graph.

        Parameters
        ----------
        graph_nb : int
             the identifier of the graph

        Returns
        -------
        nx.Graph
            the graph as an nx.Graph instance

        """
        zip_file = ZipFile(os.path.join(Param.datadir, "synthetic", "guyon.zip"))
        df = pd.read_csv(
            zip_file.open("pp1_" + str(graph_nb) + ".csv"), sep=",", header=0
        )
        df = df.drop(df.columns[0], axis=1)
        A = df.values
        G = nx.from_numpy_matrix(A)
        df = pd.read_csv(zip_file.open("pvalues.csv"), sep=",", header=0)
        df = df.drop(df.columns[0], axis=1)
        for ctr in range(G.number_of_nodes()):
            G.nodes[ctr]["weight"] = df.iloc[ctr, graph_nb - 1]

        df = pd.read_csv(zip_file.open("truehits.csv"), sep=",", header=0)
        df = df.drop(df.columns[0], axis=1)
        for ctr in range(G.number_of_nodes()):
            G.nodes[ctr]["truehit"] = df.iloc[ctr, graph_nb - 1]

        return Datasets.init_graph(G, edge_weight=None, default_groups="truehit")

    @staticmethod
    def get_gencat_graph(
        nb_nodes: int,
        nb_edges: int,
        nb_active_modules: int,
        active_module_size: List[int],
        rng_seed: int = None,
    ) -> nx.Graph:
        """
        Get a graph generated by gencat method.

        Parameters
        ----------
        nb_nodes          : int
                            number of nodes in the graph
        nb_edges          : int
                            number of edges in the graph
        nb_active_modules : int
                            number of active modules
        active_module_size: List[int]
                            size of each module
        rng_seed          : int, optional
                            seed for the random number generator.
                            If None, then global default RNG is used, by default None

        Returns
        -------
        nx.Graph
            the graph as an nx.Graph instance

        """

        rng = np.random.default_rng(rng_seed)
        # parameters
        n = nb_nodes
        m = nb_edges
        d = 1  # number of attributes
        # at this time only one attribute is allowed
        k = nb_active_modules * 6  # number of clusters
        max_deg = math.sqrt(n)  # maximum degree

        # M, class preference mean, a k x k  matrix expressing the average of connection
        # proportions from the nodes in class l1 to the nodes in class l2
        # we initalize M with a proportion of 60% intraclass links and the remaining
        # 40% over the other classe
        remaining_prop = 0.4 / (k - 1)
        M = remaining_prop * np.ones(shape=(k, k))
        np.fill_diagonal(M, 0.6)

        # D, class preference deviation, a k x k  matrix express the deviation of the number
        # of connection proportions from the nodes in class l1 to the nodes in class l2
        # we initalize D with a deviation of 0.3 for the intraclass likns and 0.1 for other links
        D = 0.1 * np.ones(shape=(k, k))
        np.fill_diagonal(D, 0.3)

        # H, attribute-class correlation is an d x k array representing the correlation between
        # a class and attribute values of nodes uses a correlation of 0.5 for active modules
        # and 0.0 for others
        H = np.zeros(shape=(d, k))
        H[:, 0:nb_active_modules] = 0.7

        phi_c = 1  # puissance of the power law for class size distribution
        omega = 0.2  # deviation of normal distribution for attributes
        r = 50  # number of iterations for edge generation
        step = (100,)  # used for inverse transform sampling
        att_type = "normal"  # attribute distribution (normal or Bernoulli)
        woAP = False  # indicates removing inverse transform sampling
        woITS = False  # indicates removing adjusting proportion
        S, X, Label = gencat.gencat_active(
            n,
            m,
            k,
            d,
            max_deg,
            M,
            D,
            H,
            nb_active=nb_active_modules,
            size_active=active_module_size,
        )

        G = nx.Graph()
        for ctr1 in range(S.shape[0] - 2):
            for ctr2 in range(ctr1, S.shape[0] - 1):
                if S[ctr1, ctr2]:
                    G.add_edge(ctr1, ctr2)
        attributes = [x[0] for x in X]
        mean = statistics.mean(attributes)
        dev = statistics.pstdev(attributes)
        for ctr, _ in enumerate(G.nodes):
            if Label[ctr] < nb_active_modules:
                G.nodes[ctr]["truehit"] = 1
            else:
                G.nodes[ctr]["truehit"] = 0

            method = "gencat"
            if method == "gencat":
                G.nodes[ctr]["weight"] = (
                    1 - norm.cdf(abs((attributes[ctr] - mean) / dev))
                ) * 2
            elif method == "guyon":
                # Guyon method
                if Label[ctr] < nb_active_modules:
                    # assign a truncated normal distribution for each node in cluster
                    mean, std = 0, 0.05
                    lower, upper = (0 - mean) / std, (1 - mean) / std
                    G.nodes[ctr]["weight"] = truncnorm.rvs(
                        loc=mean, scale=std, a=lower, b=upper, random_state=rng
                    )
                else:
                    # assign a uniform distribution for all nodes
                    G.nodes[ctr]["weight"] = rng.uniform(0, 1)
            elif method == "batra":
                # Batra method
                if Label[ctr] < nb_active_modules:
                    case = rng.normal(5, 1, 100)
                    control = rng.normal(2, 1, 10)
                    G.nodes[ctr]["weight"] = scipy.stats.mannwhitneyu(case, control)[1]
                else:
                    case = rng.normal(0, 1, 100)
                    control = rng.normal(0, 1, 10)
                    G.nodes[ctr]["weight"] = scipy.stats.mannwhitneyu(case, control)[1]
            # print(f"{G.nodes[ctr]['weight']}  {G.nodes[ctr]['truehit']}")

        return Datasets.init_graph(G, edge_weight=None, default_groups="truehit")

    @staticmethod
    def get_scale_free_graph(
        nb_nodes: int,
        nb_initial_nodes: int,
        nb_modules: int,
        module_size: int,
        p_prob: float,
        q_prob: float,
        rng_seed: int = None,
    ) -> nx.Graph:
        """
        Get a scale free graph based on the extended_barabasi_albert_graph.

        The selection of true hits is defined in Cornish et al., 2014,
        The way to set node weights is defined in Robinson et al., 2017

        the key point is that the method strives to select non-overlapping gene clusters

        Parameters
        ----------
        nb_nodes         : int
                           number of nodes in the graph
        nb_initial_nodes : int
                           number of initial nodes used by barabasi method
        nb_modules       : int
                           number of modules
        module_size      : int
                           size of each module
        p_prob           : float
                           probability value for adding an edge between existing nodes. p + q < 1
        q_prob           : float
                           probability value of rewiring of existing edges. p + q < 1
        rng_seed         : int, optional
                           seed for the random number generator.
                           If None, then global default RNG is used, by default None

        Returns
        -------
        nx.Graph
            the graph as an nx.Graph instance

        """

        def neighbors_order(G: nx.Graph, start: Any, k: int) -> int:
            """
            Get the number of neighbors of a node.

            Parameters
            ----------
            G       : nx.Graph
                      the simulated Network
            start   : Any
                      node that we want to know the neighborhood
            k       : int
                      neighborhood order

            Returns
            -------
            int
                the neighborhood of "start" at order=k

            """

            nbrs = set([start])
            for _ in range(k):
                nbrs = set((nbr for n in nbrs for nbr in G[n]))
            return nbrs

        def knbrs(G: nx.Graph, start: Any, k: int) -> Set:
            """
            Get the number of neighbors of a node.

            Parameters
            ----------
            G : nx.Graph
                Simulated Network
            start : Any
                node that we want to know the neighborhood
            k : int
                Neighborhood Order

            Returns
            -------
            Set of nodes
                neighborhood of "start"

            """

            neighbohood = set()
            while k >= 1:
                nbrs = set([start])
                for _ in range(k):
                    nbrs = set((nbr for n in nbrs for nbr in G[n]))
                neighbohood.update(nbrs)
                k -= 1
            return neighbohood

        def get_seeds(G: nx.Graph, number_seeds: int, min_distance: int) -> List:
            """
            Return possible nodes used as seeds

            Parameters
            ----------
            G : nx.Graph
                Graph containing simulated network
            number_seeds : int
                Number of centers of groups
            min_distance : int
                distance that each seed should be from each other

            Returns
            -------
            list of nodes
                the seeds

            """
            seeds = []
            list_deg = []
            for degree in G.degree():
                list_deg.append(degree[1])
            ndegree = np.quantile(list_deg, 0.8)
            # The variable "select" stores the vertices that can be chosen
            selected = []
            for degree in G.degree():
                # Excludes vertices with a high degree of connectivity
                if degree[1] <= ndegree:
                    selected.append(degree[0])

            for _ in range(number_seeds):
                sel = []
                # Chooses a node randomly in the "selected" set
                rdn = np_rng.choice(selected, size=1)[0]
                sel.append(rdn)
                # check the neighborhood of the chosen node
                neighbors = knbrs(G, rdn, min_distance)
                # Removes the neighbors of the selected seed from the set "selected"
                selected = list(set(selected).difference(neighbors))
                # The function returns false if the set "selected" is empty, that is,
                # if the chosen minimum distance is very high
                if selected.__len__() == 0:
                    return False

                seeds = list(set().union(seeds, sel))

            return seeds

        if rng_seed is None:
            rng = Random()
            np_rng = RandomState()
        else:
            rng = Random(rng_seed)
            np_rng = RandomState(rng_seed)

        G = nx.extended_barabasi_albert_graph(
            nb_nodes, nb_initial_nodes, p_prob, q_prob, rng_seed
        )
        # connect disconnected components
        comp = [x for x in nx.connected_components(G)]
        for ctr in range(len(comp) - 1):
            G.add_edge(
                random.choice(list(comp[ctr])), random.choice(list(comp[ctr + 1]))
            )

        batra_method = False
        if batra_method:
            clusters = []
            k = 3
            selected_nodes = set()
            while len(clusters) < nb_modules:
                group = set([rng.choice(list(G))])
                selected_nodes.update(group)
                while len(group) < module_size:
                    best_avg = 1000
                    best_node = None
                    for node in G.nodes:
                        if node in selected_nodes:
                            continue
                        avg = Scores.average_shortest_path(G, group | set([node]))
                        if abs(avg - k) < abs(best_avg - k):
                            best_avg, best_node = avg, node
                    group.add(best_node)
                    selected_nodes.add(best_node)
                clusters.append(group)
        else:
            # Guyon method

            # We chose a distance 10 between the seeds
            seeds = False
            distance = 10
            # If the distance is too high the algorithm decreases it
            # and tries to identify the seeds again
            while isinstance(seeds, bool):
                seeds = get_seeds(G, nb_modules, distance)
                distance -= 1
            neigh = set()
            clusters = []
            # The variable "neigh" stores the neighbors of the nodes that compose each group
            # Neighbors as well as seeds can not be selected.
            neigh.update(seeds)
            module_size -= 1
            count = 1
            for node in seeds:
                group = set()
                conect = set()
                group.add(node)
                nb_selected = 0
                while nb_selected < module_size:
                    finish = False
                    order = 0
                    rnd = rng.random()
                    # checks to neighborhood order 4 according with the function 1/pow(10, dist)
                    if rnd > 0.1 and rnd <= 1:
                        order = 1
                    elif 0.01 < rnd <= 0.1:
                        order = 2
                    elif 0.001 < rnd <= 0.01:
                        order = 3
                    elif 0.0001 < rnd <= 0.001:
                        order = 4
                    # neighbors_order returns all neighbors of the node n of a specific order
                    population = neighbors_order(G, node, order)
                    # Removes neighbors and nodes from other groups already identified
                    # besides the nodes of this same group that have already been chosen.
                    population = population - neigh
                    # After the phase above if there are no more nodes in the population
                    # the algorithm set the population using the next neighborhood order
                    while population.__len__() == 0:
                        order += 1
                        population = neighbors_order(G, node, order)
                        population = population - neigh
                        if order > 10:
                            finish = True
                            break
                    if not finish:
                        # one node is chosen randomly in the population
                        selected = np_rng.choice(list(population), size=1)[0]
                        group.add(selected)
                        # The variable neigh stores the node so that it is no more selected
                        neigh.add(selected)
                        # The entire neighborhood of the node (up to order 2) is saved so as
                        # not to be selected by a next seed
                        conect.update(knbrs(G, selected, 2))
                        nb_selected += 1
                    else:
                        print("pb: no selected module")
                        break  # nothing selected

                neigh.update(conect)
                clusters.append(group)
                count += 1
        batra_method = False
        if batra_method:
            for node in G.nodes:
                case = np_rng.normal(0, 1, 100)
                control = np_rng.normal(0, 1, 10)
                pvalue = scipy.stats.mannwhitneyu(case, control)[1]
                G.nodes[node]["weight"] = pvalue
                G.nodes[node]["truehit"] = 0

            ctr = 0
            for cluster in clusters:
                ctr += 1
                for node in cluster:
                    case = np_rng.normal(5, 1, 100)
                    control = np_rng.normal(2, 1, 10)
                    pvalue = scipy.stats.mannwhitneyu(case, control)[1]

                    G.nodes[node]["weight"] = pvalue
                    G.nodes[node]["truehit"] = ctr

        else:
            # Guyon method
            # assign a uniform distribution for all nodes
            for node in G.nodes:
                G.nodes[node]["weight"] = rng.uniform(0, 1)
                G.nodes[node]["truehit"] = 0

            # assign a truncated normal distribution for each node in cluster
            mean, std = 0, 0.05
            lower, upper = (0 - mean) / std, (1 - mean) / std
            ctr = 0
            for cluster in clusters:
                ctr += 1
                for node in cluster:
                    G.nodes[node]["weight"] = truncnorm.rvs(
                        loc=mean, scale=std, a=lower, b=upper, random_state=np_rng
                    )
                    G.nodes[node]["truehit"] = ctr

        return Datasets.init_graph(G, edge_weight=None, default_groups="truehit")

    @staticmethod
    def create_graph(elist: Iterable) -> nx.Graph:
        """
        Create a graph from edges list, removing selfloops

        Parameters
        ----------
        elist : Iterable
            list of edges

        Returns
        -------
        nx.Graph
            the created graph
        """
        G = nx.Graph(elist)
        G.remove_edges_from(nx.selfloop_edges(G))
        return Datasets.init_graph(G)

    @staticmethod
    def init_graph(
        G: nx.Graph,
        name: str = "",
        node_weight: str = "unchanged",
        edge_weight: str = "unchanged",
        default_groups: str = None,
    ) -> nx.Graph:
        """
        Initialize a graph

        Parameters
        ----------
        G : nx.Graph
            a graph
        name : str, optional
            the name of the graph, by default ""
        node_weight : str, optional
            specify the weight of nodes, by default 'unchanged'
            (None to set uniform weights equals to 1)
        edge_weight : str, optional
            specify the weight of edges, by default 'unchanged'
            (None to set uniform weights equals to 1)
        default_groups : str, optional
            the node attribute representing the known group(s), by default None

        Returns
        -------
        nx.Graph
            the initialized graph
        """

        # node labels must consist on consecutive integers starting at 0
        # if it is not the case, a mapping is necessary

        # convert nodes label to int if all labels are numbers and are
        # represented as strings

        # convert string labels to int values if possible
        try:
            if all([x.isdigit() for x in G.nodes]):
                mapping = {x: int(x) for x in G.nodes()}
                G = nx.relabel_nodes(G, mapping, copy=False)
        except AttributeError:
            pass

        mapping_needed = False

        if all(isinstance(n, int) for n in G.nodes()):
            if sorted([x for x in G.nodes()]) == list(range(G.number_of_nodes())):
                pass  # perfect
            else:
                mapping_needed = True
        else:
            mapping_needed = True

        if mapping_needed:
            nodes_mapping = dict(enumerate(sorted(G.nodes)))
            mapping = {y: x for x, y in nodes_mapping.items()}
            G = nx.relabel_nodes(G, mapping)
            for node in G.nodes:
                G.nodes[node]["label"] = nodes_mapping[node]
        for node in G.nodes:
            if "weight" not in G.nodes[node] or node_weight is None:
                G.nodes[node]["weight"] = 1
        for edge in G.edges:
            if "weight" not in G.edges[edge] or edge_weight is None:
                G.edges[edge]["weight"] = 1

        G.graph["default_groups"] = default_groups
        G.graph["nb_nodes"] = G.number_of_nodes()
        G.graph["nb_edges"] = G.number_of_edges()
        G.graph["name"] = name
        return G

    @staticmethod
    def get_groups(G: nx.Graph, groups_attribute: str = None) -> Dict[int, Set[int]]:
        """
        Return the groups identified in the graph

        Parameters
        ----------
        G : nx.Graph
            a graph
        groups_attribute : str, optional
            the node attribute representing the known group(s), by default None

        Returns
        -------
        Dict[int, Set[int]]
            a dict with the known groups

        """
        if not groups_attribute:
            groups_attribute = G.graph["default_groups"]
        if not groups_attribute:
            return {}
        groups = {}  # type: Map[int]
        for node in G.nodes:
            group_nb = int(G.nodes[node][groups_attribute])
            if group_nb == 0:
                continue
            try:
                groups[group_nb].add(node)
            except KeyError:
                groups[group_nb] = set([node])
        return groups
