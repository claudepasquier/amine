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

Parsing of command line arguments.
"""

import argparse
import pathlib


def parse_arguments(prog):
    """Parse arguments."""
    parser = argparse.ArgumentParser(prog=prog)
    parser.add_argument(
        "-x",
        "--expvalue",
        dest="expvalue_file",
        type=pathlib.Path,
        required=True,
        help="specifies the file with expression values",
    )
    parser.add_argument(
        "-g",
        "--gene_col",
        dest="gene_col",
        type=int,
        required=True,
        help=(
            "specifies the index (in expvalue_file) of the column with gene ids (first"
            " columns=0)"
        ),
    )
    parser.add_argument(
        "-l",
        "--log2fc_col",
        dest="log2fc_col",
        type=int,
        required=False,
        help=(
            "specifies the index (in expvalue_file) of the column with log2 fold change"
            " values (first columns=0)"
        ),
    )
    parser.add_argument(
        "-p",
        "--pvalue_col",
        dest="pvalue_col",
        type=int,
        required=False,
        help=(
            "specifies the index (in expvalue_file) of the column with pvalues (first"
            " columns=0)"
        ),
    )
    parser.add_argument(
        "-a",
        "--aggreg",
        dest="aggregation_method",
        type=str,
        required=False,
        choices=["pvalue", "log2fc", "pivalue", "count"],
        help="specifies the method used to aggregate scores",
    )
    parser.add_argument(
        "-s",
        "--specie",
        dest="specie",
        type=str,
        required=False,
        choices=["human", "mouse", "drosophila"],
        help="specifies the specie",
    )
    parser.add_argument(
        "-f",
        "--focus",
        dest="focus_genes_file",
        type=pathlib.Path,
        required=False,
        help="specifies the file with the list of important genes",
    )
    parser.add_argument(
        "-m",
        "--model",
        dest="precomputed",
        required=False,
        default=False,
        action="store_true",
        help="indicates if precomputed model must be used",
    )
    parser.add_argument(
        "-n",
        "--network",
        dest="network",
        type=str,
        required=False,
        help="""specifies the network to use.
                It can be a path to an edge file or a predefined protein-protein
                interaction network to choose from 'string', 'biogrid' and 'intact'""",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=False,
        type=pathlib.Path,
        help="specifies the name of the output file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        required=False,
        action="store_true",
        default=False,
        help="displays results on screen",
    )
    return parser.parse_args()
