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

Module containing the class Param.
"""

import os
import pathlib


class Param:
    """Group program parameters."""

    # name of the project
    project_name = "amine"

    # version of the STRING database
    string_version = "11.5"

    # version of the BioGRID database
    biogrid_version = "4.4.216"


    # Directory where the datasets are stored
    # by default it is the directory named data on the project root directory
    datadir = pathlib.Path(pathlib.Path(__file__).parent.parent, "data")

    precomputed_models_dir = pathlib.Path(datadir, "models")

    # Specification of all the species handled by the program with their associated id
    specie2id = {
        "human": "9606",
        "mouse": "10090",
        "drosophila": "7227",
        "Caenorhabditis elegans": "6239",
    }

    # Specification of the default specie to process
    # It must correspond to a key used in specie2id
    default_specie = "human"

    # Specification of the default network to use that can be
    # 'string', 'biogrid', 'intact' or the path to a network file in edge list format
    default_network = "string"

    #
    # Filters to be applied to the STRING database
    # this correspond to the different composants of the evidence score
    #
    ppi_string_graph ={
        "min_neighborhood": 0,   "min_fusion": 0, "min_cooccurence" : 0,
        "min_coexpression": 0, "min_experimental": 0, "min_database" : 0,
        "min_textmining": 0, "min_combined_score" : 0}

    #
    # Filter to be applied to the BioGRID database
    # There is only one choice, corresponding to the type of interaction
    # that can be 'physical', 'genetic' or None is no filter is used
    #
    ppi_biogrid_graph ={
        "inter_type": None}

    #
    # Filter to be applied to the Intact database
    # This correspond to the minimum confidence score
    #
    ppi_intact_graph = {
        "min_confidence": 0.0}

    #
    # Parameters and Filter to be used when reading a custom graph in edgelist format
    #
    ppi_custom_graph ={
        "path": None,
        "source_col": 0,
        "target_col": 1,
        "header": None # Row number to use as the column names    
        }

    verbose = False

    go_file = os.path.join(datadir, "go-basic.obo")
    goslim_file = os.path.join(datadir, "goslim_generic.obo")
    gene2go_file = os.path.join(datadir, "gene2go")

    @staticmethod
    def get_ncbi_background_genes(specie, extension: str = "txt"):
        """Get file containing the background genes for GO entichments."""
        return pathlib.Path(Param.datadir, "gene_result_" + specie + "." + extension)

    @staticmethod
    def get_string_ppi_file(specie: str):
        """Get path to interaction file."""
        filename = f"{Param.specie2id[specie]}.protein.links.detailed.v{Param.string_version}.txt.gz"
        return filename, os.path.join(Param.datadir, "real", "ppi", filename)

    @staticmethod
    def get_string_info_file(specie: str):
        """Get path to protein information file."""
        filename = (
            f"{Param.specie2id[specie]}.protein.info.v{Param.string_version}.txt.gz"
        )
        return filename, os.path.join(Param.datadir, "real", "ppi", filename)

    @staticmethod
    def get_string_alias_file(specie: str):
        """Get path to protein alias file."""
        filename = (
            f"{Param.specie2id[specie]}.protein.aliases.v{Param.string_version}.txt.gz"
        )
        return filename, os.path.join(Param.datadir, "real", "ppi", filename)

    @staticmethod
    def get_biogrid_file() -> pathlib.Path:
        """
        get the path to the biogrid network file

        Returns
        -------
        pathlib.Path
            the path to the network file
        """
        filename = f"BIOGRID-ALL-{Param.biogrid_version}.tab3.zip"
        return filename, os.path.join(Param.datadir, "real", "ppi", filename)

    @staticmethod
    def get_intact_file() -> pathlib.Path:
        """
        get the path to the intact network file

        Returns
        -------
        pathlib.Path
            the path to the network file
        """
        filename = "intact.zip"
        return filename, os.path.join(Param.datadir, "real", "ppi", filename)

    @staticmethod
    def get_node2vec_model(specie: str) -> pathlib.Path:
        """
        get the path to the model file

        Returns
        -------
        pathlib.Path
            the path to the model file
        """
        filename = f"{Param.specie2id[specie]}.{Param.string_version}.model"
        return os.path.join(Param.precomputed_models_dir, filename)
