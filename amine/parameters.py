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

    project_name = "amine"

    string_version = "11.5"
    biogrid_version = "4.4.216"

    # Datasets are located in 'datadir'
    datadir = pathlib.Path(pathlib.Path(__file__).parent.parent, "data")

    precomputed_models_dir = pathlib.Path(datadir, "models")

    specie2id = {
        "human": "9606",
        "mouse": "10090",
        "drosophila": "7227",
        "Caenorhabditis elegans": "6239",
    }

    go_file = os.path.join(datadir, "go-basic.obo")
    goslim_file = os.path.join(datadir, "goslim_generic.obo")
    gene2go_file = os.path.join(datadir, "gene2go")

    verbose = False

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
