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

Script used to import data
"""

import urllib.request
import progressbar
import yaml
from .parameters import Param


class ProgressBar:
    """
    Class encapsulating all variables and methods needed for the progress bar
    """

    def __init__(self, prefix):
        self.pbar = None
        self.widgets_max_size = [
            prefix,
            " ",
            progressbar.Bar(marker="#", left="[", right="]"),
            " ",
            progressbar.Percentage(),
            " | ",
            progressbar.DataSize(),
            " | ",
            progressbar.FileTransferSpeed(),
        ]
        self.widgets_unknown_size = [
            prefix,
            " ",
            progressbar.BouncingBar(),
            " ",
            progressbar.DataSize(),
            " | ",
            progressbar.FileTransferSpeed(),
        ]

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            if total_size > 0:
                self.pbar = progressbar.ProgressBar(
                    maxval=total_size, widgets=self.widgets_max_size
                )
            else:
                self.pbar = progressbar.ProgressBar(
                    maxval=progressbar.UnknownLength, widgets=self.widgets_unknown_size
                )
            self.pbar.start()

        downloaded = block_num * block_size
        if total_size < 0 or (total_size > 0 and downloaded < total_size):
            self.pbar.update(downloaded)
        else:
            self.pbar.finish()


# Load the program parameters from the 'config.yaml' file if it exists
try:
    with open("config.yaml", "r") as f:
        params = yaml.safe_load(f)
    for key, value in params.items():
        setattr(Param, key, value)
    print("Parameters read from 'config.yaml' file.")
except yaml.YAMLError:
    print("Problem occured when reading 'config.yaml'.")
    print("Using default parameters.")
except FileNotFoundError:
    pass

print("Downloading files")

# String
print(f"Data from String database version {Param.string_version}")
for specie_name, _ in Param.specie2id.items():
    local_filename, local_path = Param.get_string_ppi_file(specie_name)
    distant_name = f"https://stringdb-static.org/download/protein.links.detailed.v{Param.string_version}/{local_filename}"
    local_name, headers = urllib.request.urlretrieve(
        distant_name, local_path, ProgressBar(local_filename)
    )
    local_filename, local_path = Param.get_string_info_file(specie_name)
    distant_name = f"https://stringdb-static.org/download/protein.info.v{Param.string_version}/{local_filename}"
    local_name, headers = urllib.request.urlretrieve(
        distant_name, local_path, ProgressBar(local_filename)
    )
    local_filename, local_path = Param.get_string_alias_file(specie_name)
    distant_name = f"https://stringdb-static.org/download/protein.aliases.v{Param.string_version}/{local_filename}"
    local_name, headers = urllib.request.urlretrieve(
        distant_name, local_path, ProgressBar(local_filename)
    )

# BioGRID
print(f"Data from BIOGRID database version {Param.biogrid_version}")
local_filename, local_path = Param.get_biogrid_file()
distant_name = f"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-{Param.biogrid_version}/{local_filename}"
pbar = ProgressBar(local_filename)
local_name, headers = urllib.request.urlretrieve(distant_name, local_path, pbar)
pbar.pbar.finish()

# Intact
print("Data from IntAct database")
local_filename, local_path = Param.get_intact_file()
distant_name = (
    f"https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/{local_filename}"
)
local_name, headers = urllib.request.urlretrieve(
    distant_name, local_path, ProgressBar(local_filename)
)
