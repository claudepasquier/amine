"""
    Entry point

    Example of call (from the project root directory):
    > python amine -g 0 -l 2 -p 6 -s mouse \
    -x ./data/real/expression/chiou_2017/Hmga2_positive_vs_negative.csv    
"""
import pathlib
import sys
import yaml

sys.path.append(str(pathlib.Path(__file__).parent.parent))

from amine.process_real_network import process
from amine.arg_parser import parse_arguments
from amine.parameters import Param


arg = parse_arguments("amine")
Param.verbose = arg.verbose

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

modules = process(
    arg.expvalue_file,
    arg.network if arg.network else Param.default_network,
    arg.specie if arg.specie else Param.default_specie,
    arg.gene_col,
    arg.log2fc_col,
    arg.pvalue_col,
    arg.output_file,
    arg.focus_genes_file,
    arg.aggregation_method,
    arg.precomputed,
)
