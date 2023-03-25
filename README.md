# AMINE (Active Module Identification through Network Embedding)

Implementation of the active module identification method described in the manuscript:

>Pasquier, C., Guerlais, V., Pallez, D., Rapetti-Mauss, R., & Soriani, O. (2021). Identification of active modules in interaction networks using node2vec network embedding. BioRxiv.*

## Installation

Clone the repository and change to the project directory
```bash
git clone https://github.com/claudepasquier/amine.git
cd amine
```

To set up the required dependencies for the project, you should create a virtual environment using your preferred package manager. For example, you can use mamba to create the environment including all the dependencies with the following commands:
```bash
mamba env create -n amine-env -f environment.yml
mamba activate amine-env
```
Mamba is a better choice than conda for this task, as it can perform the installation faster than conda version 23.1.0, which can take more than 10 minutes.

Alternatively, you can create the environment and install manually the dependencies using conda with the following commands:
```bash
conda create -n amine-env python=3.6
conda activate amine-env
conda install -c conda-forge -c anaconda -c numba networkx scipy gensim numba pandas xlrd scikit-learn powerlaw progressbar2 openpyxl python-levenshtein
```

At this point, you can download network files with the following command:
```bash
python -m amine.download_data
```
that will download data in the ./data/real/ppi directory

## networks
The downloaded protein-protein interaction networks originating from 3 different databases and concerning 3 different species are described below.

### from STRING database (https://string-db.org/)
* human: network of 19,382 nodes and 5,968,679 edges
* mouse: 21,317 nodes and 7,248,179 edges
* drosophila: network of 13,047 nodes and 2,171,899 edges

### from IntAct database (https://www.ebi.ac.uk/intact)
* human: network of 17,721 nodes and 314,807 edges
* mouse: network of 6,998 nodes and 18,944 edges
* drosophila: network of 6,610 nodes and 29,305 edges

### from BioGRID database (https://thebiogrid.org/)
* human: network of 19,892 nodes and 780,328 edges
* mouse: network of 10,949 nodes and 60,254 edges
* drosophila: network of 9,514 nodes and 64,927 edges

### merged networks
In addition to the interactions networks stored in specific databases, we propose, for each specie, a network that contains a fusion of all interactions stored in Intact, BioGRID or STRING (for this database, only the interactions associated with a global score >= 0.7 are retained).
* for human, this represents a network with 21,936 nodes and 1,056,188 edges
* for mouse, this represents a network with 18,156 nodes and 283,649 edges
* for drosophila, this represents a network with 14,412 nodes and 211,618 edges

## Usage

### On real data
The program can be executed with the command:
```bash
python -m amine.process_real_network
```
The -h option displays explanations of the available options. A typical execution can be performed with:
```bash
python -m amine.process_real_network --expvalue ./data/real/expression/chiou_2017/Hmga2_positive_vs_negative.csv -g 0 -l 2 -p 6 -s mouse -n string -o ./data/results/Hmga2_positive_vs_negative_string_network.xlsx -v
```
The command above runs amine on the result of a differential expression analysis stored in the file specified with the parameter **--expvalue**. The parameters **-g**, **-l** and **-p** are used to specify the column with the gene names, the log2 fold changes and the p-values respectively. Numbering starts at 0, which means that zero identifies the first column. The parameter **-s** is used to indicate the specie and the parameter **-n** is used to indicate the origin of the interaction network (in the example, it is the STRING database). The parameter **-o** allows to specify the path to a file to write the results.

The file "execute_reals.sh" contains examples of commands to process data from the paper of Chiou 2017 using different networks.

### On artificial data
The method can be executed in batch mode on set of artificially generated datasets with the command:
```bash
python -m amine.process_artificial_network
```
The -h option displays explanations of the available options. Below are two examples of commands to run the program on a batch of artificial data:

* running AMINE on the dataset of 1,000 networks with 3 modules of 10 nodes used in the Robinson et al. paper and saving the results in the file "./data/results/guyon.txt":
```bash
python -m amine.process_artificial_network -g guyondata -r 1000 -m 10 -n 3 -o ./data/results/guyon.txt -v
```
* executing AMINE on 100 dense networks of 1,000 nodes with 1 modules of 20 nodes and saving the results in the file "./data/results/dense_1Knodes_moduleof20.txt":
```bash
python -m amine.process_artificial_network -r 100 -m 20 -s 1000 -o ./data/results/dense_1Knodes_moduleof20.txt -v
```
