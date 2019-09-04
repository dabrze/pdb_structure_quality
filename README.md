# The evolution of the quality of crystallographic models in the PDB

Reproducible PDB analysis accompanying the paper "On the evolution of the quality of crystallographic models in the PDB".

## Contents

The repository contains the input PDB metadata, reproducible experimental source code in the form of a Jupyter notebook, and detailed results of the analyses discussed in "On the evolution of the quality of crystallographic models in the PDB." by Brzezinski *et al.* The repository is divided into the following folders:

- the root folder contains the source code in the form of a Jupyter notebook, to start the analysis run the `PDB analysis.ipynb` notebook;
- `data` contains datasets with information about publications and structure quality connected with each PDB deposit;
- `results` detailed tabular results (in CSV format) corresponding to different parts of the analysis.

## Requirements

To run the experiments the following software has to be installed:

- experiments were run using Python 3.7.3;
- a complete list of modules if available in `requirements.txt` in the root folder.

## Contact

If you have trouble reproducing the experiments or have any comments/suggestions, feel free to write at dariusz.brzezinski (at) cs.put.poznan.pl