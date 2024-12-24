# Virus recombination simulator
WIP 24/12/24

Tool to create biologically-feasible recombinant virus genomes by randomly transposing a donor gene to an acceptor genome (copy choice by homologous template, at same position). Originally developed for testing [Castanet's](https://github.com/MultipathogenGenomics/castanet) consensus generation algorithm (and creates output summary statistics for this purpose), but may be coopted to generate synthetic recombinant genomes.

## Prerequisites
1. [Castanet](https://github.com/MultipathogenGenomics/castanet) V >= 8.0
1. [DWGSIM](https://github.com/nh13/DWGSIM)
1. Python3.x
1. Pip

## Installation
1. Clone
1. ```pip install -r requirements```
1. Start Castanet server at internal address *x*
1. Modify main namespace (end of file) in app.pipeline to specify n repeats ("n"), Castanet server endpoint *x* ("endpoint") and dwgsim path ("dwgsim_path")
1. ```$ python3 -m app.pipeline```

## Roadmap
1. Argparser
1. Instructions for generating new datasets with dev/pull_genbank_annotations.py
1. Tests