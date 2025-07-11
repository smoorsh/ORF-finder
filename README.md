# ORF-finder
A python script that searches for open reading frames (ORFs) in a FASTA-formatted DNA file.

## Introduction
This program takes a FASTA-formatted input file and creates a .gff file that contains the FASTA sequence title, start and stop coordinates, sequence length, if the ORF was found reading in the + or - reading frame. At the bottom of the gff file are the FASTA-formatted sequences for each ORF, which could be helpful for additional analysis. Each ORF is assigned an ID, which can be found in the tabular data at the top of the .gff file. This ID is used to identify FASTA sequences at the bottom of the file, making searching for the sequence associated with each ORF relatively straightforward.

## Dependencies
This program requires Biopython. Please install Anaconda or Miniconda if you do not already have it using the following link, and follow all of their install instructions:
[
](https://www.anaconda.com/docs/getting-started/anaconda/install)

After installing Anaconda, create a conda environment with necessary dependencies using the following command line:
```
conda create -n orffinder python=3.7
conda activate orffinder
conda install biopython
```

## Download Instructions

Clone this repository and change to the repository directory.
```
git clone https://github.com/smoorsh/ORF-finder.git
```

## Usage Instructions

To use the program, please use the following command line format:
```
python3 orffinder.py /path/to/fasta/file.fa /path/to/output/file.gff desired_min_bp_length
```
Be sure to replace ```/path/to/fasta/file.fa``` with the path to your FASTA file of interest.
Change ```/path/to/output/file.gff``` to the path to your output directory and desired .gff filename.
Lastly, change ```desired_min_bp_length``` to any integer divisible by 3 that you would like your minimum bp length for each ORF to be. For example, if you change change ```desired_min_bp_length``` to 300, then all ORFs must be at least 300 bp in length to be added to your results.
