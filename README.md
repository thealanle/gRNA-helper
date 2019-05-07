# gRNA-helper

- gRNA-helper.py: A primary script for this repository. Contains functionality for finding and evaluating potential sequences for knocking out a target gene in a given genome.
- mrsa_fasta.txt: The genome for a strain of MRSA. Contains multiple gene entries in FASTA format.
- mrsa_fasta_mini.txt: A truncated version of the above to reduce runtime during development.


## Usage

1. Make sure gRNA-helper.py and a valid FASTA file representing a genome are in the same directory.
2. Call gRNA-helper.py: `python3 gRNA-helper.py`
3. Enter a path to the desired FASTA file.
4. Search for and select a gene to knock out.
5. View results in output file (Default: output.csv)
