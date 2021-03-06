# gRNA-helper

- **gRNA-helper.py**: A primary script for this repository. Contains functionality for finding and evaluating potential sequences for knocking out a target gene in a given genome.
- **mrsa_fasta.txt**: Sample file containing the genome for a strain of MRSA. Contains multiple gene entries in FASTA format.
- **mrsa_fasta_mini.txt**: Sample file containing a truncated version of the above to reduce runtime during development.


## Usage

1. Make sure **gRNA-helper.py** and an annotated FASTA file representing a genome are in the same directory.
2. Run **gRNA-helper.py**: `python3 gRNA-helper.py`
3. Enter a path to the desired FASTA file.
4. Search for and select a gene to knock out.
5. View results in output file (Default: _output.csv_)
