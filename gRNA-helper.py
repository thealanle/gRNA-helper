def read_fasta(filename):
    with open(filename, 'r') as fasta_file:
        return "".join([line.strip() for line in fasta_file])[1:]
