def read_fasta(filepath):
    """
    Given a path to a FASTA file, return the header and sequence.
    """
    with open(filepath, 'r') as fasta_file:
        header = fasta_file.readline().strip()
        sequence = "".join([line.strip()
                            for line in fasta_file.readlines()])
    return (header, sequence)


print(read_fasta('mecA_fasta.txt'))
