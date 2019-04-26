import re


class Genome():
    """
    A collection of Gene objects. Contains functionality for iterating over
    itself and returning matches.
    """

    def __init__(self):
        pass


class Gene():
    """
    Given a header and sequence from a FASTA file, creates a representation of
    a gene.

    """

    def __init__(self, header, sequence):
        pass


class TargetGene(Gene):
    """
    A Gene object that contains functionality for finding eligible knockout
    targets.
    """
    pass


def read_fasta(filepath):
    """
    Given a path to a FASTA file, return the headers and sequences.

    TODO: Implement ability to intake multiple FASTA entries from a single
    file, mrsa_fasta.txt.
    """
    with open(filepath, 'r') as fasta_file:

        header = fasta_file.readline().strip()
        sequence = "".join([line.strip()
                            for line in fasta_file.readlines()])
    return (header, sequence)


def get_spacers(sequence, k=20, pam_sequence='NGG'):
    """
    Given a gene sequence, find kmers that end in a pam_sequence. Return the
    matched sequences and their indices one pair at a time.

    """
    pam_sequence = seq_to_regex(pam_sequence)

    # pam_sequence = pam_sequence.replace('N', '[ACGT]')
    expression = r'(?=([ACGT]{' + str(k) + r'}' + str(pam_sequence) + r'))'
    for result in re.finditer(expression, sequence):
        matched_sequence = result.group(1)
        matched_index = int(result.start())
        yield matched_sequence, matched_index


def seq_to_regex(sequence):
    """
    Convert a gene sequence to a regular expression snippet and return the
    result.
    """
    result = sequence
    result = result.replace('N', '[ACGT]')
    return result


MECA = """ATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAATTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGGAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAAT"""

# TEMPORARY GENOME FILE READING
with open('mrsa_fasta.txt') as genome_file:
    genome_data = ''.join([line.strip()
                           for line in genome_file if not line.startswith('>')])

matches = {}
spacers = list(get_spacers(MECA, k=20, pam_sequence='NGG'))
for spacer in spacers:
    spacer_seq = spacer[0]
    hits = re.finditer(r'(?=(' + spacer_seq + r'))', genome_data)
    for hit in hits:
        print(f"{spacer_seq} found at {hit.span(1)}")
        print(hit)
        matches.setdefault(spacer_seq,)
        matches[spacer_seq] = (hit, hit.start())

print("spacers: ", len(list(spacers)))
print("matches: ", len(matches))
