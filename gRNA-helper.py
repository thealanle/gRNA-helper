import re


class Genome():
    """
    A collection of Gene objects. Contains functionality for iterating over
    itself and returning matches.
    """

    def __init__(self):
        # DO STUFF
        pass

    def find_hits(self, target_sequence):
        result = []
        # DO STUFF
        return result


class Gene():
    """
    Given a header and sequence from a FASTA file, creates a representation of
    a gene.

    """

    def __init__(self, header, sequence):
        self.info = {}
        self.header = header
        self.sequence = sequence
        # TO-DO: Parse header values into a dict?

    def parse_header(self, header):
        # Load header data into self.info
        pass

    def get_complement(self):
        return self.sequence[::-1]


class TargetGene(Gene):
    """
    A Gene object that contains functionality for finding eligible knockout
    targets. By default, an eligible sequence is a 20-mer followed by a
    terminal "NGG".
    """

    def find_spacers(self, k=20, pam_sequence='NGG'):
        pam_sequence = self.seq_to_regex(pam_sequence)
        expression = r'(?=([ACGT]{' + str(k) + r'}' + str(pam_sequence) + r'))'
        for result in re.finditer(expression, self.sequence):
            matched_sequence = result.group(1)
            matched_index = int(result.start())
            yield matched_sequence, matched_index

    def seq_to_regex(self, sequence):
        """
        Convert a gene sequence to a regular expression snippet and return the
        result.
        """
        result = sequence
        result = result.replace('N', '[ACGT]')
        return result


def read_fasta(filepath):
    """
    Given a path to a FASTA file, return the headers and sequences.

    TODO: Implement ability to intake multiple FASTA entries from a single
    file, e.g. mrsa_fasta.txt.
    """
    header = ''
    sequence = ''

    # DO STUFF

    return header, sequence

# TODO: Move this into a separate file...


MECA_HEADER = """>lcl|KF058908.1_cds_AGU99981.1_1 [gene=mecA] [protein=penicillin-binding protein] [protein_id=AGU99981.1] [location=<1..>1983] [gbkey=CDS]"""

MECA_SEQUENCE = """ATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAATTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGGAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAAT"""

NAT_SEQUENCE = """ATGAGCATAATTACAAGATTGTTTAATAACAGTGATTTTGAAAAATTAAATCAACTATGTAAATTATATGATGATCTAGGTTATCCAACAAATGAGAATGATTTAAAAAAGAGACTAAAGAAAATAACGAATCATGATGATTACTTCCTACTGCTTTTGATAAAAGAAAATAAAATAATTGGTTTAAGTGGTATGTGTAAAATGATGTTTTACGAAAAAAATGCAGAGTATATGAGAATCCTTGCGTTTGTTATACATTCTGAATTTAGGAAAAAAGGTTATGGAAAGAGATTATTAGCTGATTCTGAAGAATTTTCTAAACGGTTGAATTGTAAAGCAATAACACTAAATAGTGGTAATAGAAATGAAAGACTATCTGCACATAAACTATATAGTGATAATGGGTATGTTAGCAATACATCTGGGTTTACTAAACAACTATAA"""

# TEMPORARY GENOME FILE READING, concatenates all sequences into a single string
with open('mrsa_fasta.txt') as f_in:
    genome_data = ''.join([l.strip() for l in f_in if not l.startswith('>')])

matches = {}

meca_target = TargetGene(MECA_HEADER, MECA_SEQUENCE)
spacers = list(meca_target.find_spacers())
for spacer in spacers:
    spacer_seq = spacer[0]
    expression = r'(?=(' + spacer_seq + r'))'
    hits = re.finditer(expression, genome_data)
    for hit in hits:
        print(f"{spacer_seq} found at {hit.span(1)}")
        matches[spacer_seq] = (hit, hit.span(1))

print("spacers: ", len(spacers))
print("matches: ", len(matches))
# for x in matches:
#     print(matches[x][0], matches[x][1])
