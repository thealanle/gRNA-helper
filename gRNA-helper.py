import re


class Genome():
    """
    A collection of Gene objects. Contains functionality for iterating over
    itself and returning matches.
    """

    def __init__(self, fasta_file=None):
        self.gene_dict = {}

        # If the Genome instance is called with a path to a FASTA file, load it.
        if fasta_file:
            self.load_fasta(fasta_file)

    def load_fasta(self, fasta_file):
        """
        Given a file in FASTA format, turn each entry into a header string and
        a sequence string, then pass those headers and sequences into new Gene
        instances. Then add each of those gene instances to self.gene_dict.
        """

        pass

    def find_hits(self, target_sequence):
        """
        Given a target sequence (representing a k-nucleotide-long sequence and a
        PAM), return a list of hits, the genes in which they occur, and the
        locations of the hits within the genome.

        SAMPLE OUTPUT: "Sequence CGATCGTAATGCTCA hits found in [GENE] at
        [START_INDEX]"
        """
        result = []
        # DO STUFF
        return result


class Gene():
    """
    Given a header and sequence from a FASTA file, creates a representation of
    a gene. The parsed header parameters are stored in the dictionary self.info.

    """

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.info = {}
        # TO-DO: Parse header values into a dict?

    def parse_header(self, header):
        """
        Given a header in FASTA format, parse the parameters and load them into
        the dict self.info.
        """
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


# TODO: Move this into a separate file...


MECA_HEADER = """>lcl|KF058908.1_cds_AGU99981.1_1 [gene=mecA] [protein=penicillin-binding protein] [protein_id=AGU99981.1] [location=<1..>1983] [gbkey=CDS]"""

MECA_SEQUENCE = """ATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAATTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGGAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAAT"""

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
