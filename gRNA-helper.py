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

        SAMPLE: "Sequence CGATCGTAATGCTCA hit found in [GENE] at [INDEX]"
        """
        hits = []
        # DO STUFF
        return hits


class Gene():
    """
    Given a header and sequence from a FASTA file, create a representation of
    a gene. Parse the header parameters and store them in the dictionary
    self.info.

    """

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.info = self.parse_header(self.header)
        self.complement = None  # This is None until get_complement is called

    def parse_header(self, header):
        """
        Given a header in FASTA format, parse the parameters and load them into
        the dict self.info.
        """

        d = {}
        expression = r'>lcl\|([\w\.]*)|\[([^=]+)=([^\]]+)]'
        results = list(re.finditer(expression, header))

        # There is only one instance of Group 1, which contains database and
        # gene ID information.
        d['gene_id'] = results[0].group(1)

        # Iterate over the remaining match groups, which contain keys and values
        # that can be added to the dictionary self.info.
        for result in results[1:]:
            d[result.group(2)] = result.group(3)

        # 'location' data must be parsed from a string into integers.
        # self.info['location'] should return a pair of ints as a tuple.
        if 'location' in d.keys():
            location = re.sub(r'[<>]', '', d['location'])
            d['location'] = tuple([int(i) for i in location.split('..')])

        return d

    def get_complement(self):
        """
        Build and return a complement strand.
        """

        if self.complement is None:
            reverse = self.sequence[::-1]
            pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            self.complement = ''.join([pairs[nuc] for nuc in reverse])
        return self.complement


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

MRSA_HEADER = """>lcl|CP015447.2_cds_ARI72360.1_79 [locus_tag=A6V38_00405] [protein=PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA] [protein_id=ARI72360.1] [location=complement(76807..78816)] [gbkey=CDS]"""

MRSA_SEQUENCE = """TTGATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAAGTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGAAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAATAAAAAATACGATATAGATGAATAA"""

# TEMPORARY GENOME FILE READING, concatenates all sequences into a single string
with open('mrsa_fasta.txt') as f_in:
    genome_data = ''.join([l.strip() for l in f_in if not l.startswith('>')])

"""
matches = {}

meca_target = TargetGene(MRSA_HEADER, MRSA_SEQUENCE)
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
"""

meca = Gene(MECA_HEADER, MECA_SEQUENCE)
