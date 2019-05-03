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
        <<WORK IN PROGRESS>>
        Given a target sequence (representing a k-nucleotide-long sequence and a
        PAM), return a list of hits, the genes in which they occur, and the
        locations of the hits within the genome.

        SAMPLE: "Sequence CGATCGTAATGCTCA hit found in [GENE] at [INDEX]"
        """

        hits = []
        meca_target = TargetGene(MRSA_HEADER, MRSA_SEQUENCE)
        spacers = list(meca_target.find_protospacers())
        for spacer in spacers:
            spacer_seq = spacer[0]
            expression = r'(?=(' + spacer_seq + r'))'
            hits = re.finditer(expression, genome_data)
            for hit in hits:
                print(f"{spacer_seq} found at {hit.span(1)}")
                matches[spacer_seq] = (hit, hit.span(1))
        return hits

    def choose_gene(self, query):
        """
        <<WORK IN PROGRESS>>
        Given a search query, return a menu of all results and return the chosen
        result.
        """
        # Create a list of Gene objects whose "protein" fields contain the query
        results = [gene for gene in self.gene_dict.values() if
                   gene.info['protein'].lower().contains(query.lower())]
        menu = '\n'.join([f"{i}) {gene.info['protein']}" for gene, i in (results, range(1, len(results) + 1))])
        print(menu)

        choice = int(input("Choose a gene.\n>")) - 1
        return results[choice]


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
            location = re.sub(r'[(complement)<>\(\)]', '', d['location'])
            if 'complement' in location:
                # location = location.replace('complement', '')
                self.on_complement_strand = True
            start, end = [int(i) for i in location.split('..')]
            d['location'] = {'start': start, 'end': end}

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

    To create a TargetGene from an existing Gene, use:
        [NEW_TARGET_GENE] = TargetGene([GENE].header, [GENE].sequence)
    """

    def find_protospacers(self, k=20, pam_sequence='NGG'):
        print("Finding protospacers...")
        p_length = k + len(pam_sequence)  # Length of protospacer + PAM
        pam_sequence = self.seq_to_regex(pam_sequence)
        expression = r'(?=([ACGT]{' + str(k) + r'}' + str(pam_sequence) + r'))'
        gene_start = self.info['location']['start']
        gene_end = self.info['location']['end']
        print("start and end: ", gene_start, gene_end, p_length)

        # Search the forward strand
        for result in re.finditer(expression, self.sequence):
            match = result.group(1)
            match_start = gene_start + result.start()
            match_end = match_start + p_length - 1
            span = (match_start, match_end)
            print(match, span)
            yield match, span

        # Search the complementary strand
        for result in re.finditer(expression, self.get_complement()):
            match = result.group(1)
            match_start = gene_end - result.start()
            match_end = match_start - p_length + 1
            span = (match_start, match_end)
            print(match, span, "complement strand")
            yield match, span

        print("All protospacers identified.")

    def seq_to_regex(self, sequence):
        """
        Convert a gene sequence to a regular expression snippet and return the
        result.
        """

        # This section lacks complete functionality, so it's overly verbose
        # for what it currently accomplishes.
        result = sequence
        result = result.replace('N', '[ACGT]')
        return result


def get_complement(sequence):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result = ''.join([pairs[nuc] for nuc in sequence[::-1]])
    return result


MRSA_HEADER = """>lcl|CP015447.2_cds_ARI72360.1_79 [locus_tag=A6V38_00405] [protein=PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA] [protein_id=ARI72360.1] [location=complement(76807..78816)] [gbkey=CDS]"""

MRSA_SEQUENCE = """TTGATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAAGTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGAAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAATAAAAAATACGATATAGATGAATAA"""

# TEMPORARY GENOME FILE READING, concatenates all sequences into a single string
with open('mrsa_fasta.txt') as f_in:
    genome_data = ''.join([l.strip() for l in f_in if not l.startswith('>')])

matches = {}

meca_target = TargetGene(MRSA_HEADER, MRSA_SEQUENCE)
spacers = list(meca_target.find_protospacers())

for spacer in spacers:
    spacer_seq = spacer[0]
    expression = r'(?=(' + spacer_seq + r'))'
    hits = re.finditer(expression, genome_data)
    for hit in hits:
        print(f"{spacer_seq} found at {hit.span(1)}")
        matches[spacer_seq] = (hit, hit.span(1))
    hits = re.finditer(expression, get_complement(genome_data))
    for hit in hits:
        print(f"{spacer_seq} found at {hit.span(1)}")
        matches[spacer_seq] = (hit, hit.span(1))

print("spacers: ", len(spacers))
print("matches: ", len(matches))
# for x in matches:
#     print(matches[x][0], matches[x][1])
