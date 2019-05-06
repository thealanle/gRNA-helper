import re


class Genome():
    """
    A collection of Gene objects. Contains functionality for iterating over
    itself and returning matches.
    """

    def __init__(self, fasta_file=None):
        self.genes = {}

        # If the Genome instance is called with a path to a FASTA file, load it.
        if fasta_file:
            self.load_fasta(fasta_file)
            print(f"{fasta_file} loaded.")

# define the function
    def load_fasta(self, fasta_file):
        """
        Given a file in FASTA format, turn each entry into a header string and
        a sequence string, then pass those headers and sequences into new Gene
        instances. Then add each of those gene instances to self.genes.
        """
        # initialize dict and lists to create elements or dict
        # header_seq_dict = {}
        header = []
        sequences = []
        my_fasta = open(fasta_file, 'r')
        fasta_content = my_fasta.read()
    # re for header to create list of headers
        header = re.findall(r'>lcl.+[.+].+', fasta_content)
    # re for sequence and fix sequence to only be string of nucleotides
        sequence = re.findall(r'(([ATCG]+\n)+)', fasta_content)
        for i in sequence:
            DNA_seq = []
            DNA_seq = i[0]
            DNA_fixed = DNA_seq.replace('\n', '')
            sequences.append(DNA_fixed)
    # creation of dictionary
        # header_seq_dict = dict(zip(header, sequences))
        # return(header_seq_dict)
        for h, s in zip(header, sequences):
            self.add_gene(h, s)

    def add_gene(self, header, sequence):
        gene = Gene(header, sequence)
        id = gene.info['gene_id']
        self.genes[id] = gene

    def find_hits(self, target_sequence):
        """
        <<WORK IN PROGRESS>>
        Given a target sequence (representing a k-nucleotide-long sequence and a
        PAM), return a list of hits, the genes in which they occur, and the
        locations of the hits within the genome.

        SAMPLE: "Sequence CGATCGTAATGCTCA hit found in [GENE] at [INDEX]"
        """

        hits = []
        ko_target = TargetGene(MRSA_HEADER, MRSA_SEQUENCE)
        spacers = list(ko_target.find_protospacers())
        for spacer in spacers:
            spacer_seq = spacer[0]
            expression = r'(?=(' + spacer_seq + r'))'
            hits = re.finditer(expression, target_genome)
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
        results = [gene for gene in self.genes.values() if query.lower() in
                   gene.info['protein'].lower()]
        menu = '\n'.join([f"{i}) {gene.info['protein']}" for gene, i in zip(results, range(1, len(results) + 1))])
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
        if 'location' in d.keys():
            # Ignore "complement" and "join" parameters
            location = re.sub(r'[(complement)(join)<>\(\)]', '', d['location'])
            if 'complement' in location:
                # location = location.replace('complement', '')
                self.on_complement_strand = True
            loc_split = location.split('..')
            start, end = int(loc_split[0]), int(loc_split[-1])
            # start, end = [int(i) for i in location.split('..')]
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
        print(f"Finding protospacers in target gene {self.info['gene_id']}.")
        p_length = k + len(pam_sequence)  # Length of protospacer + PAM
        pam_sequence = self.seq_to_regex(pam_sequence)
        expression = r'(?=([ACGT]{' + str(k) + r'}' + str(pam_sequence) + r'))'
        gene_start = self.info['location']['start']
        gene_end = self.info['location']['end']
        print("start and end: ", gene_start, gene_end)
        count = 0

        # Search the forward strand
        for result in re.finditer(expression, self.sequence):
            match = result.group(1)
            match_start = gene_start + result.start()
            match_end = match_start + p_length - 1
            span = (match_start, match_end)
            print(match, span)
            count += 1
            yield match, span

        # Search the complementary strand
        for result in re.finditer(expression, self.get_complement()):
            match = result.group(1)
            match_start = gene_end - result.start()
            match_end = match_start - p_length + 1
            span = (match_start, match_end)
            print(match, span, "complement strand")
            count += 1
            yield match, span

        print(f"{count} protospacers identified.")

    def seq_to_regex(self, sequence):
        """
        Convert a gene sequence to a regular expression snippet and return the
        result.
        """

        # This section lacks complete functionality, so it's overly verbose
        # for what it currently accomplishes.
        result = sequence

        substitution_dict = {
            'N': '[ACGT]',
            'R': '[GA]',
            'Y': '[TC]',
            'M': '[AC]',
            'S': '[GC]',
            'W': '[AT]',
            'K': '[GT]',
        }

        for sub in substitution_dict:
            result = result.replace(sub, substitution_dict[sub])

        return result


def get_complement(sequence):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result = ''.join([pairs[nuc] for nuc in sequence[::-1]])
    return result


MRSA_HEADER = """>lcl|CP015447.2_cds_ARI72360.1_79 [locus_tag=A6V38_00405] [protein=PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA] [protein_id=ARI72360.1] [location=complement(76807..78816)] [gbkey=CDS]"""

MRSA_SEQUENCE = """TTGATGAAAAAGATAAAAATTGTTCCACTTATTTTAATAGTTGTAGTTGTCGGGTTTGGTATATATTTTTATGCTTCAAAAGATAAAGAAATTAATAATACTATTGATGCAATTGAAGATAAAAATTTCAAACAAGTTTATAAAGATAGCAGTTATATTTCTAAAAGCGATAATGGTGAAGTAGAAATGACTGAACGTCCGATAAAAATATATAATAGTTTAGGCGTTAAAGATATAAACATTCAGGATCGTAAAATAAAAAAAGTATCTAAAAATAAAAAACGAGTAGATGCTCAATATAAAATTAAAACAAACTACGGTAACATTGATCGCAACGTTCAATTTAATTTTGTTAAAGAAGATGGTATGTGGAAGTTAGATTGGGATCATAGCGTCATTATTCCAGGAATGCAGAAAGACCAAAGCATACATATTGAAAATTTAAAATCAGAACGTGGTAAAATTTTAGACCGAAACAATGTGGAATTGGCCAATACAGGAACAGCATATGAGATAGGCATCGTTCCAAAGAATGTATCTAAAAAAGATTATAAAGCAATCGCTAAAGAACTAAGTATTTCTGAAGACTATATCAAACAACAAATGGATCAAAAGTGGGTACAAGATGATACCTTCGTTCCACTTAAAACCGTTAAAAAAATGGATGAATATTTAAGTGATTTCGCAAAAAAATTTCATCTTACAACTAATGAAACAGAAAGTCGTAACTATCCTCTAGAAAAAGCGACTTCACATCTATTAGGTTATGTTGGTCCCATTAACTCTGAAGAATTAAAACAAAAAGAATATAAAGGCTATAAAGATGATGCAGTTATTGGTAAAAAGGGACTCGAAAAACTTTACGATAAAAAGCTCCAACATGAAGATGGCTATCGTGTCACAATCGTTGACGATAATAGCAATACAATCGCACATACATTAATAGAGAAAAAGAAAAAAGATGGCAAAGATATTCAACTAACTATTGATGCTAAAGTTCAAAAGAGTATTTATAACAACATGAAAAATGATTATGGCTCAGGTACTGCTATCCACCCTCAAACAGGTGAATTATTAGCACTTGTAAGCACACCTTCATATGACGTCTATCCATTTATGTATGGCATGAGTAACGAAGAATATAATAAATTAACCGAAGATAAAAAAGAACCTCTGCTCAACAAGTTCCAGATTACAACTTCACCAGGTTCAACTCAAAAAATATTAACAGCAATGATTGGGTTAAATAACAAAACATTAGACGATAAAACAAGTTATAAAATCGATGGTAAAGGTTGGCAAAAAGATAAATCTTGGGGTGGTTACAACGTTACAAGATATGAAGTGGTAAATGGTAATATCGACTTAAAACAAGCAATAGAATCATCAGATAACATTTTCTTTGCTAGAGTAGCACTCGAATTAGGCAGTAAGAAATTTGAAAAAGGCATGAAAAAACTAGGTGTTGGTGAAGATATACCAAGTGATTATCCATTTTATAATGCTCAAATTTCAAACAAAAATTTAGATAATGAAATATTATTAGCTGATTCAGGTTACGGACAAGGTGAAATACTGATTAACCCAGTACAGATCCTTTCAATCTATAGCGCATTAGAAAATAATGGCAATATTAACGCACCTCACTTATTAAAAGACACGAAAAACAAAGTTTGGAAGAAAAATATTATTTCCAAAGAAAATATCAATCTATTAACTGATGGTATGCAACAAGTCGTAAATAAAACACATAAAGAAGATATTTATAGATCTTATGCAAACTTAATTGGCAAATCCGGTACTGCAGAACTCAAAATGAAACAAGGAGAAACTGGCAGACAAATTGGGTGGTTTATATCATATGATAAAGATAATCCAAACATGATGATGGCTATTAATGTTAAAGATGTACAAGATAAAGGAATGGCTAGCTACAATGCCAAAATCTCAGGTAAAGTGTATGATGAGCTATATGAGAACGGTAATAAAAAATACGATATAGATGAATAA"""

# TEMPORARY GENOME FILE READING, concatenates all sequences into a single string
# with open('mrsa_fasta.txt') as f_in:
#     target_genome = ''.join([l.strip() for l in f_in if not l.startswith('>')])

target_genome = Genome('mrsa_fasta.txt')
target = target_genome.choose_gene('peptidase')
ko_target = TargetGene(MRSA_HEADER, MRSA_SEQUENCE)
spacers = list(ko_target.find_protospacers())

matches = {}

# To-do: convert hit.span to take into account the position of the gene
for gene in target_genome.genes.values():
    for spacer in spacers:
        spacer_seq = spacer[0]
        expression = r'(?=(' + spacer_seq + r'))'
        hits = re.finditer(expression, gene.sequence)
        # seq_length = len(gene.sequence)
        gene_start = gene.info['location']['start']
        gene_end = gene.info['location']['end']

        # Iterate over the forward strand
        for hit in hits:
            print(f"{spacer_seq} found at ({hit.start(1) + gene_start}, {hit.end(1) + gene_start - 1})")
            matches[spacer_seq] = (hit, hit.span(1))
            # if gene.info['gene_id'] != ko_target.info['gene_id']:
            #     print("off-target effect")

        # Iterate over the complement strand
        hits = re.finditer(expression, get_complement(gene.sequence))
        for hit in hits:
            print(f"{spacer_seq} found at ({gene_end - hit.start(1)}, {gene_end - hit.end(1) + 1})")
            matches[spacer_seq] = (hit, hit.span(1))
            # if gene.info['gene_id'] != ko_target.info['gene_id']:
            #     print("off-target effect")

# print("number of spacers: ", len(spacers))
# print("number of matches: ", len(matches))
# for x in matches:
#     print(matches[x][0], matches[x][1])
