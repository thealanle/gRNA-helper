import re


class Genome():
    """
    A collection of Gene objects. Contains functionality for loading Genes from
    a FASTA file.
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
        id = gene.id
        self.genes[id] = gene

    def choose_gene(self, query):
        """
        Given a search query, return a menu of all results and return the chosen
        result.
        """
        # Create a list of Gene objects whose "protein" fields contain the query
        results = [gene for gene in self.genes.values() if query.lower() in
                   gene.info['protein'].lower()]
        if len(results) == 0:
            print("Search found 0 results.")
            return None
        menu = '\n'.join([f"{i}) {gene.info['protein']}" for gene, i in zip(results, range(1, len(results) + 1))])
        print(menu)
        while True:
            try:
                choice = int(input("Choose a gene.\n>")) - 1
                return results[choice]
            except TypeError as ex:
                print("Invalid choice.")
            except IndexError as ex:
                print("Selection out of range.")
            except ValueError as ex:
                print("Invalid value. Please enter a number.")


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
        self.hits = []

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
        self.id = results[0].group(1)

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

    def find_hits(self, spacers):
        """
        Given a target sequence (representing a k-nucleotide-long sequence and a
        PAM), return a list of hits and the positions at which they occur.
        """

        for spacer in spacers:
            spacer_seq = spacer[0]
            gene_start = self.info['location']['start']
            gene_end = self.info['location']['end']

            expression = r'(?=(' + spacer_seq + r'))'
            hits = re.finditer(expression, self.sequence)
            # Iterate over the forward strand
            for hit in hits:
                spacer_start = hit.start(1) + gene_start
                spacer_end = hit.end(1) + gene_start - 1
                self.hits.append((spacer_seq, spacer_start, spacer_end))

            # Iterate over the complement strand
            hits = re.finditer(expression, self.get_complement())
            for hit in hits:
                spacer_start = gene_end - hit.start(1)
                spacer_end = gene_end - hit.end(1) + 1
                self.hits.append((spacer_seq, spacer_start, spacer_end))


class TargetGene(Gene):
    """
    A Gene object that contains functionality for finding eligible knockout
    targets. By default, an eligible sequence is a 20-mer followed by a
    terminal "NGG".

    To create a TargetGene from an existing Gene, use:
        [NEW_TARGET_GENE] = TargetGene([GENE].header, [GENE].sequence)
    """

    def find_protospacers(self, k=20, pam_sequence='NGG'):
        print(f"Finding protospacers in target gene {self.id} \"{self.info['protein']}\".")
        p_length = k + len(pam_sequence)  # Length of protospacer + PAM
        pam_sequence = self.seq_to_regex(pam_sequence)
        expression = r'(?=([ACGT]{' + str(k) + r'}' + str(pam_sequence) + r'))'
        gene_start = self.info['location']['start']
        gene_end = self.info['location']['end']
        count = 0

        # Search the forward strand
        for result in re.finditer(expression, self.sequence):
            match = result.group(1)
            match_start = gene_start + result.start()
            match_end = match_start + p_length - 1
            span = (match_start, match_end)
            # print(match, span)
            count += 1
            yield match, span

        # Search the complementary strand
        for result in re.finditer(expression, self.get_complement()):
            match = result.group(1)
            match_start = gene_end - result.start()
            match_end = match_start - p_length + 1
            span = (match_start, match_end)
            # print(match, span, "complement strand")
            count += 1
            yield match, span

        print(f"{count} protospacers identified.")

    def seq_to_regex(self, sequence):
        """
        Convert a gene sequence to a regular expression snippet and return the
        result.
        """

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
            sequence = sequence.replace(sub, substitution_dict[sub])

        return sequence


def dna_to_rna(sequence):
    """
    Given a protospacer, return the matching spacer.
    """
    return sequence.replace('T', 'U')


def choose_file():
    """
    Prompt the user for a filename and return the filename if the file exists.
    """
    while True:
        try:
            response = input("Enter the name of a FASTA file:\n>")
            f_in = open(response, 'r')
            f_in.close()
            return response
        except IOError as ex:
            print(ex.strerror)


def main():
    OUTPUT_FILE = 'output.csv'

    genome = Genome(choose_file())

    target = None
    while target is None:
        target = genome.choose_gene(input('Enter gene name:\n>'))
    ko_target = TargetGene(target.header, target.sequence)

    spacers = set(ko_target.find_protospacers())

    off_target = []
    for gene in genome.genes.values():
        gene.find_hits(spacers)
        if len(gene.hits) > 0:
            print(f"{len(gene.hits)} hits found in {gene.id} \"{gene.info['protein']}\"")
            if gene.id != ko_target.id:
                off_target.extend([spacer[0] for spacer in gene.hits])

    # Truncate the PAM-matching sequence from the protospacer and convert the
    # result to an RNA sequence. In addition, append the number of off-targets.
    guides = [
        [dna_to_rna(spacer[0][:-3]),
         off_target.count(spacer[0])] for spacer in spacers]
    guides = sorted(guides, key=lambda x: x[1])

    # Header for the output report.
    log = [['Knockout Target', ko_target.id, ko_target.info['protein']], [],
           ['gRNA Candidate', 'Off-Target Hits']]

    # Print and log the spacers and their number of off-target hits.
    for guide in guides:
        log.append([guide[0], str(guide[1])])
        # print(f"{guide[0]}: {guide[1]} off-target hit(s).")

    # Open an output file and write the log contents to it.
    with open(OUTPUT_FILE, 'w') as f_out:
        for line in log:
            f_out.write(','.join(line) + '\n')

        print(f"Results have been output to {OUTPUT_FILE}.")


if __name__ == '__main__':
    main()
