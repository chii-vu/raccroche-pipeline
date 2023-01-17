class Genome:
    """
    A Class Used to Represent Genomes

    Methods
    -------
    get_genome_in_string()
        restructures and prints out chromosomes and their associated genes

    get_gene_list(genome_id)
        accesses gene_list and retrieves information corresponding to the genome_id key
    """

    def __init__(self, gene_list, all_leaves, parsed_config):
        """
        :param gene_list: dictionary
            stores the names of all genes, respective gene families along with other related information

        :param all_leaves: list
            a list containing all the information related to the inputted genomes

        :param parsed_config: list
            a list containing the file path for the output

        """

        self.gene_list = gene_list
        self.all_leaves = all_leaves

        self.output_path = parsed_config['output_path']
        self.project_name = parsed_config['project_name']

        self.genomes_output = self.output_path + self.project_name + parsed_config['output_file']['genomes']

    def get_genome_in_string(self):
        """ Groups genes by their respective chromosomes that they belong to and generates output

        For each genome, find its associated genes. Then, create a dictionary with keys being the genome's
        chromosomes and the values being their corresponding genes (in the format of gene_family_labels).
        Finally, write this dictionary to an output file.

        :return: genome_in_string: dict
            contains the chromosomes for each genome

            key: id of each genome (str)
            value: all the chromosomes (list)
        """

        genome_in_string = {}

        open(self.genomes_output, 'w').close()

        for leaf in self.all_leaves:
            with open(self.genomes_output, 'a') as f:
                f.write("#" + str(leaf[1]) + "\n")

            genes_in_chromosome = {}
            genome_id = leaf[0]
            gene_id_list = self.get_gene_list(genome_id)

            for gene in gene_id_list:
                if gene[0] in genes_in_chromosome:
                    genes_in_chromosome[gene[0]].append(gene[-1])
                else:
                    genes_in_chromosome[gene[0]] = [gene[-1]]

            chromosomes = []

            for genes in genes_in_chromosome.values():
                chromosomes.append(genes)

            genome_in_string[genome_id] = chromosomes

            # print(len(genes_in_chromosome))

            # for k, v, in genes_in_chromosome.items():
                # print(k, v)

            with open(self.genomes_output, 'a') as f:
                for chromosome, genes in genes_in_chromosome.items():
                    f.write(chromosome + " ")

                    for gene in genes:
                        f.write("%s " % genes)

                    f.write("\n")

            # print(gene_id_list)

        return genome_in_string

    def get_gene_list(self, genome_id):
        """ Retrieves information corresponding to the given genome_id

        For each gene in gene_list, extract various pieces of information if the gene exists in the genome
        (identified by genome_id), then add it into a list. Then sort the list by chromosome, then by position.

        :param genome_id: str
        :return: info_list: list
        """

        info_list = []

        for gene in self.gene_list:
            all_gene_info = self.gene_list[gene]

            if genome_id == all_gene_info[0]:
                chromosome = all_gene_info[1][0]
                position = all_gene_info[1][1]
                name = all_gene_info[1][3]
                orientation = all_gene_info[1][4]
                gene_family_label = int(all_gene_info[-1]) * int(orientation)

                single_gene_info = [chromosome, position, name, orientation, gene_family_label]

                info_list.append(single_gene_info)

        # First sorts by genes' chromosome, then by genes' start position
        info_list.sort(key=lambda x: (x[0], int(x[1])))

        return info_list
