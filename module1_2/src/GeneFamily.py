class GeneFamily:
    """
    A class used to create Gene Families
    
    Attributes
    ----------
    gene_list: dictionary
        stores the names of all genes, respective gene families along with other related information

    gene_family: dictionary
        keeps track of which genes belong to which gene families

    gene_pair: list
        stores the name of both genes, their similarity and their ks

    gene_family_label: int
        a label to identify which genes are part of which gene families

    """

    gene_list = {}
    gene_family = {}
    gene_pair = []
    gene_family_label = 1

    """
    Methods
    -------
    get_leaves_and_tree_info()
        reads and stores the information from the genomes text file

    make_gene_family(genome_leaves)
        creates a gene family

    update_genes(gene_name_1, gene_info_1, leaf_index_1, gene_name_2, gene_info_2, leaf_index_2)
        updates both the gene_list and the gene_family dictionaries

    format_gene_list()
        reformats the gene_list dictionary

    extract_gene_info()
        extracts information from the pairwise SynMap output
    """

    def __init__(self, parsed_config):
        """
        :param parsed_config: list
            a list containing all the paths for the input data and output files
        """

        self.project_name = parsed_config['project_name']
        self.input_path = parsed_config['input_path']
        self.output_path = parsed_config['output_path']
        
        self.leaf_genome_file = self.input_path + self.project_name + parsed_config['input_file']['leaf_genome_info']
        self.synteny_file_ending = parsed_config['input_file']['synteny_file_name']
        self.synteny_file_path = self.input_path + self.project_name + parsed_config['input_file']['synteny_file_path']
        self.phylo_tree_path = self.input_path + self.project_name + parsed_config['input_file']['phylo_tree_path']
        self.min_cutoff_weight = parsed_config['min_cutoff_weight']
        self.max_cutoff_weight = parsed_config['max_cutoff_weight']
        
        self.gene_list_output = self.output_path + self.project_name + parsed_config['output_file']['gene_list']
        self.gene_family_output = self.output_path + self.project_name + parsed_config['output_file']['gene_family']
        self.max_genes_in_family = parsed_config['gf1']
        self.max_genes_from_genome = parsed_config['gf2']
        self.min_genomes_in_family = parsed_config['gf3']

    def get_leaves_and_tree_info(self):
        """ Gets all the information that are related to the input genomes along with the desired tree structure

        :return: all_leaves: list
            returns a list containing all the information related to the inputted genomes

        :return: median_structure: list
            returns a list containing the desired tree structure
        """

        with open(self.leaf_genome_file) as leaf_file:
            leaf_lines = leaf_file.readlines()

        all_leaves = []
        median_structure = []

        index = 0

        for line in leaf_lines:
            if line[0] != "#" and line not in ['\n', '\r\n']:
                leaf = line.split()
                all_leaves.append(leaf)
                index += 1

        with open(self.phylo_tree_path) as tree_file:
            tree_lines = tree_file.readlines()

            for line in tree_lines:
                if line[0] != "#" and line not in ['\n', '\r\n']:
                    # median_structure.append(line.rstrip().split("\t"))
                    newick_structure = line.rstrip()
                    conversion_input = newick_structure[:]

                    # print(newick_structure)
                    # print(type(conversion_input))

                    # if newick structure is in the form "genomeName_genomeID" comment out for-loop if it's just genomeID
                    for genome in all_leaves:
                        genome_string = (str(genome[1]) + "_")
                        conversion_input = conversion_input.replace(genome_string, "")

                    # print(conversion_input)
                    median_structure = self.convert_newick(conversion_input, all_leaves)

        all_leaves.sort()
        # print(median_structure)

        return all_leaves, median_structure, newick_structure

    def convert_newick(self, newick_tree, all_leaves):
        """ Converts the input tree from Newick form to a median structure for future use

        :param newick_tree: str
            input tree structure

        :param all_leaves: list
            a list containing all the information related to the inputted genomes

        :return: median_structure: list
            returns a list containing the desired tree structure
        """

        median_structure = []

        first_right = newick_tree.find(")")
        last_left = newick_tree.rfind("(", 0, first_right)

        # print("fr "+ str(first_right))
        # print("ll " +str(last_left))

        while last_left != 0:
            #print(last_left)

            neighbours = newick_tree[last_left + 1: first_right].replace(" ", "")
            # print(neighbours)

            try:
                leaf1, leaf2 = neighbours.split(",")

            except:
                print("Binary Tree Required for this Program, Correct the Inputted Tree")
                exit(1)

            leaf_pair = leaf1 + "-" + leaf2

            median_structure.append(
                [leaf1.replace("-", ","), leaf2.replace("-", ","), self.find_other_leaves(leaf1, leaf2, all_leaves)])
            newick_tree = newick_tree.replace(newick_tree[last_left: first_right + 1], leaf_pair)

            # print(median_structure)
            #  print(newick_tree)

            first_right = newick_tree.find(")")
            last_left = newick_tree.rfind("(", 0, first_right)

            # print("FR "+str(first_right))
            # print("LL "+str(last_left))

        try:
            neighbours = newick_tree[last_left + 1: first_right].replace(" ", "")
            leaf1, leaf2 = neighbours.split(",")
            # print(leaf1, leaf2)

        except:
            print("Binary Tree Required for this Program, Correct the Inputted Tree")
            exit(1)

        return median_structure

    def find_other_leaves(self, leaf1, leaf2, all_leaves):
        """ Finds all leaves that are not a part of leaf1 or leaf2

        :param leaf1: str
            a child of a node which is a string containing leaves/a leaf to exclude

        :param leaf2: str
            the other child of the same node which is also a string containing leaves/a leaf to exclude

        :param all_leaves: list
            a list containing all the information related to the inputted genomes

        :return: result: str
            a string containing the remaining leaves' genome IDs
        """

        leaves1 = leaf1.split("-")
        leaves2 = leaf2.split("-")

        # print(leaves1)
        # print(leaves2)

        leaves_to_remove = leaves1 + leaves2

        other_leaves = [leaf[0] for leaf in all_leaves if leaf[0] not in leaves_to_remove]

        # print(leaves_to_remove)
        # print(other_leaves)

        result = ','.join(other_leaves)

        return result

    def make_gene_family(self, genome_leaves):
        """ Creates Gene Families and writes them to a text file

        For each genome, read through all the pairwise SynMap outputs. 
        First, extract information about each gene in the outputs using extract_gene_info(). 
        Then, if the similarity between the genes in the SynMap outputs falls within
        the threshold, update the gene_list and gene_family dictionaries. 
        Finally, modify the format so that the gene_family_labels are in consecutive order 
        and output these two dictionaries into files.
        
        Generate gene family file of columns in the order of "geneName","geneFamilyID","chr", "start","end",“strand","genomeCoGeID”. 

        :param genome_leaves: list
            list containing information for all the genome leaves
        """

        for i in range(0, len(genome_leaves)):
            for j in range(i, len(genome_leaves)):
                leaf_index_1 = genome_leaves[i][0]
                leaf_index_2 = genome_leaves[j][0]
                file_name = self.synteny_file_path + leaf_index_1 + "_" + leaf_index_2 + self.synteny_file_ending
                # print(file_name)

                with open(file_name) as file:
                    lines = file.readlines()

                    for line in lines:
                        if line[0] != "#":

                            gene_name_1, gene_info_1, gene_name_2, gene_info_2, ks, kn, similarity = \
                                self.extract_gene_info(line, leaf_index_1, leaf_index_2)

                            if (self.min_cutoff_weight < similarity) and (self.max_cutoff_weight >= similarity):
                                self.gene_pair.append([gene_name_1, gene_name_2, similarity, ks])
                                self.update_genes(gene_name_1, gene_info_1, leaf_index_1, gene_name_2, gene_info_2,
                                                  leaf_index_2)

        self.format_genes()


        with open(self.gene_family_output, 'w') as f:
            for family, group_of_genes in self.gene_family.items():
                f.write(str(family) + " ")

                for gene in group_of_genes:
                    f.write("%s " % gene)

                f.write("\n")

        with open(self.gene_list_output, 'w') as f:
            for gene, gene_info in sorted(self.gene_list.items(),key=self.get_GeneFamily_ID):
                f.write("%s " % gene_info[1][3])
                f.write("%s " % str(gene_info[2]))
                f.write("%s " % gene_info[1][0])
                f.write("%s " % gene_info[1][1])
                f.write("%s " % gene_info[1][2])
                f.write("%s " % gene_info[1][4])
                f.write("%s " % gene_info[0])
                f.write("\n")
                
    def get_GeneFamily_ID(self,item):
        return item[1][2]
    
    
    def update_genes(self, gene_name_1, gene_info_1, leaf_index_1, gene_name_2, gene_info_2, leaf_index_2):
        """Updates the gene_list and gene_family dictionaries accordingly

        First, determine whether each gene is present in the gene_list dictionary. 
        If they both aren't, then add both to the gene_list, and add them together into 
        the same gene_family. If one of them is present and the other isn't, add the one 
        that isn't present to the gene_list, and add it to the present genes' gene_family_label.
        Finally, if both are already present, remove one of the gene_family_labels from 
        the gene_family dictionary, and update that gene's gene_family_label value in gene_list. 
        Finally, add the gene to the other one's gene_family_label key in the gene_family dictionary.

        :param gene_name_1: str

        :param gene_info_1: list

        :param leaf_index_1: str
            leaf index for the gene related to gene_name_1

        :param gene_name_2: str

        :param gene_info_2: list

        :param leaf_index_2: str
            leaf index for the gene related to gene_name_2
        """

        gene1_in_genelist = gene_name_1 in self.gene_list
        gene2_in_genelist = gene_name_2 in self.gene_list

        if (not gene1_in_genelist) and (not gene2_in_genelist):
            self.gene_list[gene_name_1] = [leaf_index_1, gene_info_1, self.gene_family_label]
            self.gene_list[gene_name_2] = [leaf_index_2, gene_info_2, self.gene_family_label]
            self.gene_family[self.gene_family_label] = [gene_name_1, gene_name_2]
            self.gene_family_label += 1

        elif gene1_in_genelist and (not gene2_in_genelist):
            current_gene_family_1 = self.gene_list[gene_name_1][2]
            self.gene_list[gene_name_2] = [leaf_index_2, gene_info_2, current_gene_family_1]
            self.gene_family[current_gene_family_1].append(gene_name_2)

        elif (not gene1_in_genelist) and gene2_in_genelist:
            current_gene_family_2 = self.gene_list[gene_name_2][2]
            self.gene_list[gene_name_1] = [leaf_index_1, gene_info_1, current_gene_family_2]
            self.gene_family[current_gene_family_2].append(gene_name_1)

        else:
            current_gene_family_1 = self.gene_list[gene_name_1][2]
            current_gene_family_2 = self.gene_list[gene_name_2][2]

            if current_gene_family_1 != current_gene_family_2:
                genes_2 = self.gene_family[current_gene_family_2]

                for gene in genes_2:
                    self.gene_list[gene][2] = current_gene_family_1

                self.gene_family[current_gene_family_1].extend(genes_2)

                del self.gene_family[current_gene_family_2]

    # def check_genome_id(self, list_of_genes):
    #     result = []
    #
    #     for gene in list_of_genes:
    #         genome_id = self.gene_list[gene][1][0]
    #
    #         if genome_id not in result:
    #             result.append(genome_id)
    #
    #     return result

    def format_genes(self):
        """Updates the gene_family_label field for all genes in the gene_list and reformat the gene_family dictionary

        First, remove any Gene Families that have lengths over 500 and remove any genes in those families from the gene_list. 
        Then, create a new Gene Family dictionary with ordered, consecutive gene_family_labels.
        Finally, update the existing gene_list using the list of genes in each gene family.
        """

        current_index = 1
        formatted_genefamily = {}

        ids_to_remove = []
        families_to_remove = []

        for family, members in self.gene_family.items():
            genomes_seen = {}
            genomes_seen_2 = {}

            # gf1 implementation
            if len(members) > int(self.max_genes_in_family):
                # print(len(value))
                families_to_remove.append(family)

                for gene_id in members:
                    ids_to_remove.append(gene_id)

                pass

            # gf2 implementation
            if len(members) > self.max_genes_from_genome:
                for member in members:
                    member_genome_id = self.gene_list[member][0]

                    if member_genome_id not in genomes_seen.keys():
                        genomes_seen[member_genome_id] = 1

                    else:
                        genomes_seen[member_genome_id] += 1

                    if genomes_seen[member_genome_id] > int(self.max_genes_from_genome):
                        families_to_remove.append(family)

                        for gene_id in members:
                            ids_to_remove.append(gene_id)

                        break

                pass

            # gf3 implementation
            for member in members:
                member_genome_id_2 = self.gene_list[member][0]

                genomes_seen_2[member_genome_id_2] = 1

                if len(genomes_seen_2) >= self.min_genomes_in_family:
                    break

            if len(genomes_seen_2) < self.min_genomes_in_family:
                families_to_remove.append(family)

                for gene_id in members:
                    ids_to_remove.append(gene_id)

        # print(len(keys_to_remove))

        for family in list(set(families_to_remove)):
            # print(key, self.gene_family[key])
            del self.gene_family[family]

        for ids in list(set(ids_to_remove)):
            del self.gene_list[ids]

        for family, members in self.gene_family.items():
            formatted_genefamily[current_index] = members
            current_index += 1

        for family, members in formatted_genefamily.items():
            for gene in members:
                self.gene_list[gene][-1] = family

        self.gene_family = formatted_genefamily.copy()

    def extract_gene_info(self, line, leaf_index_1, leaf_index_2):
        """Extracts information from the pairwise SynMap output

        :param line: list
            a line of information from the pairwise SynMap output

        :param leaf_index_1: int
            index of the first leaf

        :param leaf_index_2: int
            index of the second leaf

        :return: gene_name_1: str
                 gene_info_1: list
                 gene_name_2: str
                 gene_info_2: list
                 ks: int
                 kn: int
                 similarity: int
        """
        info = line.split()

        gene_info_1 = info[3].split('||')
        gene_info_2 = info[7].split('||')

        gene_name_1 = leaf_index_1 + gene_info_1[3]
        gene_name_2 = leaf_index_2 + gene_info_2[3]

        ks = info[0]
        kn = info[1]
        similarity = gene_info_1[-1]

        return gene_name_1, gene_info_1, gene_name_2, gene_info_2, ks, kn, float(similarity)
