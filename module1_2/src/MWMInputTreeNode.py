class MWMInputTreeNode:
    """
        A class used to generate input for the Maximum Weight Matching algorithm

        Methods
        ----------
        get_mwm_input(a_median_structure, num_gene_families, window_size)
            applies methods to get the Maximum Weight Matching input

        combine_edges(edge1, edge2, edge3)
            sums the weights of all edges

        modify_weight(edge, weight)
            modify the weight of an edge by a specified amount

        get_edge(genome_group, window_size)
            gets edge for all genomes in a branch
        """

    def __init__(self, genomes, all_leaves):
        """
        :param genomes: file with chromosomes and their respective genes
        :param all_leaves: file containing information of the initial input genomes
        """

        self.genomes = genomes
        self.all_leaves = all_leaves

#
#     def get_tree_structure(self):
#         """read in the parameters of a file to obtain the tree structure
#
#         For each ancestor, information of its 3 branches are obtained, along with other information.
#
#         :return: median_structure: list
#             list with all 3 branches for each ancestor
#
#         :return: internal_node_number: int
#             number of ancestors
#
#         :return: window_size: int
#         :return: num_genomes: int
#         :return: num_gene_families: int
#         """
#
#         median_structure = []
#
#         with open(self.tree_node_file) as f:
#             tree_node_input = f.readlines()
#
# #? is there a way to combine this tree structure file into genome.txt?
# #? num_gene_families should be from previous constructing gene family steps
# #? in previous step, you've got 9048 gene families, why here's hard coded as 10278?
# #? how is the median structure defined?
#         internal_node_number, window_size, num_genomes, num_gene_families = tree_node_input[1].split("\t")
#
#         # print(internal_node_number, window_size, num_genomes, num_gene_families)
#
#         for line in tree_node_input[2:]:
#             if line[0] != "#":
#                 median_structure.append(line.rstrip().split("\t"))
#
#         # print(median_structure)
#
#         return median_structure, internal_node_number, window_size, num_genomes, num_gene_families
#

    def get_mwm_input(self, a_median_structure, num_gene_families, window_size, min_adj_weight, output_file):
        """applies methods to get the Maximum Weight Matching input

        First, the median structure is split up into 3 groups, each representing a descendent of the ancestor.
        Then the adjacencies for each group is calculated with get_edge().
        The adjacencies for each group is then modified and combined with modify_weight() and
        combine_edges(), resulting in the Maximum Weight Matching input.

        :param a_median_structure: list
            list with all 3 branches for one ancestor

        :param num_gene_families: int
        :param window_size: int
        :param min_adj_weight: int

        :param output_file: string
            path to output file

        :return: mwm_input: list: [node1, node2, weight]
            input for maximum weight matching
        """
     
        mwm_input = []

        group1_genome = a_median_structure[0].split(",")
        group2_genome = a_median_structure[1].split(",")
        # print(a_median_structure[2])
        group3_genome = a_median_structure[2].split(",")

        # print("a_median_structure")
        # print(group1_genome)
        # print(group2_genome)
        # print(group3_genome)

        edge_in_genome_1 = self.get_edge(group1_genome, window_size, num_gene_families)
        edge_in_genome_2 = self.get_edge(group2_genome, window_size, num_gene_families)
        edge_in_genome_3 = self.get_edge(group3_genome, window_size, num_gene_families)

        all_edges = self.combine_edges(edge_in_genome_1, edge_in_genome_2, edge_in_genome_3)

        for node1, node2 in all_edges.items():
            for node, weight in node2.items():
                if weight > min_adj_weight:
                    mwm_input.append([node1, node, weight])

        with open(output_file, 'w') as f:
            for adjacency in mwm_input:
                for item in adjacency:
                    f.write("%s" % item + " ")
                f.write("\n")

        return mwm_input

    def combine_edges(self, edge1, edge2, edge3):
        """combines all edges, sums the weight if the adjacency is present in all edges.

        First, modifies the weight of each edge according to a specified amount. Then finds the union of all edges.

        :param edge1: dict
        :param edge2: dict
        :param edge3: dict

        :return: result: dict
            union of all 3 edges
        """

        result = {}

        edge1 = self.modify_weight(edge1, 100)
        edge2 = self.modify_weight(edge2, 100)
        edge3 = self.modify_weight(edge3, 100)

        # print(edge3)

        all_nodes = set(edge1) | set(edge2) | set(edge3)

        for node in all_nodes:
            edge1.setdefault(node, {})
            edge2.setdefault(node, {})
            edge3.setdefault(node, {})

            combined_edges = {k: edge1[node].get(k, 0) + edge2[node].get(k, 0) + edge3[node].get(k, 0) for k in set(edge1[node]) | set(edge2[node])| set(edge3[node])}
            result[node] = combined_edges

        return result

    def modify_weight(self, edge, weight):
        """modify the weight of an edge by a specified amount

        :param edge: dict
        :param weight: int

        :return: edge: dict
        """

        for node in edge.values():
            for key, weights in node.items():
                if weights > 0:
                    node[key] += weight

        return edge

    def get_edge(self, genome_group, window_size, num_gene_families):
        """gets edge for all genomes in a branch

        For each genome in genome_group, for each chromosome in the genome, get all the adjacencies between genes
        according to the window size. The weight of the edge is how many times the adjacency is present.

        :param genome_group: list
            genomes in a branch

        :param window_size: int

        :return: adj: dict of dict
            key: node1
            value: {node2: weight}
        """

        adj = {}

        for genome_id in genome_group:
            for chromosome in self.genomes[genome_id]:

                for i in range(0, len(chromosome)):
                    gene1 = chromosome[i]

                    gene1_value = int(abs(gene1))

                    # gene1_node1 = 2 * gene1_value - 1
                    gene1_node2 = 2 * gene1_value

                    if gene1_value != gene1:
                        # gene1_node1 = 2 * gene1_value
                        gene1_node2 = 2 * gene1_value - 1

                    for j in range(1, int(window_size)):
                        if (i + j) >= len(chromosome):
                            break

                        gene2 = chromosome[i + j]

                        gene2_value = int(abs(gene2))

                        gene2_node1 = 2 * gene2_value - 1
                        # gene2_node2 = 2 * gene2_value

                        if gene2_value != gene2:
                            gene2_node1 = 2 * gene2_value
                            # gene2_node2 = 2 * gene2_value - 1

                        # print(gene1_node2, gene2_node1)

                        node1 = min(gene1_node2, gene2_node1)
                        node2 = max(gene1_node2, gene2_node1)

                        if (node1 != node2) and (node1 != (node2 - 1)) and (node1 % 2 == 1):
                            if node1 in adj:
                                if node2 in adj[node1]:
                                    adj[node1][node2] += 1
                                else:
                                    adj[node1][node2] = 1

                            else:
                                adj[node1] = {node2: 1}

        return adj
