class Contig:
    """
    A class used to create Contigs

    Methods
    ----------
    get_edge()
        read in the output from the Maximum Weight Matching algorithm

    find_telomeres()
        finds and creates a list of all telomeres

    find_node_1(node)
        given a node of a gene, find the other node

    get_telomere_contigs(telomere_start)
        constructs a contig from telomere_start

    get_contigs(telomere_list)
        finds and creates a list of all contigs and genes inside them
    """

    def __init__(self, mwm_output, num_gene_families):
        """
        :param mwm_output: file
            Maximum Weight Matching Output

        :param num_gene_families: int
        """

        self.mwm_output = mwm_output
        self.num_gene_families = num_gene_families

    def get_edge(self):
        """read in the output from the Maximum Weight Matching algorithm and format it as a dictionary

        edges: dictionary
            key: int
                each key corresponds to a node of a gene

            value: list
                each value is a list consisting of [node_of_neighbour_gene, weight]
        """

        edges = {}

        with open(self.mwm_output) as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            lines[i] = line.rstrip()
            lines[i] = line.split()
            edges[int(lines[i][0])] = [int(lines[i][1]), int(lines[i][2])]
            edges[int(lines[i][1])] = [int(lines[i][0]), int(lines[i][2])]

            # for old data testing purposes
            # edges[int(lines[i][0]) + 1] = [int(lines[i][1]) + 1, int(lines[i][2])]
            # edges[int(lines[i][1]) + 1] = [int(lines[i][0]) + 1, int(lines[i][2])]

        self.mwm_output = edges

        # print(list(self.mwm_output.keys())[0])
        # print(list(self.mwm_output.values())[0])
        #
        # print(list(self.mwm_output.keys())[-1])
        # print(list(self.mwm_output.values())[-1])
        #
        # print(len(self.mwm_output))

        # print(lines)
        # self.mwm_output = lines

    def find_telomeres(self):
        """finds and creates a list of all telomeres

        telomeres_list:
            If a node isn't present in the Maximum Weight Matching output, then it is added as a telomere.

        :return: telomeres_list
            a list with all the telomeres
        """

        telomeres_list = []

        for i in range(1, (2 * self.num_gene_families) + 1):
            if i not in self.mwm_output:
                telomeres_list.append(i)

        return telomeres_list

    def find_node_1(self, node):
        """given a node of a gene, find the other node

        Each gene has two nodes, when given one node, this method finds the other.

        :param node: str
            input node of a gene

        :return: other_node: int
            the other node

        :return: gene_value: int
            id of the gene

        :return: orientation: int
            orientation of the gene
        """

        int_node = int(node)

        if int_node % 2 == 1:
            return int_node + 1, int((int_node + 1) / 2), 1

        return int_node - 1, int(int_node / 2), int(-1)

    def get_telomere_contigs(self, telomere_start):
        """constructs a contig from telomere_start

        From the start of the contig, track the adjacencies between genes while appending the genes into the contig
        until the ending telomere is reached. As the adjacencies are tracked, edges are deleted from self.mwm_output.

        :param telomere_start: string
            start of the contig

        :return: contig: list
            list of genes in a contig

        :return: node1: int
            end of the contig
        """

        contig = []

        # print("Contig start Tele ", telomere_start)

        node1, gene_value, orientation = self.find_node_1(telomere_start)
        contig.append(gene_value * orientation)

        # print("first gene node1/gene/sign  ", node1, gene_value, orientation)
        # print("contig now ", contig)
        while node1 in self.mwm_output:
            node2 = self.mwm_output[node1][0]
            # print("DEL end of previous gene  ", node1)
            del self.mwm_output[node1]
            node1, gene_value, orientation = self.find_node_1(node2)
            # print("next gene startNode/endNode/gene/sign  ",node2, node1, gene_value, orientation)
            contig.append(gene_value * orientation)
            # print("contig now ", contig)
            # node1, gene_value, orientation = self.find_node_1(self.mwm_output[node1][0])
            # print("DEL start of this gene ", node2)

            del self.mwm_output[node2]

        return contig, node1

    def get_contigs(self, telomeres_list, output_file, anc):
        """finds and creates a list of all contigs and genes inside them

        First, create the contigs that can be found between telomeres in telomeres_list by iterating get_telomere_contigs().
        Then, if there are still edges in self.mwm_output, cut the edge with the lowest weight into two nodes, which
        are treated as telomeres. Repeat this process until there are no more edges in self.mwm_output.

        :param output_file: str
            path to output file

        :param telomeres_list: list
            list of all telomeres

        :param anc: str
            specifies the ancestor

        :return: contigs_list
            list of all contigs
        """

        contigs_list = []
        cut_index = 0

        while telomeres_list:
            start = telomeres_list.pop(0)
            #print(start)
            contig, node1 = self.get_telomere_contigs(start)
            #print("remove end Tel  ", node1)
            telomeres_list.remove(node1)
            contigs_list.append(contig)
        
        while len(self.mwm_output) != 0:
            # print(self.mwm_output)
            min_weight = min(self.mwm_output.items(), key=lambda x: x[-1][-1])
            # print(min_weight)
            del self.mwm_output[min_weight[0]]
            del self.mwm_output[min_weight[1][0]]
            start = min_weight[0]
            contig, node1 = self.get_telomere_contigs(start)
            assert node1 == min_weight[1][0]
            contigs_list.append(contig)
            cut_index += 1

        # print(cut_index)

        with open(output_file, "w") as f:
            contigs_list.sort(key=lambda x: (len(x), x[0]), reverse=True)
            counter = 1

            f.write("> Ancestor " + str(anc) + "\n")

            for contig in contigs_list:
                # f.write("contig " + str(counter) + "\n")
                for gene in contig:
                    f.write("%s" % gene + " ")

                counter += 1
                f.write(" | \n")
        
        return contigs_list

