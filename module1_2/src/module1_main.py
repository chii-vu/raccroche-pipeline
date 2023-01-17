"""module1_main

This script allows the user to parse and load a configuration file for the pipeline.

This tool accepts .yaml files.

This script requires that `yaml` be installed within the user's Python environment.
"""


# import statements
import argparse
import sys
import yaml
from GeneFamily import GeneFamily
from Genome import Genome
from MWMInputTreeNode import MWMInputTreeNode
from mwmatching import *
from Contig import *
import os


def parse_input(argv):
    """Parses input files

    :param argv:
    :return parser.parse_args():
    """

    parser = argparse.ArgumentParser(description="")

    # command-line option for configuration file
    parser.add_argument('--config', help="", dest="config_file", required=True)

    return parser.parse_args()


def load_config(config_file):
    """Loads the configuration file

    :param config_file:
    :return loaded_config:
    """

    try:
        with open(config_file, 'r') as stream:
            loaded_config = yaml.load(stream, Loader=yaml.FullLoader)

    except Exception as msg:
        print("Error While Loading Configuration File: {}".format(msg))
        sys.exit(0)

    return loaded_config


def main():
    """Executes functions to parse and load a configuration file

    :param input_args:
    """

    if len(sys.argv) < 3:
        print('Usage: module1_main.py --config path_to_config_file')
        sys.exit(0)

    parsed_input = parse_input(sys.argv)
    parsed_config = load_config(parsed_input.config_file)

    # print(data_file)
    # print(leaf_file)

    print("Start")

    directory = parsed_config["output_path"] + parsed_config["project_name"]
    os.makedirs(directory, exist_ok=True)

    # STEP 1: Reading in Genomes and Newick Tree Structure
    print("STEP 1: Reading in Genomes and Newick Tree Structure")
    gene_family = GeneFamily(parsed_config)
    all_leaves, median_structure, newick_structure = gene_family.get_leaves_and_tree_info()
    # print(all_leaves)
    # print(median_structure)
    print("Step 1 DONE")

    # STEP 2: Creating Gene Families
    print("STEP 2: Creating Gene Families")
    gene_family.make_gene_family(all_leaves)
    print("Step 2: DONE")

    # STEP 3: Representing Genomes by Gene Family Labels
    print("STEP 3: Representing Genomes by Gene Family Labels")
    initialize_genome = Genome(gene_family.gene_list, all_leaves, parsed_config)
    genome = initialize_genome.get_genome_in_string()
    print("Step 4: DONE")

    # STEPS 4-6: Generating Maximum Weight Matching Input, Performing MWMatching, Constructing Ancestral Contigs
    print("STEPS 4-6: Generating Maximum Weight Matching Input, Performing MWMatching, Constructing Ancestral Contigs")

    # Parameters for Generating Maximum Weight Matching Inputs
    window_size = parsed_config['ws']
    min_adj_weight = parsed_config['min_mwm_weight']
    num_gene_families = len(gene_family.gene_family)

    # num_genomes = parsed_config['gn']
    # tree_node = num_genomes - 2
    # print(num_genomes, tree_node, window_size, num_gene_families)

    input_tree_node = MWMInputTreeNode(genome, all_leaves)

    anc = 1
    dcj_files = []

    for structure in median_structure:
        # print(structure)

        # output file paths for Maximum Weight Matching Input, Maximum Weight Matching Output and Ancestral Contigs
        outfile_mwmin = directory + parsed_config["output_file"]["mwm_input_template"] + str(anc) + ".txt"
        outfile_mwmout = directory + parsed_config["output_file"]["mwm_output_template"] + str(anc) + ".txt"
        outfile_contig = directory + parsed_config["output_file"]["contig_template"] + str(anc) + ".txt"

        # Adding File Paths for later DCJ Computations
        dcj_files.append(outfile_contig)

        # Generating Maximum Weight Matching Input
        print("Starting MWM Input Generation")
        mwm_input = input_tree_node.get_mwm_input(structure, num_gene_families, window_size, min_adj_weight, outfile_mwmin)
        # print("MWM Input Generation DONE")

        # Performing Maximum Weight Matching
        print("Starting Maximum Weight Matching")
        maxWeightMatching(mwm_input, outfile_mwmout)
        # print("Maximum Weight Matching DONE")

        # Constructing Ancestral Contigs
        print("Starting Construction of Ancestral Contig")
        contig = Contig(outfile_mwmout, num_gene_families)
        contig.get_edge()
        list_telomeres = contig.find_telomeres()
        contig_list = contig.get_contigs(list_telomeres, outfile_contig, anc)
        # print("Contig Construction DONE")

        anc += 1

    # STEP 7: Calculating DCJ Distance Between Ancestral Contigs
    print("STEP 7: Calculating DCJ Distance Between Ancestral Contigs")
    jar_path = parsed_config['input_file']['jar_path']
    dcj_output = []

    for i in range(1, anc):
        for j in range(i + 1, anc):
            dcj_output_path = directory + parsed_config['output_file']['dcj_output_path'] + str(i) + "_" + str(j) + ".txt"
            command = "java -jar " + jar_path + "/UniMoG-java11.jar " + str(dcj_files[i - 1]) + " " + str(dcj_files[j - 1]) + " -m=1 -p >" + dcj_output_path
            dcj_output.append(dcj_output_path)
            os.system(command)

    print("Step 7 DONE")

    dcj_summary = directory + parsed_config['output_file']['dcj_summary_path']

    print("Generating Summary of DCJ Calculations")
    with open(dcj_summary, 'w') as dcj_summary_file:
        for i in range(len(median_structure)):
            dcj_summary_file.write("median structure for Ancestor  "+str((i + 1))+":"+ "%s" % median_structure[i] + "\n")

        for file in dcj_output:
            path = file

            with open(path, 'r') as dcj_file:
                dcj_info = dcj_file.readlines()

            dcj_summary_file.write(dcj_info[0])

    print("Summary Generated")

    # print(java_command)
    # print(dcj_files)


if __name__ == '__main__':
    main()