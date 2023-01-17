
# raccroche-v2

This repository contains the code for Module 1 and 2 of the RACCROCHE pipeline.


## Overview
Given the phylogenetic relationships of several extant species, the reconstruction of their ancestral genomes at the gene and chromosome level is made difficult by the cycles of whole genome doubling followed by fractionation in plant lineages. Fractionation scrambles the gene adjacencies that enable existing reconstruction methods. We propose an alternative approach that postpones the selection of gene adjacencies for reconstructing small ancestral segments and instead accumulates a very large number of syntenically validated candidate adjacencies to produce long ancestral contigs through maximum weight matching. Likewise, we do not construct chromosomes by successively piecing together contigs into larger segments, but instead count all contig co-occurrences on the input genomes and cluster these, so that chromosomal assemblies of contigs all emerge naturally ordered at each ancestral node of the phylogeny. These strategies result in substantially more complete reconstructions than existing methods. We deploy a number of quality measures: contig lengths, continuity of contig structure on successive ancestors, coverage of the reconstruction on the input genomes, and rearrangement implications of the chromosomal structures obtained. The reconstructed ancestors can be functionally annotated and are visualized by painting the ancestral projections on the descendant genomes, and by highlighting syntenic ancestor-descendant relationships. Our methods can be applied to genomes drawn from a broad range of clades or orders.

### Module 1: construct gene families and list candidate adjacencies
- Step 1: Pre-process gene families
- Step 2: List generalized adjacencies
- Step 3: List candidate adjacencies

### Module 2: construct ancestral contigs through Maximum Weight Matching
- Step 4: Construct contigs through Maximum Weight Matching
## Run Locally

Clone the project

```bash
  git clone https://git.cs.usask.ca/nma904/module1_2.git
```

Go to the project directory

```bash
  cd module1_2
```

Add input data to the inputData folder

Edit config.yaml to desired specifications

### Run Modules Through Command Line

Run the module1_main.py file
```bash
  python module1_main.py --config path_to_config_file
```

### Run Modules Using Module1_2.ipynb

Run each cell block in the Jupyter Notebook.





## Authors

- [Qiaoji Xu](https://github.com/Qiaojilim)
- [Lingling Jin](https://github.com/jin-repo/RACCROCHE)
- Chunfang Zheng
- James H. Leeben-Mack
- David Sankoff
- [Alex Liu](https://github.com/alexl0806)

## Emails
- limqiaojixu@gmail.com
- lingling.jin@cs.usask.ca


## License

[BSD](https://en.wikipedia.org/wiki/BSD_licenses/)


## Citations

If you use the raccroche pipline for ancestral contigs reconstruction, please cite:
[RACCROCHE: ancestral flowering plantchromosomes and gene ordersbased on generalized adjacenciesand chromosomal gene co-occurrences](https://link.springer.com/chapter/10.1007/978-3-030-79290-9_9)
