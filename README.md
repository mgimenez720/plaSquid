# plaSquid: PLAsmid SeQuences IDentification in metagenomic assemblies.  

### Description

- **plaSquid** is a Nextflow pipeline for plasmid detection and classification in genomic and metagenomic data. This pipeline uses either genomic or metagenomic *assemblies* as input. It uses two different approaches to look for plasmids: to search against a plasmidic database and to look for plasmid specific genes (i.e. RIP and *Relaxases*).

- **plaSquid** also classifies plasmids into replicon types and MOB groups by comparing RIPs or Relaxases against custom HMMs. 

- It summarises the information gathered by the two complementary approaches in a single output table and allows further analysis on plasmidic contigs as ir outputs plasmidic contigs in a multifasta file (*"Result.fasta"*)

### Installation

    git clone https://github.com/mgimenez720/plaSquid/
    cd plaSquid/
 
Using the option (*"-profile conda"*) when running plaSquid will build a conda environment within the base directory. This environment can be reused in subsequent runs.   

#### Dependencies




### Usage 


    nextflow run plaSquid.nf --contigs 'data/*.fasta'

    arguments:

    --contigs       Path to input data (must be surrounded with quotes).
    --mmi           Path to Minimap2 indexed (.mmi) or fasta (.fasta/.fna) plsdb database.
    --outdir        Path to output directory where results are written.
    --help          Print help message and exit

    subworkflows:

    --minidist      Run mapping of contigs against plsdb database.
    --repsearch     Run search and classification of RIP and MOB (Rel) genes.

    profiles:

    -profile conda  Installs dependencies using a conda environment
    -profile server runs using 15 cpus and 50 Gb
    -profile test   tests dependencies and normal functioning


    Authors:

    Matías Giménez
    Ignacio Ferrés
    Gregorio Iraola


    Microbial Genomics Laboratory
    Institut Pasteur Montevideo (Uruguay)

