## Sequence Annotation Pipeline

### Description

A modular annotation pipeline for prokaryotic nucleotide sequences.

### Local installation

You can run it as a docker container, providing you have [installed](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce) the docker first:

Pull the latest version of the image:

```docker pull nestorzaburannyi/annotate```

Then run the image as a container in an interactive mode:

```docker run -it nestorzaburannyi/annotate```

Or, also mount your working directory as /pwd within the container.

```docker run -it -v "$(pwd)":/pwd nestorzaburannyi/annotate```

### Usage

Update the databases before the first use:

```annotate -update```

Start the annotation

```Usage: annotate -input <input_file> [OPTIONS]
Example: annotate -input my_genome.fasta -taxid 13

        Parameters              Type            Default                                 Description

Mandatory parameters:
        -input                  string          -                                       input file name (FASTA or GenBank format)

Optional parameters:
        -rna                    boolean         true                                    prediction of RNA genes
        -rna-r                  boolean         true                                    prediction of ribosomal RNA genes
        -rna-t                  boolean         true                                    prediction of transport RNA genes
        -rna-tm                 boolean         true                                    prediction of transport-messenger RNA genes
        -rna-nc                 boolean         true                                    prediction of noncoding RNA genes
        -cds                    boolean         true                                    prediction of CDS genes
        -cds-i                  boolean         true                                    prediction of CDS genes (ab initio)
        -cds-h                  boolean         true                                    prediction of CDS genes (homology-based)
        -cds-a                  boolean         true                                    annotation of CDS genes
        -taxid                  integer         -                                       NCBI taxonomy id of the organism (e.g. 52)
        -strain                 string          -                                       strain name
        -bioproject             string          -                                       a BioProject identifier (NCBI-compatibility)
        -prefix                 string          -                                       prefix for locus numbering (NCBI-compatibility)
        -step                   integer         10                                      increment of locus numbering (NCBI-compatibility)
        -sbt                    string          -                                       file with a filled sbt form
        -transparent            boolean         false                                   transparent mode allows to keep already annotated features
        -help                   boolean         false                                   print this help
        -quiet                  boolean         false                                   suppress messages to the console
        -verbose                boolean         false                                   print more messages to the console
        -output                 string          current working directory               write output files to: output/UUID

Advanced parameters:
        -rna-r-program          string          rnammer                                 rnammer, infernal
        -rna-t-program          string          trnascanse                              trnascanse, aragorn, infernal
        -rna-tm-program         string          aragorn                                 aragorn, infernal
        -rna-nc-program         string          infernal                                infernal
        -cds-i-program          string          prodigal                                prodigal, glimmer, genemarks
        -cds-i-score            string          5                                       score threshold for ab initio predictions
        -cds-h-program          string          blast                                   blast, diamond
        -cds-h-evalue           string          1e-5                                    e-value threshold for homology searches
        -cds-a-program          string          pannzer                                 pannzer, emapper
```
