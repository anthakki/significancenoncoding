
# Noncoding variant significance analysis 

This repository contains the updated source code for our noncoding variant discovery analysis. The analysis is based on the work by [Dietlein et al.](https://doi.org/10.1126/science.abg5601) but has been extended to be better suited for large-scale datasets and be more robust to dataset variations, e.g. by more detailed modeling of variant correlations. The original source code was openly provided by [Dietlein et al.](https://doi.org/10.5281/zenodo.5913867) but has been heavily modified since. If you use the modified pipeline, please cite our [work](#) detailing the newer analysis.

# Installation 

Pre-built [Docker images](https://hub.docker.com/r/anthakki/significancenoncoding) are available along with the [Docker build files](https://github.com/anthakki/significancenoncoding-docker) and instructions, which takes care of the dependencies and will download the resources for offline use. The latest Docker image can be fetched from Docker Hub with (this will download ~3.4 GB of files, typically taking few to several minutes, but can varying significantly):
```bash
	docker pull anthakki/significancenoncoding
```

If you prefer to compile the program yourself, e.g. for developement purposes, `make` will build the Java binaries. This requires a [Java JDK 8](https://jdk.java.net/java-se-ri/8-MR6) or newer in `$PATH` and [GNU make](https://www.gnu.org/software/make/). You can simply issue (takes few seconds):
```bash
	make
```

If you have [cURL](https://curl.se/), the runtime resources can be automatically downloaded using (again, will download ~3.4 GB of files, which can take several minutes):
```bash
	make AnnotationFilesComplete.zip
```

For development, OpenJDK 1.8.0_452, GNU make 3.81, and cURL 8.7.1 on macOS 14 arm64 have been used. The Docker images have been built using OpenJDK 1.8.0_452, GNU make 4.3, and cURL 8.5.0 on Ubuntu Linux 24.04 amd64 with podman 4.9.3 (substituting for Docker).

A POSIX shell, [Python 3](https://www.python.org/), and [AWS CLI](https://aws.amazon.com/cli/) are required to run the demo below. Python 3.9.6 and AWS CLI 2.27.17 were used on macOS 14 arm64.

# Usage

A [user manual](https://zenodo.org/records/5913867/preview/UserManual.pdf) for the original method is available. This very much still applies for the inputs and outputs of the method as well running it.

In addition, a wrapper shell script [SignificanceNoncoding](SignificanceNoncoding.1.md) is provided, which accepts both analysis and Java options, and will automatically set the Java classpath and locate the support files if they are already downloaded in the same path. For example, following the manual:
```bash
	./SignificanceNoncoding -Xmx55G Breast ~/MutationTestFiles/MutationFiles.txt ~/SignificanceAnalysis/
```
will configure Java with 55 GB of heap, and run the analysis for the entity `Breast` using the input files specified in `~/MutationTestFiles/MutationFiles.txt` and output the results in `~/SignificanceAnalysis/`. The available options are lxhaustively documented in its [man page](SignificanceNoncoding.1.md). The actual amount of heap memory required varies by the dataset, but ranges typically from 10 to 100 GB.

When running through Docker, the script is the entry point, so simply run:
```bash
	docker run [docker-options] anthakki/significancenoncoding -Xmx55G Breast ~/MutationTestFiles/MutationFiles.txt ~/SignificanceAnalysis/
```

Note that Docker options must be set to allow access to the data and output files from within the container, depending where they are located.

Please refer to the [manual](https://zenodo.org/records/5913867/preview/UserManual.pdf) how to format your input files and interpret the outputs.

Briefly, the mutations are arranged by cancer type in a variant of tab-separated text [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files, specifying their coordinate (`Chromosome`, `Position`), the alteration (`Reference_Allele`, `Tumor_Seq_Allele`), and the involved sample (`Tumor_Sample_Barcode`). The reference files are in hg19, as is assumed for the mutation coordinates as well. Chromosomes `1` to `22`, and `X` and `Y` are used, and listed without the `chr` prefix. Positions are 1-based, as in VCF files, and refer to the position of the leftmost reference nucleotide. The variants should be listed in a VCF-style flanked format for indels, such that the reference and alternative alleles are never empty. The sample barcodes are not interpreted, only checked for equality. An example follows:
 
| Chromosome | Position  | Reference_Allele | Tumor_Seq_Allele | Tumor_Sample_Barcode |
|:-----------|:----------|:-----------------|:-----------------|:---------------------|
| 1          |    565949 | C                | T                | DO218440             |
| 4          |  35520715 | G                | GA               | DO218440             |
| 6          |   2212887 | CCTG             | C                | DO218695             |
| 2          | 146790325 | C                | T                | DO218695             |
| 13         |  83189877 | TA               | T                | DO218695             |
| ...        | ...       | ...              | ...              | ...                  |

A tab-separated text file listing the corresponding entity names and the paths to these MAF files is used used to refer to the input data collectively (e.g. `~/MutationTestFiles/MutationFiles.txt` in above). There is no header in this file (the headers in cursive are merely illustrative) and any relative paths are relative to the file itself.

| _entity_   | _filename_        |
|:-----------|:------------------|
| Biliary    | Biliary.maf.gz    |
| Bladder    | Bladder.maf.gz    |
| Brain      | Brain.maf.gz      |
| ...        | ...               |

If the entity name starts with `Biliary`, `Bladder`, `Brain`, `Breast`, `Cervix`, `Colorectal`, `Endocrine`, `Endometrium`, `Esophagus`, `Gastric`, `HeadNeck`, `Kidney`, `Leukemia`, `Liver`, `Lung`, `Lymphoma`, `Myeloid`, `Ovary`, `Pancreas`, `Pleura`, `Prostate`, `Sarcoma`, `Skin`, or `Thyroid`, the corresponding expression data are used for annotating the significant regions with nearby genes.

Finally, the output files are written into the output directory (e.g. `~/SignificanceAnalysis/` in the above), with one `FDR_Weighted_Combined_<entity>.txt` per entity. This file lists the genomic coordinate (chromosome, a BED-style 0-based start position of the 2,500 kb genomic window), the final pooled FDR combined from all the tests, and possible gene annotations. In the annotations, genes with high cancer type specific expression are annotated with `(*)` if the above standard entity names were used. There is no header in this file (again, the headers in cursive are just illustrative):

| _chrom_ | _chromStart_ | _fdr_                 | _annotation_   |
|:--------|:-------------|:----------------------|:---------------|
| 1       |    203267500 |  0.034177982009055946 | BTG2(\*), CHIT1 |
| 1       |    203227000 |  0.034177982009055946 | BTG2(\*), CHIT1 |
| 1       |    203272500 |  0.034177982009055946 | BTG2(\*)        |
| 1       |    203275000 |  0.034177982009055946 | BTG2(\*)        |
| 1       |    229577500 | 3.1888888997861688E-4 | ACTA1(\*)       |
| 2       |     89152500 |               1.0E-10 | RPIA(\*)        |
| 2       |     89155000 |               1.0E-10 | RPIA(\*)        |
| 2       |     89157500 |               1.0E-10 | RPIA(\*)        |
| 2       |     89160000 |               1.0E-10 | RPIA(\*)        |
| 2       |     89160000 |               1.0E-10 | CXCR4(\*)       |
| ...     | ...          | ...                   | ...             |

The above file is an excerpt of the results for `Leukemia` in the public ICGC dataset (see below).

# Example analysis

We have prepared a full example analysis which uses the public PCAWG whole genomes from the [ICGC 25k dataset](https://docs.icgc-argo.org/docs/data-access/icgc-25k-data). Since the data is served from an object storage, [AWS Command Line Interface](https://aws.amazon.com/cli/) are required for the download. Please consult the ICGC link for more details and instructions of downloading these files manually.

You can use the following script to download the relevant files (a total of ~900 MB is downloaded, which took <1 min for the author):
```bash
	(cd demo/ && ./demo-10-download.sh)
```

The following script can be used to convert the downloaded files into the appropriate format (taking ~3 min for the author), as detailed above, writing ~150 MB of results under `demo/result_MutationFiles/`:
```bash
	(cd demo/ && ./demo-20-prepare.sh)
```

Finally, the following script runs the model, writing results under `demo/result_Significance/`:
```bash
	(cd demo/ && ./demo-30-significance.sh)
```

This runs the analysis sequentially on each entity (cancer type), so it is not optimal for practical analysis. The analysis has been run on a 10-core computer, taking ~6 h per cancer type, but fewer (or more) CPUs can directly affect the runtime. The script uses 60 GB of heap memory (`-Xmx60G`) for Java, and up to 30 GB of temporary space might be needed.

For running a single entity, here `Leukemia`, the following can be used:
```bash
	(cd demo/ && ./demo-30-significance.sh Leukemia)
```

The following tabulates the full list of significant regions at 0.05 for `Leukemia` in this dataset, with the trivially similar adjacent regions merged:

| chrom | chromStart | chromEnd   | fdr                    | annotation                                    |
|:------|:-----------|:-----------|:-----------------------|:----------------------------------------------|
| 1     |  203267500 |  203277500 |   0.034177982009055946 | BTG2(\*), CHIT1                               |
| 1     |  229567500 |  229570000 |  0.0003188888997861688 | ACTA1(\*)                                     |
| 2     |   89152500 |   89162500 |                  1e-10 | RPIA(\*)                                      |
| 2     |  136875000 |  136877500 |                  1e-10 | CXCR4(\*)                                     |
| 2     |  198262500 |  198272500 | 0.00037691083004286473 | SF3B1(\*)                                     |
| 3     |  186782500 |  186785000 |                  1e-10 | ST6GAL1(\*)                                   |
| 3     |  187462500 |  187465000 |                  1e-10 | BCL6(\*)                                      |
| 3     |  187660000 |  187662500 | 1.0701461680172585e-05 | BCL6(\*), LPP                                 |
| 6     |   91002500 |   91012500 |                  1e-10 | BACH2(\*), MAP3K7(\*), BACH2                  |
| 9     |   37325000 |   37425000 |  2.792588827454938e-05 | ZCCHC7(\*), GRHPR(\*), ZCCHC7                 |
| 9     |  139390000 |  139400000 | 1.5063241821267523e-07 | NOTCH1(\*)                                    |
| 12    |  122425000 |  122525000 |                  1e-10 | WDR66(\*), BCL7A(\*), WDR66, MLXIP(\*), BCL7A |
| 14    |   96177500 |   96187500 |                  1e-10 | TCL1A(\*), C14ORF132                          |
| 14    |  106107500 |  106117500 |                  1e-10 | ADAM6(\*), TMEM121                            |
| 14    |  106170000 |  106180000 |  0.0007221165068917303 | ADAM6(\*), TMEM121                            |
| 14    |  106212500 |  106222500 |                  1e-10 | ADAM6(\*), TMEM121                            |
| 14    |  106240000 |  106242500 |                  1e-10 | ADAM6(\*), TMEM121                            |
| 14    |  106325000 |  106327500 |                  1e-10 | ADAM6(\*), TMEM121                            |
| 14    |  106822500 |  106832500 | 0.00010293548053696062 | ADAM6(\*)                                     |
| 14    |  107257500 |  107260000 |  9.806587301466229e-08 | ADAM6(\*)                                     |
| 22    |   23027500 |   23030000 |                  1e-10 | GGTLC2(\*), RSPH14                            |
| 22    |   23222500 |   23225000 |  6.802780409947805e-05 | GGTLC2(\*), RSPH14                            |
| 22    |   23225000 |   23232500 |                  1e-10 | GGTLC2(\*), RSPH14                            |
| 22    |   23245000 |   23247500 |                  1e-10 | GGTLC2(\*), RSPH14                            |

# Copying

All files are distributed under the [3-clause BSD license](LICENSE.txt). The original source code is Copyright (c) 2021 Felix Dietlein and the added portions are Copyright (c) 2023 Antti Hakkinen. 
