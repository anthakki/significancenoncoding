
# Noncoding variant significance analysis 

This repository contains the updated source code for our noncoding variant discovery analysis. The analysis is based on the work by [Dietlein et al.](https://doi.org/10.1126/science.abg5601) but has been extended to be better suited for large-scale datasets and be more robust to dataset variations, e.g. by more detailed modeling of variant correlations. The original source code was openly provided by [Dietlein et al.](https://doi.org/10.5281/zenodo.5913867) but has been heavily modified since. If you use the modified pipeline, please cite our [work](#) detailing the newer analysis.

# Installation 

Pre-built [Docker images](https://hub.docker.com/r/anthakki/significancenoncoding) are available along with the [Docker build files](https://github.com/anthakki/significancenoncoding-docker) and instructions, which takes care of the dependencies and will download the resources for offline use.

If you prefer to compile the program yourself, e.g. for developement purposes, `make` will build the Java binaries. This requires a [Java JDK 8](https://jdk.java.net/java-se-ri/8-MR6) or newer in `$PATH` and [GNU make](https://www.gnu.org/software/make/). You can simply issue:
```bash
	make
```

# Usage

A [user manual](https://zenodo.org/records/5913867/preview/UserManual.pdf) for the original method is available. This very much still applies for the inputs and outputs of the method as well running it.

In addition, a wrapper shell script is provided, which accepts both analysis and Java options, and will automatically set the Java classpath and locate the support files if they are already downloaded in the same path. For example, following the manual:
```bash
	./SignificanceNoncoding -Xmx55G Breast ~/MutationTestFiles/MutationFiles.txt ~/SignificanceAnalysis/
```
will configure Java with 55 GB of heap, and run the analysis for the entity `Breast` using the input files specified in `~/MutationTestFiles/MutationFiles.txt` and output the results in `~/SignificanceAnalysis/`.

Please refer to the [manual](https://zenodo.org/records/5913867/preview/UserManual.pdf) how to format your input files and interpret the outputs. A set of [test files](http://storage.googleapis.com/noncoding_analysishg19/MutationTestFiles.zip) is also available.

# Copying

All files are distributed under the [3-clause BSD license](LICENSE.txt). The original source code is Copyright (c) 2021 Felix Dietlein and the added portions are Copyright (c) 2023 Antti Hakkinen. 
