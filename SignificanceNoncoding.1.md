% SignificanceNoncoding(1) Version 0.9 | General Commands

# NAME

**SignificanceNoncoding** - Noncoding variant significance analysis

# SYNOPSIS

| **SignificanceNoncoding** \[_java_options_\] \[-k\] \[-z\] \[_options_\] _entity_ _path_file_ _output_folder_ \[_annotation_folder_\]

# DESCRIPTION

## Common Java options

-cp, -classpath _path_
: Add _path_ to Java class path

-Xmx*size*
: Set maximum heap size for Java

## Options

-default_annotation_folder _annotation_folder_
: Specifies the path to the annotation files as _annotation_folder_ (default: automatic)

-disable_test _test_
: Disables the specified test _test_

-no_mutationfiles, -always_mutationfiles
: Skip/ always read and filter the input mutations

-no_readall, -always_readall
: Skip/ always compute the mutation statistics for all cancer types to be used in the reference models

-no_significance, -always_significance
: Skip/ always compute significance models

-no_combine, -always_combine
: Skip/ always compute model pooling

-k
: Keep intermediate files (default: delete)

-jar_path _path_
: Specifies the path to the main Jar (default: automatic)

-seed=_seed_
: Seed the PRNG

-z, -z*level*
: Compress output files using gzip. Optionally, compression level _level_ (`1` to `9`) can be specified.
