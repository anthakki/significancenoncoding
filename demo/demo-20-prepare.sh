#!/bin/sh

set -e

./pcawg-entities.py 'PCAWG/donors_and_biospecimens/pcawg_sample_sheet.tsv' 'pcawg_entity_map.tsv' >'entities.tsv'

mkdir -p 'result_MutationFiles/'

./pcawg-maf2snmaf.py 'PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz' |
	./snmaf-scatter.py -z -l '' -o 'result_MutationFiles/{}.maf.gz' '/dev/stdin' 'entities.tsv'
