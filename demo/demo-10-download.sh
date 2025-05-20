#!/bin/sh

set -e

# see https://docs.icgc-argo.org/docs/data-access/icgc-25k-data

s3_bucket='s3://icgc25k-open'
s3_endpoint='https://object.genomeinformatics.org'

filenames="$(cat <<-EOT
 	PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz
	PCAWG/donors_and_biospecimens/pcawg_sample_sheet.tsv
EOT)"

shell() {
	echo '+' "$@" >&2
	"$@"
}

echo "$filenames" | while read filename; do
	if ! [ -e "$filename" ]; then
		shell aws s3 cp --no-sign-request --endpoint-url "$s3_endpoint" "$s3_bucket/$filename" "$filename"
	fi
done
