#!/bin/sh

set -e

shell() {
	echo '+' "$@"
	"$@"
}

run_Significance() {
	entity="$1"
	shell ../SignificanceNoncoding -Xmx60G -seed=123 -z "$entity" 'result_MutationFiles/MutationFiles.txt' 'result_Significance/'
}

mkdir -p 'result_Significance/'

if [ "$#" -gt 1 ]; then
	for entity in "$@"; do
		run_Significance "$entity"
	done
else
	for filename in 'result_MutationFiles/'*'.maf.gz'; do
		entity="$(basename "$filename")"
		entity="${entity%%.maf.gz}"

		run_Significance "$entity"
	done
fi
