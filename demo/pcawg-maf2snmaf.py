#!/usr/bin/env python3

import gzip
import io

def gzopen_r(filename, mode = 'r', magic = b'\x1f\x8b'):
	file = open(filename, mode = 'rb')

	if file.peek(len(magic))[:len(magic)] == magic:
		file = gzip.GzipFile(fileobj = file, mode = mode)
	if 'b' not in mode:
		file = io.TextIOWrapper(file)

	return file

def chomp(line):
	if line.endswith('\n'):
		line = line[:-1]
	if line.endswith('\r'):
		line = line[:-1]

	return line

class TSVReader:
	def __init__(self, stream, header = True):
		self._stream = stream
		self.head = next(iter(self)) if header else None 

	def __iter__(self):
		for line in (chomp(line) for line in self._stream):
			yield line.split('\t')

if __name__ == '__main__':
	import getopt
	import sys

	def usage():
		sys.stderr.write(f'Usage: {sys.argv[0]} input.maf [...]' '\n')
		sys.exit(1)

	try:
		opts, args = getopt.getopt(sys.argv[1:], '')
	except getopt.GetoptError:
		usage()

	for key, val in opts:
		pass

	if len(args) < 1:
		usage()

	chrom_values = []
	for idx in range(22):
		chrom_values.append(f'{idx + 1}')
	for chrom in ('X', 'Y'):
		chrom_values.append(chrom)

	out_head = ('Chromosome', 'Position', 'Reference_Allele', 'Tumor_Seq_Allele', 'Tumor_Sample_Barcode')
	sys.stdout.write('\t'.join(key for key in out_head) + '\n')

	for arg in args:
		with gzopen_r(arg, 'r') as stream:
			parser = TSVReader(stream)

			(*ixs, ) = (parser.head.index(key) for key in ('Chromosome', 'Strand',
				'Start_position', 'End_position',
				'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
				'ref_context',
				'Tumor_Sample_Barcode'))

			for data in parser:
				chrom, strand, pos, end_pos, ref, germ_ref, alt, ref_ctx, sample = (data[ix] for ix in ixs)

				# NB. PCAWG specific assumptions
				assert chrom in chrom_values
				assert strand == '+'
				assert germ_ref == ref

				def parse_ale(seq):
					if seq == '-':
						assert seq == '-'
						return ''
					else:
						assert len(seq) > 0 and all(base in 'ACGT' or base in 'N' for base in seq)
						return seq

				# NB. allele is either "-" for empty or non-empty sequence 
				ref, alt = (parse_ale(seq) for seq in (ref, alt))
				pos, end_pos = (int(val) for val in (pos, end_pos))
				pad = 10

				# NB. PCAWG uses 1-based coordinates
				 # for "ins" range includes the flanking nucleotides
				 # and pos does point to the left flanking base in the "ins" case
				assert end_pos - pos == ((len(ref) - 1) if ref != '' else 1)

				if ref == '':
					assert end_pos == pos + 1   # NB. fix so they are 1-based closed ranges
					pos, end_pos = end_pos, pos
					pad = pad + 1   # NB. extra padding is included if flanks included 

				assert end_pos - pos == len(ref) - 1   # sanity check..

				ref_ctx = ref_ctx.upper()
				assert all(base in 'ACGT' or base in 'N' for base in ref_ctx)

				assert len(ref_ctx) == pad + len(ref) + pad
				assert ref_ctx[pad:(pad + len(ref))] == ref

				# NB. add VCF-style left flank if necessary
				if ref == '' or alt == '':
					assert pad > 0
					lbase = ref_ctx[pad - 1]

					ref, alt = lbase + ref, lbase + alt
					pos, end_pos = pos - 1, end_pos

				assert len(ref) > 0   # sanity check..
				assert len(alt) > 0
				assert end_pos - pos == len(ref) - 1

				# NB. there are a few "insN" mutations, and left flank "N"s
				 # skip those
				if 'N' in ref or 'N' in alt:
					sys.stderr.write(f'warning: \'N\' in sequence, dropping variant..' '\n')
					continue

				out_data = (chrom, pos, ref, alt, sample)
				assert len(out_data) == len(out_head)

				sys.stdout.write('\t'.join(str(val) for val in out_data) + '\n')
