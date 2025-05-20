#!/usr/bin/env python3

import gzip
import io
import re

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
		sys.stderr.write(f'Usage: {sys.argv[0]} sample_sheet.tsv entity_map.tsv' + '\n')
		sys.exit(1)

	try:
		opts, (sample_sheet_fn, entity_map_fn) = getopt.getopt(sys.argv[1:], '')
	except getopt.GetoptError:
		usage()
	except ValueError:
		usage()

	entity_by_project = {}

	with gzopen_r(entity_map_fn, 'r') as stream:
		parser = TSVReader(stream)

		for data in parser:
			data = dict(zip(parser.head, data))
			project, entity = data['dcc_project_code'], data['ENTITY']

			assert project not in entity_by_project
			entity_by_project[project] = entity

	out_head = ('SAMPLE', 'DONOR', 'ENTITY', 'FILTER')
	sys.stdout.write('\t'.join(key for key in out_head) + '\n')

	with gzopen_r(sample_sheet_fn, 'r') as stream:
		parser = TSVReader(stream)

		for data in parser:
			data = dict(zip(parser.head, data))
			sample, donor, project = data['aliquot_id'], data['icgc_donor_id'], data['dcc_project_code']

			entity = entity_by_project[project]
			filter = 'PASS'

			exclude = data['donor_wgs_exclusion_white_gray']
			if exclude not in ('Whitelist',):
				entity = '-'
				filter = 'BL' if exclude in ('Blacklist',) else 'GL'   # NB. Blacklist/Graylist

			sample_typ, sample_lib = data['dcc_specimen_type'], data['library_strategy']
			sample_typ_cat, sample_typ_com = re.match('^([\w\s]+)(?: - (.*))?$', sample_typ).groups()

			# NB. skip normals/ cell lines/ RNA-seq etc. not relevant to us

			if sample_typ_cat in ('Primary tumour', 'Metastatic tumour', 'Recurrent tumour'):
				pass
			elif sample_typ_cat in ('Normal',):
				continue
			elif sample_typ_cat in ('Cell line',):
				entity = '-'
				filter = 'CL'   # NB. Cell line
			else:
				assert False

			if sample_lib in ('WGS',):
				pass
			elif sample_lib in ('RNA-Seq',):
				continue
			else:
				assert False

			assert (entity != '-') == (filter == 'PASS')

			out_data = (sample, donor, entity, filter)
			assert len(out_data) == len(out_head)

			sys.stdout.write('\t'.join(str(val) for val in out_data) + '\n')
