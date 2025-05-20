#!/usr/bin/env python3

import gzip
import io
import os.path
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
		sys.stderr.write(f'Usage: {sys.argv[0]} [-z] [-l list.txt] -o output.maf input.maf entities.tsv' + '\n')
		sys.exit(1)

	write_gzip = False
	out_fn = ''
	list_fn = ''

	try:
		opts, (maf_fn, entities_fn) = getopt.getopt(sys.argv[1:], 'z l::o:')
	except getopt.GetoptError:
		usage()
	except ValueError:
		usage()

	for key, val in opts:
		if key == '-z':
			write_gzip = True
		if key == '-o':
			out_fn = val
		if key == '-l':
			list_fn = True

	if not out_fn:
		usage()

	entities_by_sample = {}
	entities = set()

	with gzopen_r(entities_fn, 'r') as stream:
		parser = TSVReader(stream)

		for data in parser:
			data = dict(zip(parser.head, data))

			assert (data['ENTITY'] != '-') == (data['FILTER'] == 'PASS')

			if data['ENTITY'] != '-':
				assert data['SAMPLE'] not in entities_by_sample
				entities_by_sample[data['SAMPLE']] = data

				entities.add(data['ENTITY'])

	streams = {}

	def open_stream(filename, head):
		assert entity not in streams

		if write_gzip:
			streams[entity] = io.TextIOWrapper(gzip.GzipFile(filename, 'w', compresslevel = 9, mtime = 0))
		else:
			streams[entity] = open(filename, 'w')

		streams[entity].write('\t'.join(head) + '\n')

	with gzopen_r(maf_fn, 'r') as stream:
		parser = TSVReader(stream)

		sample_ix = parser.head.index('Tumor_Sample_Barcode')

		for data in parser:
			sample = data[sample_ix]

			if sample not in entities_by_sample:
				continue

			entity = entities_by_sample[sample]['ENTITY']
			donor = entities_by_sample[sample]['DONOR']

			out_data = list(data)
			out_data[sample_ix] = donor

			if entity not in streams:
				open_stream(out_fn.replace('{}', entity), parser.head)

			streams[entity].write('\t'.join(out_data) + '\n')

		for entity in entities:
			if entity not in streams:
				open_stream(out_fn.replace('{}', entity), parser.head)

	for entity in streams:
		streams[entity].close()

	if list_fn is True:
		list_fn = os.path.join(re.sub(r'[^/]*[{][}].*$', '', out_fn), 'MutationFiles.txt')

	if list_fn:
		with open(list_fn, 'w') as stream:
			for entity in sorted(entities):
				stream.write('\t'.join(str(val) for val in (entity, os.path.relpath(out_fn.replace('{}', entity), os.path.dirname(list_fn)))) + '\n')
