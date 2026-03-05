import argparse
import json
import ast
import os
import importlib.util
import sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import vaex
import pyarrow
from Bio.SeqRecord import SeqRecord
from decimal import *


def create_module(file_path, module_name):
	spec = importlib.util.spec_from_file_location(module_name, file_path)
	new_module = importlib.util.module_from_spec(spec)
	sys.modules[module_name] = new_module
	spec.loader.exec_module(new_module)

def create_modules_from_directory(directory, root_name="tutils"):
	stack = []

	for foldername, subfolders, filenames in os.walk(directory):
		for filename in filenames:
			if filename.endswith(".py"):
				file_path = os.path.join(foldername, filename)
				rel_path = os.path.relpath(file_path, directory)
				module_parts = rel_path.split(os.path.sep)
				module_parts[-1] = os.path.splitext(module_parts[-1])[0]
				if module_parts[-1] == "__init__":
					module_parts = module_parts[:-1]
				if module_parts:
					module_name = ".".join([root_name] + module_parts)
				else:
					module_name = root_name  # Root __init__.py

				stack.append((file_path, module_name))

	while stack:
		file_path, module_name = stack.pop(0)
		try:
			print(file_path, module_name, flush=True)
			create_module(file_path, module_name)
		except (ModuleNotFoundError, ImportError) as e:
			print(e, flush=True)
			print(f"Module {module_name} not found. Adding to retry stack... for {file_path}")
			stack.append((file_path, module_name))

os.environ['PYTHONPATH'] = '/groups/vaksler_group/Tal/python:' + os.environ.get('PYTHONPATH', '')
directory_path = "/groups/vaksler_group/Tal/python/tutils"  # Replace with the path to your directory
create_modules_from_directory(directory_path)

from tutils.databases import *
from tutils.sequence_alignment import get_local_score


def _main():

	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--output_folder', help='output file location', type=str, required=True)
	parser.add_argument('--start', help='starting node index', type=int, required=True)
	parser.add_argument('--step', help='step from strating node', type=int, required=True)
	parser.add_argument('--alignment_name', help='the name for the alignment weight "_score" will be added at the end', type=str, default='alignment')
	parser.add_argument('--gap_open_penalty', help='the cost of opening a gap', type=int, default=0)
	parser.add_argument('--gap_extend_penalty', help='the cost of extending an open gap', type=int, default=-1)

	options = parser.parse_args()

	database = None
	if options.database:
		 database = UniprotDB(options.database)

		
	stoping_index = options.start + options.step if options.start + options.step < len(database.dataframe.index) else len(database.dataframe.index)
	hdf5s = []
	for idx_source in range(options.start, stoping_index, 1):
		alignment_scores = []
		source_name = database.dataframe.index[idx_source]
		output_path = '{}/{}_{}.hdf5'.format(options.output_folder, source_name, idx_source)
		
		if idx_source + 1 < len(database.dataframe.index):
			hdf5s.append(output_path)
		if os.path.exists(output_path):
			continue

		for idx_target in range(idx_source + 1, len(database.dataframe.index)):
			
			target_name = database.dataframe.index[idx_target]
			
			source_sequence = database.dataframe.loc[source_name,'Sequence']
			target_sequence = database.dataframe.loc[target_name,'Sequence']
			source_record = SeqRecord(source_sequence, id=source_name)
			target_record = SeqRecord(target_sequence, id=target_name)
			
			alignment = get_local_score(source_record, target_record)
			alignment_identity = alignment['%local_id']
			normalized_alignment_score = alignment['normalized_local_score']
			alignment_score = alignment['local_score']
			alignment_scores.append({'source': source_name, 'target': target_name, '{}_score'.format(options.alignment_name): alignment_score, '%{}_id'.format(options.alignment_name): alignment_identity , 'normalized_{}_score'.format(options.alignment_name): normalized_alignment_score})


		df = pd.DataFrame(alignment_scores)
		vdf = vaex.from_pandas(df, name='{}_{}_{}'.format(options.protein, source_name, idx_source), copy_index=False)
		vdf.export_hdf5(output_path)
		vdf.close()
		
	vdf = vaex.open_many(hdf5s)
	vdf.export_hdf5('{}/{}_{}.hdf5'.format(options.output_folder, options.start, stoping_index))
	vdf.close()
	#json_object = json.dumps(alignment_scores, indent = 4)

	#with open(options.output + '_' + str(options.start) + '_' + str(stoping_index), "w") as outfile:
		#outfile.write(json_object)

if __name__ == '__main__':
	_main()


