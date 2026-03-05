import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import orjson
import ast
import os
import glob
import itertools
from decimal import *
from collections import Counter
from functools import partial
import importlib.util
import sys


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
			create_module(file_path, module_name)
		except (ModuleNotFoundError, ImportError) as e:
			print(e)
			print(f"Module {module_name} not found. Adding to retry stack... for {file_path}")
			stack.append((file_path, module_name))


directory_path = "/groups/vaksler_group/Tal/python/tutils"
create_modules_from_directory(directory_path)

from tutils.databases import DBFactory
from tutils.adapters import GraphLoader

def array_converter(value, sep=';', cast_type = str):
	result = []
	
	if pd.notnull(value) and value != '':
		if sep in value:
			result = [ v.lstrip().strip() for v in str(value).split(sep)[:-1] ]
		else:
			result.append(cast_type(value))
			
	return list(result)

def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_weight', help='output file location', type=str, required=False)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)
	parser.add_argument('--database_type', help='database file location for filtering', type=str, required=False, default='uniprot')


	options, unknown_args = parser.parse_known_args()

	print("protein = {}".format(options.protein))
	print("xgmml_file = {}".format(options.xgmml_file))
	print("output = {}".format(options.output))

	result_dictionary = {}
	result_dictionary['protein'] = options.protein
	result_dictionary['xgmml_file'] = options.xgmml_file
	
	db = None
	if options.database:
		print("database = {}".format(options.database))
		db = DBFactory.create(options.database_type, path = options.database)
		result_dictionary['database'] = options.database
		result_dictionary['database_type'] = options.database_type
		
	graph = GraphLoader.load_adapter('igraph', options.xgmml_file, db)['adapter'].current_state
	
	
	result_dictionary['n_vertices'] = graph.vcount()
	result_dictionary['n_edges'] = graph.ecount()
	
	degrees = graph.vs.degree()
	result_dictionary['min_degree'] = min(degrees)
	result_dictionary['max_degree'] = max(degrees)
	
	if options.alignment_weight:
		print("alignment_weight = {}".format(options.alignment_weight))
		result_dictionary['alignment_weight'] = options.alignment_weight
		result_dictionary['min_alignment_weight'] = min(graph.es[options.alignment_weight])
		result_dictionary['max_alignment_weight'] = max(graph.es[options.alignment_weight])

	
	json_bytes = orjson.dumps(result_dictionary, option=orjson.OPT_SERIALIZE_NUMPY | orjson.OPT_INDENT_2)

	with open(options.output, "wb") as f:
		f.write(json_bytes)

if __name__ == '__main__':
	_main()


