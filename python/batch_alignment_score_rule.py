import argparse
import networkx as nx
import networkxgmml
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import json
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

from tutils.databases import UniprotDB

def array_converter(value, sep=';', cast_type = str):
	result = []
	
	if pd.notnull(value) and value != '':
		if sep in value:
			result = [ v.lstrip().strip() for v in str(value).split(sep)[:-1] ]
		else:
			result.append(cast_type(value))
			
	return list(result)

def load_db(location):

	uniprot_db = UniprotDB.read(location)

	return uniprot_db


def load_network(location, db = None):
	f = open(location, 'rb')
	nx_g = networkxgmml.XGMMLReader(f)
	f.close()
	graph = ig.Graph.from_networkx(nx_g)
	nx_g.clear()
	
	graph.to_undirected(mode='each')
	graph.vs['name'] = graph.vs['_nx_name']
	
	
	if db is not None:
		db_missing_entries = set(graph.vs['name']) - set(db.index)
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in db_missing_entries))
		graph_fragments_entries = set(graph.vs['name']).intersection(set(db[(db['Fragment'])].index))
		graph.delete_vertices(graph.vs.select(lambda v: v['name'] in graph_fragments_entries))

	return graph

def get_edges_le(graph, degree):
	return graph.vs.select(_degree_le = degree)

def get_edges_ge(graph, degree):
	return graph.vs.select(_degree_ge = degree)
	
def edges_function(rule):
	if rule == 'GE':
		return get_edges_ge
	elif rule == 'LE':
		return get_edges_le
		
		
def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--alignment_scores', help='degree file location', type=int, nargs='+', required=True)
	parser.add_argument('--alignment_weight', help='degree file location', type=str, required=True)
	parser.add_argument('--relevances', help='the taxonomy relevancy group', type=str, nargs='*', required=True, default=[])
	parser.add_argument('--delete_vertices', help='the taxonomy relevancy group', type=bool, required=False, default=False, action=argparse.BooleanOptionalAction)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=False)
	parser.add_argument('--rule', help='name of edge attribute', type=str, required=False, default='GE')
	parser.add_argument('--analyses', help='', type=str, nargs='*', required=False, default=[])

	options, unknown_args = parser.parse_known_args()

	print("protein = {}".format(options.protein))
	print("graph = {}".format(options.xgmml_file))
	print("alignment_scores = {}".format(options.alignment_scores))
	print("rule = {}".format(options.rule))
	print("alignment_weight = {}".format(options.alignment_weight))
	print("relevances = {}".format(options.relevances))
	print("analyses = {}".format(options.analyses))
	
	kwargs = {}
	while len(unknown_args) > 0:
		k = args.pop(0)
		v = args.pop(0)
		kwargs[k] = v

	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)
	
	result_dictionary = {}

	for alignment_score in options.alignment_scores:
		result_dictionary[alignment_score] = {}
		result_dictionary[alignment_score]['function_clustering'] = {}
		result_dictionary[alignment_score]['taxonomy_clustering'] = {}
		
		graph_filter = FilterFactory.create_filter('alignment', alignment_weight = options.alignment_weight, alignment_score_rule = options.rule, alignment_score = alignment_score, delete_vertices = options.delete_vertices)
		
		taxonomy_clustering_result = get_graph_taxonomy_modularity(selected_analyses = ['gini', 'custom', 'chi', 'rand_index'], graph = graph, alignment_weight = options.alignment_weight, graph_filter = graph_filter, db = uniprot_db)

		for taxonomy_level, analysis_result in taxonomy_clustering_result.items():
			result_dictionary[alignment_score]['taxonomy_clustering'][taxonomy_level] = analysis_result.to_dict()
		
		for analysis in options.analyses:
			analysis = analysis.lower()
			algoritm_name, str_analysis_names = analysis.split(':') if ':' in analysis else (analysis, '')
			algorithm_function = ClusteringAnalysisDictionaryFactory.create_clustering_algorithm(algoritm_name)
			
			list_clustering_dictionary_decorator = []

		
		
			analysis_names = [a for a in str_analysis_names.split(',') if a != '']
			if analysis_names is not None:
				for analysis_name in analysis_names:
					clustering_dictionary_decorator = ClusteringAnalysisDictionaryFactory.create_analysis(analysis_name)
					list_clustering_dictionary_decorator.append(clustering_dictionary_decorator)
			
			algoritm_analysis = algorithm_function(graph, alignment_weight = options.alignment_weight, graph_filter = graph_filter, **kwargs)			
			
			
			
			result_dictionary[alignment_score]['function_clustering'][algoritm_name] = algoritm_analysis.to_dict()
			for function_relevancy_str in options.relevances:
				relevance_name, relevance_config = function_relevancy_str.split(':')
				relevance_xlsx, index_column_name, relevancy_colum_name, relevancy_type, reviewed = relevance_config.split(';')
				relevancy = Relevancy(relevance_xlsx, index_column_name, relevancy_colum_name, RelevancyFactory.get_method(relevancy_type))

				relevancy_analysis = AnalysisDictionary({})
				for selected_analysis in list_clustering_dictionary_decorator:
					relevancy_analysis = selected_analysis(relevancy_analysis, algoritm_analysis.vertex_clustering.graph, algoritm_analysis.vertex_clustering, relevancy, reviewed.lower() == 'reviewed', uniprot_db)
					
				result_dictionary[alignment_score]['function_clustering'][algoritm_name][relevance_name] = relevancy_analysis.to_dict()
				result_dictionary[alignment_score]['function_clustering'][algoritm_name][relevance_name]['relevance_xlsx'] = relevance_xlsx
				result_dictionary[alignment_score]['function_clustering'][algoritm_name][relevance_name]['reviewed'] = reviewed.lower() == 'reviewed'
			
			
				
		simple_result = simple_analysis(graph = graph, graph_filter = graph_filter)
	
		
		result_dictionary[alignment_score]['database'] = options.database
		result_dictionary[alignment_score]['xgmml_file'] = options.xgmml_file
		result_dictionary[alignment_score]['graph_analysis'] = simple_result.to_dict()		
	
	
	json_object = json.dumps(result_dictionary, indent = 4)

	with open(options.output, "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


