import argparse
import json
import ast
import os
import importlib.util
import sys
import itertools
import igraph as ig
import numpy as np


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


from tutils.databases import UniprotDB
from tutils.graphs import IGraphAdapter

def cantor_pair(x: int, y: int) -> int:
		"""Cantor pairing function: uniquely maps a pair of integers to one integer."""
		return (x + y) * (x + y + 1) // 2 + y

def precompute_edges_by_threshold_numpy_index(g: ig.Graph, weight_attr: str, relevant_weights = None, undirected=True):
	"""
	Precompute edges above each threshold using vertex indices and Cantor pairing.
	Returns dict: threshold -> set of edge IDs (integers).
	"""
	edges_by_threshold = {}

	# Edge endpoints and weights
	sources = np.array([e.source for e in g.es])
	targets = np.array([e.target for e in g.es])
	weights = np.array(g.es[weight_attr])

	# Create edge IDs using Cantor pairing
	if undirected:
		# sort endpoints for undirected edges
		sorted_pairs = np.sort(np.column_stack((sources, targets)), axis=1)
		edge_ids = np.array([cantor_pair(u, v) for u, v in sorted_pairs])
	else:
		edge_ids = np.array([cantor_pair(u, v) for u, v in zip(sources, targets)])

	# Counting sort by weight
	max_weight = weights.max()
	min_weight = weights.min()
	weight_buckets = [[] for _ in range(max_weight + 1)]

	
	for eid, w in zip(edge_ids, weights):
		weight_buckets[w].append(eid)

	current_edges = set()
	min_relevant_weights = min(relevant_weights)
	
	for w in range(max_weight, min_weight - 1, -1):
		if w < min_relevant_weights:
			break
			
		current_edges.update(weight_buckets[w])
		if w in relevant_weights:
			edges_by_threshold[w] = current_edges.copy()

	return edges_by_threshold
	

def compute_jaccard_over_weight(g1, g2, weight_attr, g1_weights_window, g2_weights_window):
	"""
	Compute Jaccard over all thresholds using sliding windows of relevant_weights in parallel.
	"""
	

	jaccard_results = {}

	
	edges_g1_dict = precompute_edges_by_threshold_numpy_index(g1, weight_attr, relevant_weights=g1_weights_window)
	edges_g2_dict = precompute_edges_by_threshold_numpy_index(g2, weight_attr, relevant_weights=g2_weights_window)

				

	threshold_pairs = [(t1, t2) for t1, t2 in itertools.product(g1_weights_window, g2_weights_window)]

	for t1, t2 in threshold_pairs:
		intersection = edges_g1_dict[t1] & edges_g2_dict[t2]
		union = edges_g1_dict[t1] | edges_g2_dict[t2]
		jaccard = len(intersection) / len(union) if len(union) > 0 else -1.0
		jaccard_results[(t1, t2)] = jaccard
			
	return jaccard_results
		
		



def _main():

	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--xgmml_file_g1', help='output file location', type=str, required=True)
	parser.add_argument('--xgmml_file_g2', help='output file location', type=str, required=True)
	parser.add_argument('--output_folder', help='output file location', type=str, required=True)
	parser.add_argument('--start_g1', help='starting node index', type=int, required=True)
	parser.add_argument('--start_g2', help='starting node index', type=int, required=True)
	parser.add_argument('--step_g1', help='step from strating node', type=int, required=True)
	parser.add_argument('--step_g2', help='step from strating node', type=int, required=True)
	parser.add_argument('--alignment_weight', help='the name for the alignment weight "_score" will be added at the end', type=str, default='alignment')


	options = parser.parse_args()

	database = None
	if options.database:
		 database = UniprotDB(options.database)

		
	ga1 = IGraphAdapter(options.xgmml_file_g1, database)
	ga1.load()

	ga2 = IGraphAdapter(options.xgmml_file_g2, database)
	ga2.load()
	
	g1 = ga1.graph
	g2 = ga2.graph


	g1_min, g1_max = min(g1.es[options.alignment_weight]), max(g1.es[options.alignment_weight])
	g2_min, g2_max = min(g2.es[options.alignment_weight]), max(g2.es[options.alignment_weight])
   
	starting_index_g1 = options.start_g1 if options.start_g1 > g1_min else g1_min
	starting_index_g2 = options.start_g2 if options.start_g2 > g2_min else g2_min
	
	stoping_index_g1 = starting_index_g1 + options.step_g1 if starting_index_g1 + options.step_g1 < g1_max + 1 else g1_max + 1
	stoping_index_g2 = starting_index_g2 + options.step_g2 if starting_index_g2 + options.step_g2 < g2_max + 1 else g2_max + 1
	
	g1_weights_window = list(range(starting_index_g1, stoping_index_g1))
	g2_weights_window = list(range(starting_index_g2, stoping_index_g2))
	
	jaccard_results = compute_jaccard_over_weight(g1, g2, options.alignment_weight, g1_weights_window, g2_weights_window)
	jaccard_results_list = [{"{}_g1".format(options.alignment_weight): t1, "{}_g2".format(options.alignment_weight): t2, "jaccard": val} for (t1, t2), val in jaccard_results.items()]
	
	data = {}
	data['g1'] = options.xgmml_file_g1
	data['g2'] = options.xgmml_file_g2
	data['database'] = options.database
	data['g1_weights_window'] = g1_weights_window
	data['g2_weights_window'] = g2_weights_window
	data['alignment_weight'] = options.alignment_weight
	data['jaccard_results'] = jaccard_results_list
	
	json_object = json.dumps(data, indent = 4)
	
	with open('{}/g1-{}-{}_g2-{}-{}.json'.format(options.output_folder, starting_index_g1, stoping_index_g1, starting_index_g2, stoping_index_g2), "w") as outfile:
		outfile.write(json_object)

if __name__ == '__main__':
	_main()


