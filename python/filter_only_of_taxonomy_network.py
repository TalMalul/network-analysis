import argparse
import ast
import os
import itertools
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

os.environ['PYTHONPATH'] = '/groups/vaksler_group/Tal/python:' + os.environ.get('PYTHONPATH', '')
directory_path = "/groups/vaksler_group/Tal/python/tutils"  # Replace with the path to your directory
create_modules_from_directory(directory_path)


from xgmmlgraphadapter import XGMMLReader, XGMMLWriter
from tutils.databases import UniprotDB

def only_in_or_filter_network_by_attributes(graph, db, filter_dict):

	if db is not None:
		vs_filtered = set()
		for attribute, values in filter_dict.items():
			for value in values:
				vs_filtered.update(graph.vs.select(name_in = db.loc[db[attribute] == value].index))
		graph = graph.subgraph(vs_filtered)
		return graph

def only_in_and_filter_network_by_attributes(graph, filter_dict):
	for attribute, values in filter_dict.items():
		vs_filtered = set()
		for value in values:
			vs_filtered.update(graph.vs.select(lambda v: all(a == value for a in v[attribute])))
		graph = graph.subgraph(vs_filtered)
	return graph
	
def array_converter(value):
	result = []

	if pd.notnull(value) and value != '':
		if ';' in value:
			result = str(value).strip().split(';')[:-1]
		else:
			result = ast.literal_eval(value)

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

def _main():

	
	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--protein', help='the name of the protien', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--output', help='output file location', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--taxonomy_level', help='the taxonomy level', type=str, required=True)
	parser.add_argument('--taxonomy', nargs='+', help='the taxonomy', type=str, required=True)


	options = parser.parse_args()
		

	print("protein = " + options.protein)
	print("output = " + options.output)
	print("xgmml_file = {}".format(options.xgmml_file))
	print("taxonomy_level = {}".format(options.taxonomy_level))
	print("taxonomy = {}".format(options.taxonomy))
	
	uniprot_db = None
	if options.database:
		print("database = {}".format(options.database))
		uniprot_db = load_db(options.database)

	graph = load_network(options.xgmml_file, uniprot_db)

	graph = only_in_or_filter_network_by_attributes(graph, uniprot_db, {options.taxonomy_level: options.taxonomy})
	nx_graph = graph.to_networkx()
	

	f = open(options.output, 'w')
	networkxgmml.XGMMLWriter(f ,nx_graph, '',directed=False)
	f.close()

if __name__ == '__main__':
	_main()


