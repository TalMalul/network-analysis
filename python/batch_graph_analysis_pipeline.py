import argparse
import json
import os
import glob
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
			print(file_path, module_name, flush=True)
			create_module(file_path, module_name)
		except (ModuleNotFoundError, ImportError) as e:
			print(e, flush=True)
			print(f"Module {module_name} not found. Adding to retry stack... for {file_path}")
			stack.append((file_path, module_name))

os.environ['PYTHONPATH'] = '/groups/vaksler_group/Tal/python:' + os.environ.get('PYTHONPATH', '')
directory_path = "/groups/vaksler_group/Tal/python/tutils"  # Replace with the path to your directory
create_modules_from_directory(directory_path)


from tutils.pipelines import GraphAnalysisPipelineTemplate
from tutils.services import GraphRegistry
from tutils.adapters import IGraphAdapter, NetworkXAdapter, GraphToolAdapter

GraphRegistry.register_package(
	package_name='igraph', 
	method_name='visit_igraph', 
	adapter_class=IGraphAdapter
)

GraphRegistry.register_package(
	package_name='networkx', 
	method_name='visit_networkx', 
	adapter_class=NetworkXAdapter
)

GraphRegistry.register_package(
	package_name='graphtool', 
	method_name='visit_graphtool', 
	adapter_class=GraphToolAdapter
)

class SimplePipeline(GraphAnalysisPipelineTemplate):
	def __init__(self, xgmml_file, database_file, database_type, settings_file, output_folder, graph_filter_iterration, weights):
		super().__init__(xgmml_file, database_file, database_type, settings_file, output_folder, graph_filter_iterration, weights)
		
def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--output_folder', help='output folder location', type=str, required=True)
	parser.add_argument('--settings_file', help='settings_file', type=str, required=True)
	parser.add_argument('--xgmml_file', help='output file location', type=str, required=True)
	parser.add_argument('--database', help='database file location for filtering', type=str, required=True)
	parser.add_argument('--database_type', help='database file location for filtering', type=str, required=False, default='uniprot')
	parser.add_argument('--itteration', help='database file location for filtering', type=int, required=False, default=0)
	parser.add_argument('--alignment_weights', help='database file location for filtering', type=str, required=False, default=None)
	

	options, unknown_args = parser.parse_known_args()
	
	
		

	print("output_folder = {}".format(options.output_folder), flush=True)
	print("settings_file = {}".format(options.settings_file), flush=True)

	kwargs = {}
	while len(unknown_args) > 0:
		k = args.pop(0)
		v = args.pop(0)
		kwargs[k] = v
	
	weights = None
	
	if options.alignment_weights is not None:
		with open(options.alignment_weights) as json_file:
			weights = json.load(json_file)
	
	print('start pipeline', flush=True)
	sys.stdout.flush()
	simple_pipeline = SimplePipeline(options.xgmml_file, options.database, options.database_type, options.settings_file, output_folder = options.output_folder, graph_filter_iterration = options.itteration, weights = weights)	
	simple_pipeline.run()
	
if __name__ == '__main__':
	_main()


