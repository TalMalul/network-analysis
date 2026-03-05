import argparse
import json
import os
import glob
import importlib.util
import sys
from subprocess import run, PIPE


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


directory_path = "/groups/vaksler_group/Tal/python/tutils"  # Replace with the path to your directory
create_modules_from_directory(directory_path)

from tutils.settings import SettingParser

		
def _main():


	parser = argparse.ArgumentParser(description="Find motifs of a graph")
	parser.add_argument('--output_folder', help='output folder location', type=str, required=True)
	parser.add_argument('--settings_file', help='settings_file', type=str, required=True)
	parser.add_argument('--network_general', help='networkx general file', type=str, required=True)
	parser.add_argument('--network_config', help='network config file', type=str, required=True)
	

	options, unknown_args = parser.parse_known_args()

	print("output_folder = {}".format(options.output_folder))
	print("settings_file = {}".format(options.settings_file))
	print("network_general = {}".format(options.network_general))
	print("network_config = {}".format(options.network_config))
	
	
	kwargs = {}
	while len(unknown_args) > 0:
		k = args.pop(0)
		v = args.pop(0)
		kwargs[k] = v

	
	runner_batch_settings = '/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/bash/runner_batch_settings.sh'
	with open(options.settings_file) as json_file:
		json_config = json.load(json_file)
		parser = SettingParser(json_config)

		filters = parser.parse_filters()

		command = ['bash', runner_batch_settings, '--network_general', options.network_general, '--network_config', options.network_config, '--output_folder', options.output_folder]
		result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
		print(result.stdout, result.stderr)

	

	
if __name__ == '__main__':
	_main()


