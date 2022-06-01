#!/usr/bin/env python3

import sys
import pyroe
from pyroe import load_fry

def get_mode_num_genes(frydir, quiet=True):
	import json
	import sys
	import os
	# since alevin-fry 0.4.1 the generic "meta_info.json"
	# has been replaced by a more informative name for each
	# sub-command. For quantification, it is "quant.json".
	# we check for both files here, in order.
	meta_info_files = ["quant.json", "meta_info.json"]

	fpath = os.path.sep.join([frydir, meta_info_files[0]])
	# first, check for the new file, if we don't find it, check
	# for the old one.
	if not os.path.exists(fpath):
		if not quiet:
			print(f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}.", file=sys.stderr)
		fpath = os.path.sep.join([frydir, meta_info_files[1]])
		# if we don't find the old one either, then return None
		if not os.path.exists(fpath):
			raise IOError(f"Found no {meta_info_files[1]} file either; cannot proceed.")

	# if we got here then we had a valid json file, so 
	# use it to get the number of genes, and if we are 
	# in USA mode or not.
	meta_info = json.load(open(fpath))
	ng = meta_info['num_genes']
	usa_mode = meta_info['usa_mode']
	if not quiet:
		print(f"USA mode: {usa_mode}", file=sys.stderr)
	return usa_mode, ng

def compare_quants(args):
	import json

	ref_quant_dir = args.ref_quant	
	test_quant_dir = args.test_quant

	mode_ng_ref = get_mode_num_genes(ref_quant_dir)
	mode_ng_test = get_mode_num_genes(test_quant_dir)
	if  mode_ng_ref != mode_ng_test:
		print(f"Cannot compare a quantification result of type {mode_ng_ref} to one of type {mode_ng_test}", file=sys.stderr)
		sys.exit(1)

	usa_mode = mode_ng_ref[0]
	odict = None
	if usa_mode:
		a = load_fry(ref_quant_dir, output_format='raw')
		b = load_fry(test_quant_dir, output_format='raw')
		
		odict = { "nobs_ref" : a.n_obs, "nobs_test" : b.n_obs }

		odict['diff_U'] = float(abs(a.layers['unspliced'] - b[ a.obs_names, : ].layers['unspliced']).sum())
		odict['diff_S'] = float(abs(a.layers['spliced'] - b[ a.obs_names, : ].layers['spliced']).sum())
		odict['diff_A'] = float(abs(a.layers['ambiguous'] - b[ a.obs_names, : ].layers['ambiguous']).sum())

		odict['obs_ref-obs_test'] = list(a.obs_names.difference(b.obs_names))
		odict['obs_test-obs_ref'] = list(b.obs_names.difference(a.obs_names))
	else:
		a = load_fry(ref_quant_dir, output_format='raw')
		b = load_fry(test_quant_dir, output_format='raw')
		
		odict = { "nobs_ref" : a.n_obs, "nobs_test" : b.n_obs }

		odict['diff_X'] = float(abs(a.X - b[ a.obs_names, : ].X).sum())

		odict['obs_ref-obs_test'] = list(a.obs_names.difference(b.obs_names))
		odict['obs_test-obs_ref'] = list(b.obs_names.difference(a.obs_names))

	with open(args.output, 'w') as ofile:
		json.dump(odict, ofile, sort_keys=True, indent=4)


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Compare two quantification results')
	parser.add_argument('ref_quant', type=str, 
                help='The reference quantification result')
	parser.add_argument('test_quant', type=str, 
                help='The test quantification result to compare against the reference')
	parser.add_argument('output', type=str,
		help='Where to write the output report')
	args = parser.parse_args()
	compare_quants(args)
