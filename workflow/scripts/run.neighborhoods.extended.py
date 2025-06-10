import argparse
import os

import pandas as pd
from neighborhoods import read_bed, count_features_for_bed

def parseargs(required_args=True):
	class formatter(
		argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
	):
		pass

	epilog = ""
	parser = argparse.ArgumentParser(
		description="Run neighborhood for a given cell type",
		epilog=epilog,
		formatter_class=formatter,
	)
	readable = argparse.FileType("r")

	parser.add_argument(
		"--candidate_regions",
		required=required_args,
		help="Bed file containing candidate_enhancer_regions",
	)
	parser.add_argument(
		"--regions_id",
		required=required_args,
		help="description of input regions",
	)
	parser.add_argument(
		"--outdir",
		required=required_args,
		help="Directory to write Neighborhood files to.",
	)
	parser.add_argument(
		"--assay_name",
		default="",
		nargs="?",
		help="name of assay",
	)
	parser.add_argument(
		"--assay_files",
		default="",
		nargs="?",
		help="comma delimited string of bam files for this assay",
	)
	parser.add_argument(
		"--chrom_sizes",
		required=required_args,
		help="Genome file listing chromosome sizes",
	)
	parser.add_argument(
		"--chrom_sizes_bed",
		required=required_args,
		help="Associated .bed file of chrom_sizes",
	)
	parser.add_argument(
		"--outfile",
		required=required_args,
		help="output file with assay read counts",
	)


	# replace textio wrapper returned by argparse with actual filename
	args = parser.parse_args()
	for name, val in vars(args).items():
		if hasattr(val, "name"):
			setattr(args, name, val.name)
	print(args)
	return args

def count_reads_for_assay(args):
	features = {}
	features[args.assay_name] = args.assay_files.split(" ")

	candidate_regions = read_bed(args.candidate_regions, extra_colnames=["name"])
	candidate_regions = count_features_for_bed(
		df = candidate_regions,
		bed_file = args.candidate_regions,
		genome_sizes = args.chrom_sizes,
		genome_sizes_bed = args.chrom_sizes_bed,
		features = features,
		directory = args.outdir,
		filebase = args.regions_id,
		skip_rpkm_quantile = True,
		use_fast_count = True
	)

	candidate_regions.to_csv(args.outfile, sep = "\t", index = False, header = True, float_format="%.6f")

def main(args):
	count_reads_for_assay(args)


if __name__ == "__main__":
	args = parseargs()
	main(args)