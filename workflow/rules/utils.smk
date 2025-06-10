import glob
import os
import pandas as pd
MAX_MEM_MB = 250 * 1000  # 250GB

def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

# return list of bam accessions and runtypes for a given experiment accession
def get_metadata_for_rpm_assay(metadata_df, experiment_accession):
	OUTPUT_KEY = {"paired-ended": "alignments", "single-ended": "unfiltered alignments"}

	meta = metadata_df[metadata_df['Experiment accession'] == experiment_accession]

	if len(meta) == 0:
		raise Exception(f"No metadata provided for experiment accession {experiment_accession}.")

	# get runtypes from fastq files
	runtype_dict = dict(zip(
		meta[meta['File type'] == "fastq"]['Technical replicate(s)'],
		meta[meta['File type'] == "fastq"]['Run type']
	))

	# get bam file accessions for each run type (same replicate order)
	file_accessions = []
	for rep, rt in runtype_dict.items():
		file_row = meta[(meta['Output type'] == OUTPUT_KEY[rt]) & (meta['Technical replicate(s)'] == rep)]

		if file_row.empty:
			raise Exception(f"No BAM file found for replicate {rep} with run type {rt} for experiment {experiment_accession}.")

		this_acc = file_row['File accession'].iloc[0]
		file_accessions.append(this_acc)

	return file_accessions, list(runtype_dict.values())

def get_metadata_for_peak_assay(metadata_df, experiment_accession):
    OUTPUT_TYPE_PRIORITY = [
        "Optimal IDR thresholded peaks",
        "IDR thresholded peaks",
        "Replicated peaks",
        "Pseudoreplicated peaks",
        "Peaks"
    ]

    peak_meta = metadata_df[
        (metadata_df["Experiment accession"] == experiment_accession) &
        (metadata_df["File format"] == "bed narrowPeak")
    ].copy()

    if peak_meta.empty:
        raise Exception(f"No narrowPeak file found for experiment {experiment_accession}")

    peak_meta["rep_count"] = [tech_reps.count(",") + 1 for tech_reps in peak_meta["Technical replicate(s)"]]
    peak_meta = peak_meta[peak_meta["rep_count"] == peak_meta["rep_count"].max()]

    for output_type in OUTPUT_TYPE_PRIORITY:
        match = peak_meta[peak_meta["Output type"] == output_type]
        if not match.empty:
            accession = match.iloc[0]["File accession"]
            return accession, output_type

    raise Exception(f"No suitable output type found for experiment {experiment_accession}")

def process_encode_metadata(metadata_file, config, scratch_dir):
    rpm_assays = config["RPM_assays"] + config["RPM_expanded_assays"]
    peak_assays = config["peak_overlap_assays"]

    meta = pd.read_csv(metadata_file, sep="\t")
    meta = meta[meta["Biosample term name"].isin(config["pred_cell_types"])]

    processed_bam_dict = {}
    bam_accessions_dict = {}
    bam_runtypes_dict = {}
    processed_peaks_dict = {}
    peak_accessions_dict = {}

    def snakecase(s):
        return s.lower().replace(" ", "_")

    for biosample in config["pred_cell_types"]:
        this_processed_reads = {}
        this_bam_accessions = {}
        this_bam_runtypes = {}
        this_processed_peaks = {}
        this_peak_accessions = {}

        # RPM assays
        for assay in rpm_assays:
            file_accessions = ""
            runtypes = ""

            if assay in config[biosample].get("processed_files", {}):
                reads = config[biosample]["processed_files"][assay].get("reads")
                if reads:
                    this_processed_reads[assay] = reads
                    this_bam_accessions[assay] = ""
                    this_bam_runtypes[assay] = ""
                    continue

            if assay in config[biosample].get("experiment_accessions", {}):
                experiment_accession = config[biosample]["experiment_accessions"][assay]
                file_accessions, runtypes = get_metadata_for_rpm_assay(meta, experiment_accession)
                target_bam_files = [
                    os.path.join(scratch_dir, biosample, f"{acc}.{RUNTYPE_KEY[rt]}.filtered.sorted.bam")
                    for rt, acc in zip(runtypes, file_accessions)
                ]
                this_processed_reads[assay] = target_bam_files
                this_bam_accessions[assay] = file_accessions
                this_bam_runtypes[assay] = runtypes
            else:
                raise Exception(f"No data specified to calculate {biosample} {assay} RPM.")

        # Peak assays
        for assay in peak_assays:
            if assay in config[biosample].get("processed_files", {}):
                peaks = config[biosample]["processed_files"][assay].get("peaks")
                if peaks:
                    this_processed_peaks[assay] = peaks
                    this_peak_accessions[assay] = ""  # No accession needed
                    continue

            if assay in config[biosample].get("experiment_accessions", {}):
                experiment_accession = config[biosample]["experiment_accessions"][assay]
                accession, output_type = get_metadata_for_peak_assay(meta, experiment_accession)
                peak_filename = f"{accession}.{snakecase(output_type)}.peaks.bed.gz"
                this_processed_peaks[assay] = os.path.join(scratch_dir, biosample, peak_filename)
                this_peak_accessions[assay] = accession
            else:
                raise Exception(f"No data specified to locate peaks for {biosample} {assay}")

        processed_bam_dict[biosample] = this_processed_reads
        bam_accessions_dict[biosample] = this_bam_accessions
        bam_runtypes_dict[biosample] = this_bam_runtypes
        processed_peaks_dict[biosample] = this_processed_peaks
        peak_accessions_dict[biosample] = this_peak_accessions

    return (
        processed_bam_dict,
        bam_accessions_dict,
        bam_runtypes_dict,
        processed_peaks_dict,
        peak_accessions_dict
    )


