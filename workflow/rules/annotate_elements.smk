### RULES
# return list of candidate elements overlapping peaks for this assay 
rule make_chr_sizes_bed:
	input:
		config["chr_sizes"]
	output: 
		os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	resources:
		mem_mb=1*1000,
		runtime = 10
	shell:
		"""
			awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input} > {output}
		"""

# split elements by cell type and get candidate elements for each
# resize for baseline counts if necessary (but keep element_name original)
rule split_elements_by_cell_type:
	input:
		elements_file = lambda wildcards: config["elements_to_annotate"][wildcards.data_cat]
	params:
		chr_col = lambda wildcards: config["chr_columns"][wildcards.data_cat],
		start_col = lambda wildcards: config["start_columns"][wildcards.data_cat],
		end_col = lambda wildcards: config["end_columns"][wildcards.data_cat],
		ct_col = lambda wildcards: config["cell_type_columns"][wildcards.data_cat],
		trim_size = config["element_trim_size"]
	output:
		element_subset = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "cell_type_subset.tsv.gz"),
		elements_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements.bed4")
	resources:
		mem_mb=32*1000,
		runtime = 60
	shell:
		"""
			set +o pipefail;

			# Determine the command to read the input file
			if [[ {input.elements_file} == *.gz ]]; then
				READ_CMD="zcat {input.elements_file}"
			else
				READ_CMD="cat {input.elements_file}"
			fi

			# Write header
			eval "$READ_CMD" | head -n 1 | gzip > {output.element_subset}

			# Cell-type-specific or generic subset
			if [[ "{params.ct_col}" == "NONE" ]]; then
				eval "$READ_CMD" | sed 1d | gzip >> {output.element_subset}
			else
				col_number=$(eval "$READ_CMD" | head -n 1 | tr '\t' '\n' | nl -v 1 | grep -w "{params.ct_col}" | awk '{{print $1}}')
				eval "$READ_CMD" | awk -F'\t' -v col=$col_number '$col == "{wildcards.this_cell_type}"' | \
					gzip >> {output.element_subset}
			fi

			# Save BED file of filtered subset (chr, start, end, element_name)
			zcat {output.element_subset} | \
				csvtk cut -t -f {params.chr_col},{params.start_col},{params.end_col} | \
				sed 1d | sort -k1,1 -k2,2n | uniq | \
				awk -F'\t' 'BEGIN {{OFS="\t"}} {{print $0, $1 ":" $2 "-" $3}}' | \
				awk -F'\t' 'BEGIN {{OFS="\t"}} {{ col1=$1; col2=$2 + {params.trim_size}; col3=$3 - {params.trim_size};
					print $1, col2, col3, $4}}' > {output.elements_bed}
    """


# outputs original name of prediction ovelrapping peak
rule get_peak_overlaps:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements.bed4"),
		bed_file = lambda wildcards: PROCESSED_PEAKS_DICT[wildcards.this_cell_type][wildcards.assay]
	params:
		peak_ext = lambda wildcards: config["peak_ext_size"][wildcards.assay],
		chr_sizes = config["chr_sizes"]
	output:
		interm_peaks = os.path.join(SCRATCH_DIR,"elements", "{data_cat}", "{this_cell_type}", "{assay}", "peaks.bed"), 
		this_overlap = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "{assay}", "peak_overlap.txt")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb = 32*1000,
		runtime = 30
	shell:
		"""
			if [[ {input.bed_file} == *.gz ]]; then
				zcat {input.bed_file} | cut -f 1-3 | bedtools slop -b {params.peak_ext} -g {params.chr_sizes} | sort -k1,1 -k2,2n > {output.interm_peaks}
			else
				cat {input.bed_file} | cut -f 1-3 | bedtools slop -b {params.peak_ext} -g {params.chr_sizes} | sort -k1,1 -k2,2n > {output.interm_peaks}
			fi
			
			bedtools intersect -wa -a {input.enh_bed} -b {output.interm_peaks} | \
				cut -f4 > {output.this_overlap}
		"""

# output is bed4 with trim + expansion applied; name is original
rule expand_candidate_elements: 
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements.bed4")
	params:
		pred_ext = config["element_ext_size"],
		chr_sizes = config["chr_sizes"]
	output:
		enh_expanded = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements_expanded.bed4")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb = 32*1000,
		runtime = 30
	shell:
		"""
			cat {input.enh_bed} | bedtools slop -b {params.pred_ext} -g {params.chr_sizes} > {output.enh_expanded}
		"""

rule count_reads_baseline_elements_predictions:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements.bed4"),
		assay_bam = lambda wildcards: PROCESSED_BAM_DICT[wildcards.this_cell_type][wildcards.assay],
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		assay_counts = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "{assay}", "baseline_regions.RPM.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""	
			mkdir -p {params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay}		
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "candidate_regions" \
				--outdir {params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay} \
				--assay_name {wildcards.assay} \
				--assay_files "{input.assay_bam}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_counts}
		"""

rule calculate_fold_change_signal_baseline_elements:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements.bed4"),
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		assay_bw = lambda wildcards: config[wildcards.this_cell_type]["fold_change_bws"][wildcards.assay],
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		bw_download = temp(os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "{assay}", "fold_change_bw.bigWig")),
		assay_fc = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "{assay}", "baseline_regions.fold_change.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""
		    curl -L {params.assay_bw} > {output.bw_download}
			mkdir -p {params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay}
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "candidate_regions" \
				--outdir "{params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay}" \
				--assay_name "{wildcards.assay}_fc" \
				--assay_files "{output.bw_download}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_fc}
		"""

rule count_reads_expanded_elements:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements_expanded.bed4"),
		assay_bam = lambda wildcards: PROCESSED_BAM_DICT[wildcards.this_cell_type][wildcards.assay],
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		assay_counts = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "{assay}", "expanded_regions.RPM.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""
			mkdir -p {params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay}_expanded
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "expanded_candidate_regions" \
				--outdir {params.scratch_dir}/{wildcards.this_cell_type}/{wildcards.assay}_expanded \
				--assay_name "{wildcards.assay}" \
				--assay_files "{input.assay_bam}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_counts}

		"""


rule gather_annotations_for_cell_type:
	input:
		enh_list = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "cell_type_subset.tsv.gz"),
		peak_overlaps = expand(os.path.join(SCRATCH_DIR, "elements", "{{data_cat}}", "{{this_cell_type}}", "{assay}", "peak_overlap.txt"), assay = config["peak_overlap_assays"]),
		rpm_original = expand(os.path.join(SCRATCH_DIR, "elements", "{{data_cat}}", "{{this_cell_type}}", "{assay}", "baseline_regions.RPM.tsv"), assay = config["RPM_assays"]),
		rpm_expanded = expand(os.path.join(SCRATCH_DIR, "elements", "{{data_cat}}", "{{this_cell_type}}", "{assay}", "expanded_regions.RPM.tsv"), assay = config["RPM_expanded_assays"]),
		fc_original = expand(os.path.join(SCRATCH_DIR, "elements", "{{data_cat}}", "{{this_cell_type}}", "{assay}", "baseline_regions.fold_change.tsv"), assay = config["FC_assays"])
	params:
		peak_overlap_assays = config["peak_overlap_assays"],
		rpm_assays = config["RPM_assays"],
		rpm_expanded_assays = config["RPM_expanded_assays"],
		fc_assays = config["FC_assays"],
		chr_col = lambda wildcards: config["chr_columns"][wildcards.data_cat],
		start_col = lambda wildcards: config["start_columns"][wildcards.data_cat],
		end_col =  lambda wildcards: config["end_columns"][wildcards.data_cat]
	output:
		enh_plus = os.path.join(SCRATCH_DIR, "elements", "{data_cat}", "{this_cell_type}", "elements_annotated.tsv.gz")
	resources:
		mem_mb=32*1000
	script:
		os.path.join(SCRIPTS_DIR, "gather_annotations_for_one_cell_type.R")

rule merge_annotations_across_cell_types:
	input:
		enh_plus = lambda wildcards: expand(os.path.join(SCRATCH_DIR, "elements", "{{data_cat}}", "{this_cell_type}", "elements_annotated.tsv.gz"),
			this_cell_type = config["element_cell_types"][wildcards.data_cat])
	params:
		cell_types = lambda wildcards: config["element_cell_types"][wildcards.data_cat],
		cell_type_col = lambda wildcards: config["cell_type_columns"][wildcards.data_cat]
	output:
		elements_final = os.path.join(RESULTS_DIR, "{data_cat}.chromatin_annotations.tsv")
	resources:
		mem_mb=32*1000
	script:
		os.path.join(SCRIPTS_DIR, "gather_across_cell_types.R")
