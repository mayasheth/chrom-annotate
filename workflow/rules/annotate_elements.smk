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

# input is bed4; output is has resized coordinates BUT original name 
rule resize_prediction_elements:
	input:
		enh_bed = lambda wildcards: os.path.join(config[wildcards.cell_type]["E2G_results"], "Neighborhoods", "EnhancerList.bed")
	params:
		pred_trim_size = config["pred_trim_size"]
	output:
		enh_resized_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "resized_candidate_elements.bed")
	resources:
		mem_mb=8*1000,
		runtime = 30
	shell:
		"""
			awk 'BEGIN {{OFS="\t"}} {{ 
				col1=$1; col2=$2 + {params.pred_trim_size}; col3=$3 - {params.pred_trim_size};
				print $1, col2, col3, $4
			}}' {input.enh_bed} > {output.enh_resized_bed}

		"""

# outputs original name of prediction ovelrapping peak
rule get_peak_overlap_predictions:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "resized_candidate_elements.bed"),
		bed_file = lambda wildcards: PROCESSED_PEAKS_DICT[wildcards.cell_type][wildcards.assay]
	params:
		peak_ext = lambda wildcards: config["peak_ext_size"][wildcards.assay],
		chr_sizes = config["chr_sizes"]
	output:
		interm_peaks = os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "peaks.bed"), 
		this_overlap = os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "peak_overlap.txt")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb = 32*1000,
		runtime = 30
	shell:
		"""
			if [[ {input.bed_file} == *.gz ]]
			then
				zcat {input.bed_file} | cut -f 1-3 | bedtools slop -b {params.peak_ext} -g {params.chr_sizes} | sort -k1,1 -k2,2n > {output.interm_peaks}
			else
				cat {input.bed_file} | cut -f 1-3 | bedtools slop -b {params.peak_ext} -g {params.chr_sizes} | sort -k1,1 -k2,2n > {output.interm_peaks}
			fi
			
			bedtools intersect -wa -a {input.enh_bed} -b {output.interm_peaks} | \
				cut -f4 > {output.this_overlap}
		"""

# output is bed4 with trim + expansion applied; name is original
rule expand_candidate_elements_predictions: 
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "resized_candidate_elements.bed")
	params:
		pred_ext = config["pred_ext_size"],
		chr_sizes = config["chr_sizes"]
	output:
		enh_expanded = os.path.join(SCRATCH_DIR, "{cell_type}", "EnhancerList.expanded.bed4")
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
		enh_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "resized_candidate_elements.bed"),
		assay_bam = lambda wildcards: PROCESSED_BAM_DICT[wildcards.cell_type][wildcards.assay],
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		assay_counts = os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "baseline_regions.RPM.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""			
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "candidate_regions" \
				--outdir {params.scratch_dir}/{wildcards.cell_type}/{wildcards.assay} \
				--assay_name {wildcards.assay} \
				--assay_files "{input.assay_bam}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_counts}
		"""

rule calculate_fold_change_signal_baseline_elements_predictions:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "resized_candidate_elements.bed"),
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		assay_bw = lambda wildcards: config[wildcards.cell_type]["fold_change_bws"][wildcards.assay],
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		bw_download = temp(os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "fold_change_bw.bigWig")),
		assay_fc = os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "baseline_regions.fold_change.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""
		    curl -L {params.assay_bw} > {output.bw_download}
			
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "candidate_regions" \
				--outdir "{params.scratch_dir}/{wildcards.cell_type}/{wildcards.assay}" \
				--assay_name "{wildcards.assay}_fc" \
				--assay_files "{output.bw_download}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_fc}
		"""

rule count_reads_expanded_elements_predictions:
	input:
		enh_bed = os.path.join(SCRATCH_DIR, "{cell_type}", "EnhancerList.expanded.bed4"),
		assay_bam = lambda wildcards: PROCESSED_BAM_DICT[wildcards.cell_type][wildcards.assay],
		chr_sizes = config["chr_sizes"],
		chr_sizes_bed = os.path.join(SCRATCH_DIR, "chr_sizes.bed")
	params:
		scratch_dir = SCRATCH_DIR,
		scripts_dir = SCRIPTS_DIR
	output: 
		assay_counts = os.path.join(SCRATCH_DIR, "{cell_type}", "{assay}", "expanded_regions.RPM.tsv")
	conda:
		config["envs"]["ABC"]
	resources:
		mem_mb=32*1000
	shell:
		"""
			mkdir -p {params.scratch_dir}/{wildcards.cell_type}/{wildcards.assay}_expanded
			python {params.scripts_dir}/run.neighborhoods.extended.py \
				--candidate_regions {input.enh_bed} \
				--regions_id "expanded_candidate_regions" \
				--outdir {params.scratch_dir}/{wildcards.cell_type}/{wildcards.assay}_expanded \
				--assay_name "{wildcards.assay}" \
				--assay_files "{input.assay_bam}" \
				--chrom_sizes {input.chr_sizes} \
				--chrom_sizes_bed {input.chr_sizes_bed} \
				--outfile {output.assay_counts}

		"""


rule gather_annotations_for_predictions:
	input:
		enh_list = lambda wildcards: os.path.join(config[wildcards.cell_type]["E2G_results"], "Neighborhoods", "EnhancerList.txt"),
		peak_overlaps = expand(os.path.join(SCRATCH_DIR, "{{cell_type}}", "{assay}", "peak_overlap.txt"), assay = config["peak_overlap_assays"]),
		rpm_original = expand(os.path.join(SCRATCH_DIR, "{{cell_type}}", "{assay}", "baseline_regions.RPM.tsv"), assay = config["RPM_assays"]),
		rpm_expanded = expand(os.path.join(SCRATCH_DIR, "{{cell_type}}", "{assay}", "expanded_regions.RPM.tsv"), assay = config["RPM_expanded_assays"]),
		fc_original = expand(os.path.join(SCRATCH_DIR, "{{cell_type}}", "{assay}", "baseline_regions.fold_change.tsv"), assay = config["FC_assays"])
	params:
		peak_overlap_assays = config["peak_overlap_assays"],
		rpm_assays = config["RPM_assays"],
		rpm_expanded_assays = config["RPM_expanded_assays"],
		fc_assays = config["FC_assays"],
		chr_col = config["chr_columns"]["predictions"],
		start_col = config["start_columns"]["predictions"],
		end_col = config["end_columns"]["predictions"]
	output:
		enh_plus = os.path.join(RESULTS_DIR, "{cell_type}", "EnhancerList.extended.tsv")
	resources:
		mem_mb=32*1000
	script:
		os.path.join(SCRIPTS_DIR, "gather_annotations_for_one_cell_type.R")