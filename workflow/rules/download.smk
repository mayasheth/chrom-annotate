
# download bam files for one biosample
rule download_and_process_bam:
    output:
        bam_raw = temp(os.path.join(SCRATCH_DIR, "{biosample}", "{accession}.{rt}.bam")),
        bam_out = os.path.join(SCRATCH_DIR, "{biosample}", "{accession}.{rt}.filtered.sorted.bam"),
        ind_out = os.path.join(SCRATCH_DIR, "{biosample}", "{accession}.{rt}.filtered.sorted.bam.bai")
    resources:
        mem_mb = 32*1000,
        runtime = 75
    shell:
        """
            curl -L https://www.encodeproject.org/files/{wildcards.accession}/@@download/{wildcards.accession}.bam > {output.bam_raw}
            
            if [[ {wildcards.rt} = 'se' ]]; then
                samtools view -F 780 -q 30 -u {output.bam_raw} | samtools sort -o {output.bam_out}
            elif [[ {wildcards.rt} == 'pe' ]]; then
                samtools sort {output.bam_raw} -o {output.bam_out}
            fi

            samtools index {output.bam_out}       
        """



# download bam files for one biosample
rule download_peaks:
    output:
        peaks_out = os.path.join(SCRATCH_DIR, "{biosample}", "{accession}.{output_type}.peaks.bed.gz")
    shell:
        """
            curl -L https://www.encodeproject.org/files/{wildcards.accession}/@@download/{wildcards.accession}.bed.gz > {output.peaks_out}
    
        """