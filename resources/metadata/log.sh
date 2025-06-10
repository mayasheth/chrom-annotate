# downloaded K562 ChIP-seq metadata

curl -L -o ENCODE_K562_ChIP_metadata.tsv  "https://www.encodeproject.org/metadata/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&biosample_ontology.term_name=K562&files.file_type=bam&files.file_type=fastq&files.file_type=bed+narrowPeak&assay_title=TF+ChIP-seq&assay_title=Histone+ChIP-seq&control_type%21=%2A"

# manually selected relevant columns and filtered to GRCh38, no perturbations, etc. to produce `filtered_ENCODE_K562_ChIP_metadata.tsv`
