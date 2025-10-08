version 1.0

import "tasks/pbtools.wdl" as PB

workflow CombinedSkeraSingleCell {
    meta {
        description: "Combined workflow that runs pbSkerawQC and pbSingleCell together, capturing all outputs"
    }
    
    input {
        # Required inputs
        File input_bam
        File mas_adapters_fasta
        File primer_fasta
        File barcodes_list
        String sample_id
        String read_design
        String gcs_output_dir
        Int arraysize
        
        # Optional inputs
        Int num_threads = 16
        Boolean trimPolyA = true
        Boolean clipAdapters = true
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? boot_disk_size_gb
    }
    
    # Call pbSkerawQC first
    call PB.pbSkerawQC {
        input:
            hifi_bam = input_bam,
            sample_id = sample_id,
            arraysize = arraysize,
            mas_adapters_fasta = mas_adapters_fasta,
            num_threads = num_threads,
            gcs_output_dir = gcs_output_dir,
            mem_gb = mem_gb,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb,
            boot_disk_size_gb = boot_disk_size_gb
    }
    
    # Call pbSingleCell using the output from pbSkerawQC
    call PB.pbSingleCell {
        input:
            skera_bam = pbSkerawQC.skera_bam,
            sample_id = sample_id,
            primer_fasta = primer_fasta,
            barcodes_list = barcodes_list,
            read_design = read_design,
            trimPolyA = trimPolyA,
            clipAdapters = clipAdapters,
            num_threads = num_threads,
            gcs_output_dir = gcs_output_dir,
            mem_gb = mem_gb,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb,
            boot_disk_size_gb = boot_disk_size_gb
    }
    
    output {
        # Outputs from pbSkerawQC
        File skera_bam = pbSkerawQC.skera_bam
        String QC_plots = pbSkerawQC.QC_plots
        
        # Outputs from pbSingleCell
        File corrected_reads = pbSingleCell.corrected_reads
    }
}
