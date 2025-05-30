version 1.0

import "tasks/pbtools.wdl" as PB

workflow combinedSkeraSingleCell {
    meta {
        description: "Combined workflow that runs pbSkerawQC followed by pbSingleCell"
    }

    input {
        # Inputs for pbSkerawQC
        File hifi_bam
        String? sample_id
        File mas_adapters_fasta
        Int num_threads
        Int arraysize = 8
        String gcs_output_dir

        # Additional inputs for pbSingleCell
        File primer_fasta
        File barcodes_list
        String read_design
        Boolean trimPolyA = true
        Boolean clipAdapters = true

        # Optional runtime parameters
        Int? mem_gb_skera
        Int? mem_gb_singlecell
        Int? preemptible_attempts
        Int? disk_space_gb_skera
        Int? disk_space_gb_singlecell
        Int? boot_disk_size_gb
    }

    # Determine sample ID if not provided
    String resolved_sample_id = select_first([sample_id, sub(basename(hifi_bam,".bam"),".hifi_reads","")])

    # Step 1: Run pbSkerawQC
    call PB.pbSkerawQC {
        input:
            hifi_bam = hifi_bam,
            sample_id = resolved_sample_id,
            mas_adapters_fasta = mas_adapters_fasta,
            num_threads = num_threads,
            arraysize = arraysize,
            gcs_output_dir = gcs_output_dir,
            mem_gb = mem_gb_skera,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb_skera,
            cpu = num_threads,
            boot_disk_size_gb = boot_disk_size_gb
    }

    # Step 2: Run pbSingleCell using the output from pbSkerawQC
    call PB.pbSingleCell {
        input:
            skera_bam = pbSkerawQC.skera_bam,
            sample_id = resolved_sample_id,
            primer_fasta = primer_fasta,
            barcodes_list = barcodes_list,
            read_design = read_design,
            trimPolyA = trimPolyA,
            clipAdapters = clipAdapters,
            num_threads = num_threads,
            gcs_output_dir = gcs_output_dir,
            mem_gb = mem_gb_singlecell,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb_singlecell,
            cpu = num_threads,
            boot_disk_size_gb = boot_disk_size_gb
    }

    output {
        # Outputs from pbSkerawQC
        File skera_bam = pbSkerawQC.skera_bam
        String skera_qc = pbSkerawQC.QC_plots
        
        # Outputs from pbSingleCell
        File corrected_reads = pbSingleCell.corrected_reads
    }
}