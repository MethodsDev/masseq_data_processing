version 1.0

task filterLRReadsUsingSRCBs {
    meta {
        description: "Filters long reads using short read cell barcodes"
    }

    input {
        File input_bam
        File barcode_file
        File filter_script
        String? sample_id

        # Optional runtime inputs
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }

    # Defaults
    Int default_mem_gb = 16
    Int default_cpu = 4
    Int default_boot_disk_size_gb = 50

    # Disk estimation
    Float input_files_size_gb = 3.0 * size(input_bam, "GiB") + size(barcode_file, "GiB") + size(filter_script, "GiB")
    Int calculated_disk_space_gb = ceil(input_files_size_gb + 100)

    # Determine sample ID
    String extracted_sample_id = sub(basename(input_bam, ".bam"), "\\.(CB_sorted|corrected)$", "")
    String final_sample_id = select_first([sample_id, extracted_sample_id])

    command <<<
        set -euxo pipefail

        echo "Starting LR filtering using SR cell barcodes..."
        echo "Sample ID: ~{final_sample_id}"
        echo "Input BAM: ~{input_bam}"
        echo "Barcode file: ~{barcode_file}"
        echo "Threads: ~{select_first([cpu, default_cpu])}"

        python ~{filter_script} \
            -i ~{input_bam} \
            -o ~{final_sample_id}.filtered.corrected.bam \
            -b ~{barcode_file} \
            -r ~{final_sample_id}.SR_CB_qc_report.txt \
            -s ~{final_sample_id} \
            -t ~{select_first([cpu, default_cpu])}

        echo "Filtering completed!"

        echo "Output files:"
        ls -lh ~{final_sample_id}.*
    >>>

    output {
        File filtered_bam   = "~{final_sample_id}.filtered.corrected.bam"
        File qc_report      = "~{final_sample_id}.SR_CB_qc_report.txt"
        File barcode_plot   = "~{final_sample_id}.filtered.corrected.bam.barcode_plot.png"
        File barcode_counts = "~{final_sample_id}.filtered.corrected.bam.barcode_counts.tsv"
    }

    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:latest"
        memory: select_first([mem_gb, default_mem_gb]) + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, calculated_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        cpu: select_first([cpu, default_cpu])
        preemptible: select_first([preemptible_attempts, 0])
    }
}


workflow FilterLRReadsWorkflow {
    meta {
        description: "Workflow to filter long reads using short read cell barcodes"
    }

    input {
        File input_bam
        File barcode_file
        File filter_script
        String? sample_id

        # Optional runtime parameters
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }

    call filterLRReadsUsingSRCBs {
        input:
            input_bam = input_bam,
            barcode_file = barcode_file,
            filter_script = filter_script,
            sample_id = sample_id,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts,
            boot_disk_size_gb = boot_disk_size_gb
    }

    output {
        File filtered_bam   = filterLRReadsUsingSRCBs.filtered_bam
        File qc_report      = filterLRReadsUsingSRCBs.qc_report
        File barcode_plot   = filterLRReadsUsingSRCBs.barcode_plot
        File barcode_counts = filterLRReadsUsingSRCBs.barcode_counts
    }
}
