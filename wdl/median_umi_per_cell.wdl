version 1.0

task medianUmiPerCell {
    meta {
        description: "Runs median_umi_per_cell.py on a BAM file"
    }

    input {
        File input_bam
        File script_file
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
    Int default_cpu = 8
    Int default_boot_disk_size_gb = 50

    # Disk estimation
    Float input_files_size_gb = 2.0 * size(input_bam, "GiB") + size(script_file, "GiB")
    Int calculated_disk_space_gb = ceil(input_files_size_gb + 50)

    # Sample ID / prefix
    String extracted_sample_id = sub(basename(input_bam, ".bam"), "\\.(CB_sorted|corrected|filtered)$", "")
    String final_sample_id = select_first([sample_id, extracted_sample_id])

    command <<<
        set -euxo pipefail

        echo "Starting median UMI per cell calculation..."
        echo "Sample ID: ~{final_sample_id}"
        echo "Input BAM: ~{input_bam}"
        echo "Script: ~{script_file}"
        echo "Threads: ~{select_first([cpu, default_cpu])}"

        /usr/bin/time -v python ~{script_file} \
            -i ~{input_bam} \
            -t ~{select_first([cpu, default_cpu])}

        echo "Run completed."
        echo "Files in working directory:"
        ls -lh
    >>>

    output {
        # Adjust these names if your script writes different outputs
        File median_umi_tsv = "~{final_sample_id}.median_umi_per_cell.tsv"
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

workflow MedianUmiPerCellWorkflow {
    meta {
        description: "Workflow to run median_umi_per_cell.py on a BAM file"
    }

    input {
        File input_bam
        File script_file
        String? sample_id

        # Optional runtime parameters
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }

    call medianUmiPerCell {
        input:
            input_bam = input_bam,
            script_file = script_file,
            sample_id = sample_id,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts,
            boot_disk_size_gb = boot_disk_size_gb
    }

    output {
        File median_umi_tsv = medianUmiPerCell.median_umi_tsv
    }
}
