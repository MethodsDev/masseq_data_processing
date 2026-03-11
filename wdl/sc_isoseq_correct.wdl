version 1.0

workflow isoseqCorrect {
    input {
        File refine_bam
        String? sample_id
        File barcodes_list
        Int num_threads
        String gcs_output_dir
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? boot_disk_size_gb
    }

    call pbCorrect {
        input:
            refine_bam = refine_bam,
            sample_id = sample_id,
            barcodes_list = barcodes_list,
            num_threads = num_threads,
            gcs_output_dir = gcs_output_dir,
            mem_gb = mem_gb,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb,
            boot_disk_size_gb = boot_disk_size_gb
    }

    output {
        File corrected_reads = pbCorrect.corrected_reads
    }
}

task pbCorrect {
    meta {
        description: "Given a refine.bam, runs isoseq correct."
    }

    input {
        File refine_bam
        String? sample_id
        File barcodes_list
        Int num_threads
        String gcs_output_dir
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }

    String resolved_sample_id = select_first([sample_id, sub(basename(refine_bam, ".bam"), ".refine", "")])

    Float input_files_size_gb = 2.5 * (size(refine_bam, "GiB"))
    Int default_ram = 32
    Int default_disk_space_gb = ceil((input_files_size_gb * 5) + 1024)
    Int default_boot_disk_size_gb = 50
    Int machine_mem = if defined(mem_gb) then select_first([mem_gb]) else default_ram

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    command <
        set -euxo pipefail

        echo "Running isoseq correct..."
        isoseq correct --barcodes ~{barcodes_list} -j ~{num_threads} ~{refine_bam} ~{resolved_sample_id}.corrected.bam
        echo "isoseq correct completed."

        echo "Uploading corrected bams..."
        gsutil -m cp *.corrected* ~{outdir}correct/
        echo "Copying corrected reads completed!"
    >>>

    output {
        File corrected_reads = "~{resolved_sample_id}.corrected.bam"
    }

    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }
}
