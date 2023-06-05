version 1.0

task pbSkera {
    meta {
        description: "Given hifi reads, spilts MAS 8x or 16x array structure using provided adapters"
    }
# ------------------------------------------------
#Inputs required
    input {
        # Required:
        File hifi_bam
        File mas_adapters_fasta
        Int num_threads
        String gcs_output_dir
        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(hifi_bam, "GiB"))
    Int default_ram = 8
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    command <<<
        set -euxo pipefail

        echo ~{outdir}/skera/$sampleid.skera.bam
        sampleid=`echo ~{hifi_bam} | awk -F '/' '{print $NF}' | awk -F '.hifi_reads' '{print $1}'`
        skera split -j ~{num_threads} ~{hifi_bam} ~{mas_adapters_fasta} $sampleid.skera.bam
        echo "Copying skera out to gcs path provided..."
        gsutil -m cp $sampleid.skera.bam ~{outdir}/skera/

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String skera_out        = "~{outdir}/skera/$sampleid.skera.bam"
    }

# ------------------------------------------------
# Runtime settings:
    runtime {
    docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag1"
    memory: machine_mem + " GiB"
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
    preemptible: select_first([preemptible_attempts, 0])
    cpu: select_first([cpu, 2])
    }

}


