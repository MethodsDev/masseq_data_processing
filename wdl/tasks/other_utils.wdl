version 1.0

task mergeBulk {
    meta {
    description: "For a bulk workflow, given de-mutliplexed reads, collapses reads across barcodes together uing samtools merge."
}
# ------------------------------------------------
#Inputs required
input {
    # Required:
    String input_dir
    Boolean rename = false
    File? rename_key # json with barcode and sample name pairing
    Int num_threads
    String gcs_output_dir
    File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"

    # Optional:
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int cpu = num_threads
    Int? boot_disk_size_gb
}
# Computing required disk size
Int default_ram = 16
Int default_disk_space_gb = 40
Int default_boot_disk_size_gb = 40

# Mem is in units of GB
Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

command <<<
    set -euxo pipefail

    bash ~{monitoringScript} > monitoring.log &

    gsutil -m cp ~{input_dir}/* .
    echo "Copying input bams completed!"

    echo "Merging across barcodes..."
    declare -a bc_indexes=('bc01' 'bc02' 'bc03' 'bc04' 'bc05' 'bc06' 'bc07' 'bc08' 'bc09' 'bc10' 'bc11' 'bc12')
    wd=~{input_dir}
    for  bc in ${bc_indexes[@]}
    do
        samtools merge $wd/merge/$bc.bam `ls $wd/*$bc*.bam`
        samtools sort -@ ~{cpu} $wd/merge/$bc.bam > $wd/merge/$bc.merged.sorted.bam
        samtools index -@ ~{cpu} $wd/merge/$bc.merged.sorted.bam
    done
    echo "Uploading merged bams..."
    gsutil -m cp $wd/merge/* ~{outdir}merge/

>>>
# ------------------------------------------------
# Outputs:
output {
    # Default output file name:
    String merge_out        = "~{outdir}merge"
    File monitoringLog      = "monitoring.log"
}

# ------------------------------------------------
# Runtime settings:
runtime {
    docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag4"
    memory: machine_mem + " GiB"
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
    bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
    preemptible: select_first([preemptible_attempts, 0])
    cpu: cpu
} 
