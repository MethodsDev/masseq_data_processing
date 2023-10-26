version 1.0

import "tasks/pbtools.wdl" as PB

workflow pb_masseq_workflow {
    input{
        File input_bam
        File mas_adapters_fasta
        String gcs_output_dir
        String sample_id
        Int arraysize

        # Optional:
        Int num_threads = 20
        #Int? mem_gb
        #Int? preemptible_attempts
        #Int? disk_space_gb
        #Int? cpu
        #Int? boot_disk_size_gb
    }
    call PB.pbSkerawQC{
        input:
            hifi_bam            = input_bam,
            sample_id           = sample_id,
            arraysize           = arraysize,
            mas_adapters_fasta  = mas_adapters_fasta,
            num_threads         = num_threads,
            gcs_output_dir      = gcs_output_dir
    }
    output{
        String out_path         = pbSkerawQC.skera_out
    }
}