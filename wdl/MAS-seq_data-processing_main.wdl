version 1.2

import "tasks/pbtools.wdl" as PB

workflow pb_masseq_workflow {
    input{
        File input_bam
        File mas_adapters_fasta
        String gcs_output_dir
        String sample_id
        Int arraysize

        # Optional:
        Int num_threads = 16

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
