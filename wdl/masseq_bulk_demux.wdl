version 1.0

import "tasks/pbtools.wdl" as PB

workflow pb_masseq_workflow {
    input{
        File input_bam
        String sample_id
        File bulk_barcodes_fasta
        Boolean trimPolyA
        String gcs_output_dir

        # Optional:
        Int num_threads = 20
        #Int? mem_gb
        #Int? preemptible_attempts
        #Int? disk_space_gb
        #Int? cpu
        #Int? boot_disk_size_gb
    }
    call PB.pbLimaBulk{
        input:
            skera_bam               = input_bam,
            sample_id               = sample_id,
            trimPolyA               = trimPolyA,
            bulk_barcodes_fasta     = bulk_barcodes_fasta,
            num_threads             = num_threads,
            gcs_output_dir          = gcs_output_dir
    }
    output{
        String out_path         = pbLimaBulk.demux_out
    }
}