version 1.0

import "tasks/pbtools.wdl" as PB

workflow masseq_bulk_demux_main {
    input{
        File input_bam
        String sample_id
        File bulk_barcodes_fasta
        Boolean trimPolyA
        Boolean clipAdapters
        String gcs_output_dir

        # Optional:
        Int num_threads = 16

    }
    call PB.pbLimaBulk{
        input:
            skera_bam               = input_bam,
            sample_id               = sample_id,
            bulk_barcodes_fasta     = bulk_barcodes_fasta,
            trimPolyA               = trimPolyA,
            clipAdapters            = clipAdapters,
            num_threads             = num_threads,
            gcs_output_dir          = gcs_output_dir
    }
    output{
        String out_path         = pbLimaBulk.demux_out
    }
}