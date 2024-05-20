version 1.0

import "tasks/pbtools.wdl" as PB

workflow sc_dedup {
    input{
        File input_bam
        String sample_id
        Boolean keep_non_real_cells = true
        String gcs_output_dir

        # Optional:
        Int num_threads = 16

    }
    call PB.pbGroupdedup{
        input:
            input_bam            = input_bam,
            sample_id            = sample_id,
            keep_non_real_cells  = keep_non_real_cells,
            num_threads          = num_threads,
            gcs_output_dir       = gcs_output_dir
    }
    output{
        String dedup_path        = pbGroupdedup.dedup_out
        File deduped_bam         = pbGroupdedup.deduped_bam
    }
}
