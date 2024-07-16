version 1.0

import "tasks/pbtools.wdl" as PB

workflow merge_replicates {
    input{
        String refine_out
        String lima_out
        File barcode_to_sample
        String? datasetId
        Boolean mergeBams
        String gcs_output_dir

        # Optional:
        Int num_threads = 16

    }
    call PB.bulkMerge{
        input:
            refine_bampath      = refine_out,
            lima_dir            = lima_out,
            barcode_to_sample   = barcode_to_sample,
            datasetId           = datasetId,
            mergeBams           = mergeBams,
            num_threads         = num_threads,
            gcs_output_dir      = gcs_output_dir
    }
    output{
        String merge_out        = bulkMerge.merge_out 
    }
}
