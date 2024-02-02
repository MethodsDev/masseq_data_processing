version 1.0

import "tasks/other_utils.wdl" as utils

workflow mergeBulkByBarcodes {
    input {
        String input_dir_path
        Boolean rename
        String gcs_output_dir

        # Optional:
        Int num_threads = 16

    }
    call utils.mergeBulk {
        input:
            input_path           = input_dir_path,
            rename               = rename,
            num_threads          = num_threads,
            gcs_output_dir       = gcs_output_dir
    }
    output {
        String gcs_output_dir       = mergeBulk.merge_out
 #       File monitoringLog          = "monitoring.log"

    }
}
