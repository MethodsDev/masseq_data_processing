version 1.0

import "tasks/pbtools.wdl" as PB

workflow pb_masseq_refine {
    input{
        String input_path
        File primer_fasta
        Boolean trimPolyA
        String gcs_output_dir

        # Optional:
        Int num_threads = 16
        #Int? mem_gb
        #Int? preemptible_attempts
        #Int? disk_space_gb
        Int cpu = 16
        #Int? boot_disk_size_gb
    }
    call PB.pbRefine{
        input:
            input_path              = input_path,
            trimPolyA               = trimPolyA,
            primer_fasta            = primer_fasta,
            num_threads             = num_threads,
            gcs_output_dir          = gcs_output_dir
    }
    output{
        String out_path         = pbRefine.refine_out
    }
}