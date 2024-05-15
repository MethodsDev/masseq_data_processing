version 1.0

import "tasks/pbtools.wdl" as PB

workflow pb_sc_lima_isoseq {
    input{
        String skera_bam
        String sample_id
        File primer_fasta
        String 10x_barcodes_list
        String read_design
        Boolean trimPolyA = true
        Boolean clipAdapters = true 
        String gcs_output_dir

        # Optional:
        Int num_threads = 16

    }
    call PB.pbSingleCell{
        input:
            skera_bam               = skera_bam,
            sample_id               = sample_id,
            trimPolyA               = trimPolyA,
            clipAdapters            = clipAdapters, 
            primer_fasta            = primer_fasta,
            10x_barcodes_list       = 10x_barcodes_list,
            read_design             = read_design, 
            num_threads             = num_threads,
            gcs_output_dir          = gcs_output_dir
    }
    output{
        String corrected_reads        = pbSingleCell.corrected_reads
    }
}

