version 1.0

import "tasks/pbtools.wdl" as PB

workflow pb_sc_lima_isoseq {
    input{
        File skera_bam
        String sample_id
        File primer_fasta
        File barcodes_list
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
            barcodes_list           = barcodes_list,
            read_design             = read_design, 
            num_threads             = num_threads,
            gcs_output_dir          = gcs_output_dir
    }
    output{
        String corrected_reads        = pbSingleCell.corrected_reads
    }
}

