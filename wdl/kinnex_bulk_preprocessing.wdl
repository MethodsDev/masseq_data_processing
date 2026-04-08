version 1.0

import "tasks/pbtools.wdl" as PB

workflow kinnex_bulk_preprocessing {
    meta {
        description: "End-to-end bulk MASseq pipeline: Skera split + QC plots -> Lima demux + Refine -> Merge replicates"
    }

    input {
        # --- Skera inputs ---
        File   input_bam
        File   mas_adapters_fasta
        Int    arraysize
        String? sample_id

        # --- Lima / Demux inputs ---
        File    bulk_barcodes_fasta
        Boolean trimPolyA
        Boolean clipAdapters

        # --- Merge inputs ---
        File    barcode_to_sample
        String? datasetId
        Boolean mergeBams

        # --- Shared ---
        String gcs_output_dir
        Int    num_threads = 16
    }

    # Step 1: Skera split + QC plots
    call PB.pbSkerawQC {
        input:
            hifi_bam           = input_bam,
            sample_id          = sample_id,
            arraysize          = arraysize,
            mas_adapters_fasta = mas_adapters_fasta,
            num_threads        = num_threads,
            gcs_output_dir     = gcs_output_dir
    }

    # Step 2: Lima demux + Isoseq Refine
    call PB.pbLimaBulk {
        input:
            skera_bam           = pbSkerawQC.skera_bam,
            sample_id           = sample_id,
            bulk_barcodes_fasta = bulk_barcodes_fasta,
            trimPolyA           = trimPolyA,
            clipAdapters        = clipAdapters,
            num_threads         = num_threads,
            gcs_output_dir      = gcs_output_dir
    }

    # Step 3: Merge replicates
    call PB.bulkMerge {
        input:
            refine_bampath    = pbLimaBulk.refine_out,
            lima_dir          = pbLimaBulk.lima_out,
            barcode_to_sample = barcode_to_sample,
            datasetId         = datasetId,
            mergeBams         = mergeBams,
            num_threads       = num_threads,
            gcs_output_dir    = gcs_output_dir
    }

    output {
        # From Step 1
        File   skera_bam  = pbSkerawQC.skera_bam
        String QC_plots   = pbSkerawQC.QC_plots

        # From Step 2
        String refine_out = pbLimaBulk.refine_out
        String lima_out   = pbLimaBulk.lima_out

        # From Step 3
        String merge_out  = bulkMerge.merge_out
    }
}
