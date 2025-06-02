version 1.0

task mergeAndbcstats {
    meta {
        description: "Merges corrected BAM files and runs isoseq bcstats analysis"
    }
    
    input {
        # Required inputs
        Array[File] corrected_reads
        Array[String] movie_names
        String sample_id
        
        # Optional inputs
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
        Int num_threads = 8
    }
    
    # Set defaults
    Int default_mem_gb = 16
    Int default_cpu = 8
    Int default_boot_disk_size_gb = 50
    
    # Calculate input file sizes for disk space estimation
    Float input_files_size_gb = 3.0 * size(corrected_reads, "GiB")
    Int calculated_disk_space_gb = ceil((3 * input_files_size_gb) + 200)
    
    command <<<
        set -euxo pipefail
    
        echo "Starting merge and bcstats analysis..."
        echo "Sample ID: ~{sample_id}"
        echo "Number of input BAMs: ~{length(corrected_reads)}"
        echo "Movie names: ~{sep=',' movie_names}"
    
        # Create output directories
        mkdir -p /mnt/data/saturation_index/bcstats_out
        mkdir -p /mnt/data/saturation_index/corrected_merged_bams
    
        # List and verify input files
        echo "Input BAM files:"
        for bam_file in ~{sep=' ' corrected_reads}; do
            if [[ -f "$bam_file" ]]; then
                echo "Found: $bam_file ($(du -h "$bam_file" | cut -f1))"
            else
                echo "Missing: $bam_file"
                exit 1
            fi
        done
    
        # Merge BAM files using samtools
        echo "Merging BAM files with samtools..."
        if [ ~{length(corrected_reads)} -eq 1 ]; then
            echo "Single BAM file detected, copying instead of merging..."
            cp ~{corrected_reads[0]} temp_merged.bam
        else
            echo "Multiple BAM files detected, merging..."
            samtools merge -@ ~{num_threads} temp_merged.bam ~{sep=' ' corrected_reads}
        fi
        
        # Sort by cell barcode (CB tag)
        echo "Sorting BAM by cell barcode (CB tag)..."
        samtools sort -@ ~{num_threads} -t CB -o /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam temp_merged.bam
        
        # Index the sorted BAM
        echo "Indexing sorted BAM..."
        samtools index /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam
        
        # Run isoseq bcstats
        echo "Running isoseq bcstats..."
        isoseq bcstats \
            --method percentile \
            --percentile 95 \
            --json /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.json \
            -o /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.tsv \
            /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam 2>&1
        
        echo "bcstats analysis completed!"
        
        # Run UMI counting and saturation analysis scripts
        echo "Running UMI counting and saturation analysis (VLOGdatatransform)..."
        python /usr/local/src/masseq_data_processing/sc_scripts/Isoseq_corrected_bam_umi_counting_and_saturation_VLOGdatatransform.py \
            --bam_file /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam \
            --tsv_output ~{sample_id}.umi_counts_vlog.tsv \
            --bcstats_file /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.tsv \
            --saturation_index_plot ~{sample_id}.saturation_index_vlog.png
        
        echo "Running UMI counting and saturation analysis (vComprehensive)..."
        python /usr/local/src/masseq_data_processing/sc_scripts/Isoseq_corrected_bam_umi_counting_and_saturation_vComprehensive.py \
            --bam_file /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam \
            --tsv_output ~{sample_id}.umi_counts_comprehensive.tsv \
            --bcstats_file /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.tsv \
            --saturation_index_plot ~{sample_id}.saturation_index_comprehensive.png
        
        echo "UMI counting and saturation analysis completed!"
        
        # Copy outputs to working directory for WDL to capture
        cp /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam ./
        cp /mnt/data/saturation_index/corrected_merged_bams/~{sample_id}.CB_sorted.bam.bai ./
        cp /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.json ./
        cp /mnt/data/saturation_index/bcstats_out/~{sample_id}.bcstats.tsv ./
        
        echo "All outputs copied to working directory"
        ls -lh
    >>> 
    
    output {
        File merged_sorted_bam = "~{sample_id}.CB_sorted.bam"
        File merged_sorted_bam_index = "~{sample_id}.CB_sorted.bam.bai"
        File bcstats_json = "~{sample_id}.bcstats.json"
        File bcstats_tsv = "~{sample_id}.bcstats.tsv"
        File umi_counts_vlog_tsv = "~{sample_id}.umi_counts_vlog.tsv"
        File saturation_index_vlog_plot = "~{sample_id}.saturation_index_vlog.png"
        File umi_counts_comprehensive_tsv = "~{sample_id}.umi_counts_comprehensive.tsv"
        File saturation_index_comprehensive_plot = "~{sample_id}.saturation_index_comprehensive.png"
    }
    
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:latest"
        memory: select_first([mem_gb, default_mem_gb]) + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, calculated_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        cpu: select_first([cpu, default_cpu])
        preemptible: select_first([preemptible_attempts, 0])
    }
}

workflow mergeAndbcstatsWorkflow {
    meta {
        description: "Workflow to merge corrected BAM files and perform bcstats analysis"
    }
    
    input {
        Array[File] corrected_reads
        Array[String] movie_names
        String sample_id
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
        Int? num_threads
    }
    
    call mergeAndbcstats {
        input:
            corrected_reads = corrected_reads,
            movie_names = movie_names,
            sample_id = sample_id,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts,
            boot_disk_size_gb = boot_disk_size_gb,
            num_threads = num_threads
    }
    
    output {
        File merged_sorted_bam = mergeAndbcstats.merged_sorted_bam
        File merged_sorted_bam_index = mergeAndbcstats.merged_sorted_bam_index
        File bcstats_json = mergeAndbcstats.bcstats_json
        File bcstats_tsv = mergeAndbcstats.bcstats_tsv
        File umi_counts_vlog_tsv = mergeAndbcstats.umi_counts_vlog_tsv
        File saturation_index_vlog_plot = mergeAndbcstats.saturation_index_vlog_plot
        File umi_counts_comprehensive_tsv = mergeAndbcstats.umi_counts_comprehensive_tsv
        File saturation_index_comprehensive_plot = mergeAndbcstats.saturation_index_comprehensive_plot
    }
}