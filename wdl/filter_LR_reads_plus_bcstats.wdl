version 1.0

task filterLRReadsUsingSRCBs {
    meta {
        description: "Filters long reads using short read cell barcodes"
    }
    
    input {
        File input_bam
        File barcode_file  # File containing short read cell barcodes
        File filter_script  # Python script for filtering
        String? sample_id
        
        # Optional inputs
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }
    
    # Set defaults
    Int default_mem_gb = 16
    Int default_cpu = 4
    Int default_boot_disk_size_gb = 50
    
    # Calculate input file sizes for disk space estimation
    Float input_files_size_gb = 3.0 * size(input_bam, "GiB") + size(barcode_file, "GiB") + size(filter_script, "GiB")
    Int calculated_disk_space_gb = ceil(input_files_size_gb + 100)

    # Extract sample ID from BAM filename if not provided
    String extracted_sample_id = sub(basename(input_bam, ".bam"), "\\.(CB_sorted|corrected)$", "")
    String final_sample_id = select_first([sample_id, extracted_sample_id])

    command <<<
        set -euxo pipefail
        
        echo "Starting LR filtering using SR cell barcodes..."
        echo "Sample ID: ~{final_sample_id}"
        echo "Input BAM: ~{input_bam}"
        echo "Barcode file: ~{barcode_file}"
        echo "Filter script: ~{filter_script}"
        
        # Run the filtering script
        echo "Running filtering script..."
        python ~{filter_script} \
            -i ~{input_bam} \
            -o ~{final_sample_id}.filtered.corrected.bam \
            -b ~{barcode_file} \
            -r ~{final_sample_id}.SR_CB_qc_report.txt \
            -s ~{final_sample_id}
        
        echo "Filtering completed!"
        
        echo "Output files:"
        ls -lh ~{final_sample_id}.*
    >>>
    
    output {
        File filtered_bam = "~{final_sample_id}.filtered.corrected.bam"
        File qc_report = "~{final_sample_id}.SR_CB_qc_report.txt"
        File barcode_plot = "~{final_sample_id}.filtered.corrected.bam.barcode_plot.png"
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

task isoseqBcstats {
    meta {
        description: "Performs isoseq bcstats analysis on CB-sorted BAM files"
    }
    
    input {
        File corrected_reads_bam
        String? sample_id
        
        # Optional inputs
        Int percentile = 95
        String method = "percentile"
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }
    
    # Set defaults
    Int default_mem_gb = 8
    Int default_cpu = 4
    Int default_boot_disk_size_gb = 50
    
    # Calculate input file sizes for disk space estimation
    Float input_files_size_gb = 2.0 * size(corrected_reads_bam, "GiB")
    Int calculated_disk_space_gb = ceil(input_files_size_gb + 100)

    # Extract sample ID from BAM filename if not provided
    String extracted_sample_id = sub(basename(corrected_reads_bam, ".bam"), "\\.(CB_sorted|corrected|filtered)$", "")
    String final_sample_id = select_first([sample_id, extracted_sample_id])

    command <<<
        set -euxo pipefail
        
        echo "Starting isoseq bcstats analysis..."
        echo "Sample ID: ~{final_sample_id}"
        echo "Input BAM: ~{corrected_reads_bam}"
        echo "Method: ~{method}"
        echo "Percentile: ~{percentile}"
        
        # Check isoseq version
        echo "Checking isoseq version..."
        isoseq --version 2>&1
        
        # Run isoseq bcstats
        echo "Running isoseq bcstats..."
        isoseq bcstats \
            --method ~{method} \
            --percentile ~{percentile} \
            --json ~{final_sample_id}.bcstats.json \
            -o ~{final_sample_id}.bcstats.tsv \
            ~{corrected_reads_bam} 2>&1
        
        echo "bcstats analysis completed!"
        
        echo "Output files copied to working directory:"
        ls -lh *.bcstats.*
    >>>
    
    output {
        File bcstats_json = "~{final_sample_id}.bcstats.json"
        File bcstats_tsv = "~{final_sample_id}.bcstats.tsv"
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

workflow CombinedFilterAndBcstatsWorkflow {
    meta {
        description: "Workflow to filter long reads using short read cell barcodes and perform bcstats analysis"
    }
    
    input {
        File input_bam
        File barcode_file
        File filter_script
        String? sample_id
        
        # Filter task parameters
        Int? filter_mem_gb
        Int? filter_disk_space_gb
        Int? filter_cpu
        Int? filter_preemptible_attempts
        Int? filter_boot_disk_size_gb
        
        # Bcstats task parameters
        Int? bcstats_percentile
        String? bcstats_method
        Int? bcstats_mem_gb
        Int? bcstats_disk_space_gb
        Int? bcstats_cpu
        Int? bcstats_preemptible_attempts
        Int? bcstats_boot_disk_size_gb
    }
    
    # Step 1: Filter long reads using short read cell barcodes
    call filterLRReadsUsingSRCBs {
        input:
            input_bam = input_bam,
            barcode_file = barcode_file,
            filter_script = filter_script,
            sample_id = sample_id,
            mem_gb = filter_mem_gb,
            disk_space_gb = filter_disk_space_gb,
            cpu = filter_cpu,
            preemptible_attempts = filter_preemptible_attempts,
            boot_disk_size_gb = filter_boot_disk_size_gb
    }
    
    # Step 2: Run bcstats analysis on the filtered BAM
    call isoseqBcstats {
        input:
            corrected_reads_bam = filterLRReadsUsingSRCBs.filtered_bam,
            sample_id = sample_id,
            percentile = bcstats_percentile,
            method = bcstats_method,
            mem_gb = bcstats_mem_gb,
            disk_space_gb = bcstats_disk_space_gb,
            cpu = bcstats_cpu,
            preemptible_attempts = bcstats_preemptible_attempts,
            boot_disk_size_gb = bcstats_boot_disk_size_gb
    }
    
    output {
        File filtered_bam = filterLRReadsUsingSRCBs.filtered_bam
        File qc_report = filterLRReadsUsingSRCBs.qc_report
        File barcode_plot = filterLRReadsUsingSRCBs.barcode_plot
        File bcstats_json = isoseqBcstats.bcstats_json
        File bcstats_tsv = isoseqBcstats.bcstats_tsv
    }
}