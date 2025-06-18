version 1.0

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
    String extracted_sample_id = sub(basename(corrected_reads_bam, ".bam"), "\\.(CB_sorted|corrected)$", "")
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

workflow IsoseqBcstatsWorkflow {
    meta {
        description: "Workflow to perform isoseq bcstats analysis on CB-sorted BAM files"
    }
    
    input {
        File corrected_reads_bam
        String sample_id
        Int? percentile
        String? method
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }
    
    call isoseqBcstats {
        input:
            corrected_reads_bam = corrected_reads_bam,
            sample_id = sample_id,
            percentile = percentile,
            method = method,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts,
            boot_disk_size_gb = boot_disk_size_gb
    }
    
    output {
        File bcstats_json = isoseqBcstats.bcstats_json
        File bcstats_tsv = isoseqBcstats.bcstats_tsv
    }
}