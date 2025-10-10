version 1.0

import "isoseq_bcstats.wdl" as BCSTATS
import "tasks/pbtools.wdl" as PB

# Task for merging without bcstats (extracted from mergeAndbcstats)
task MergeOnly {
    meta {
        description: "Merges corrected BAM files"
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
    
    # Calculate input file sizes for disk space estimation
    Float input_files_size_gb = 3.0 * size(corrected_reads, "GiB")
    Int calculated_disk_space_gb = ceil((3 * input_files_size_gb) + 200)
    Int default_boot_disk_size_gb = ceil((2.5 * input_files_size_gb)) 

    command <<<
        set -euxo pipefail
    
        echo "Starting merge analysis..."
        echo "Sample ID: ~{sample_id}"
        echo "Number of input BAMs: ~{length(corrected_reads)}"
        echo "Movie names: ~{sep=',' movie_names}"
    
        # Create output directories
        mkdir ./tmp
        mkdir -p ./data/corrected_merged_bams
    
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
        samtools sort -@ ~{num_threads} -t CB -T ./tmp -o ./data/corrected_merged_bams/~{sample_id}.CB_sorted.bam temp_merged.bam
        
        # Index the sorted BAM
        echo "Indexing sorted BAM..."
        samtools index ./data/corrected_merged_bams/~{sample_id}.CB_sorted.bam
        
        echo "Merge completed!"
        
        # Copy outputs to working directory
        cp ./data/corrected_merged_bams/~{sample_id}.CB_sorted.bam ./
        cp ./data/corrected_merged_bams/~{sample_id}.CB_sorted.bam.bai ./
        
        echo "All outputs copied to working directory"
        ls -lh
    >>> 
    
    output {
        File merged_sorted_bam = "~{sample_id}.CB_sorted.bam"
        File merged_sorted_bam_index = "~{sample_id}.CB_sorted.bam.bai"
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

# Task for counting read lengths
task CountReadLengthsTask {
    input {
        File input_bam
        File? input_bam_index
        File python_script
        String output_basename
        Int workers
        
        # Runtime parameters
        Int cpu = 4
        Int memory_gb = 8
        Int? disk_space_gb
        String docker = "quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3"
    }

    Int default_disk_space_gb = ceil(1.2 * size(input_bam, "GB")) 
    String output_filename = "${output_basename}.pkl"

    command <<<
        set -euox pipefail
        
        # Copy the input Python script to the working directory
        cp "~{python_script}" count_read_lengths.py
        
        # Make the script executable
        chmod +x count_read_lengths.py
        
        # Run the script with mapped reads only (no --include-unmapped flag)
        python3 count_read_lengths.py \
            --bam "~{input_bam}" \
            --out "~{output_filename}" \
            --workers ~{workers}
    >>>

    output {
        File read_length_counts = output_filename
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_gb} GB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
    }
}

# Task for counting UMIs per barcode
task CountUMI {
    input {
        File bam
        File bai
        File process_barcodes_py # Python script
        String sample_id
        String cb_tag = "CB"
        String umi_tag = "XM"
    }

    command <<<
        set -euox pipefail

        # Install required Python packages
        pip install pysam pandas

        cp ~{process_barcodes_py} process_barcodes_py.py

        python3 process_barcodes_py.py \
            ~{bam} \
            --sample_id ~{sample_id} \
            --cb_tag ~{cb_tag} \
            --umi_tag ~{umi_tag} \
            -o ~{sample_id}.counts.tsv
    >>>

    output {
        File counts = "~{sample_id}.counts.tsv"
    }

    runtime {
        docker: "python:3.10-slim"
        cpu: 4
        memory: "32G"
        disks: "local-disk 200 SSD"
    }
}

workflow EnhancedMergeWorkflow {
    meta {
        description: "Enhanced merging workflow with conditional bcstats, groupdedup, read length counting, and UMI counting"
    }
    
    input {
        # Required inputs
        Array[File] corrected_reads
        Array[String] movie_names
        String sample_id
        
        # Conditional flags
        Boolean bcstat = false
        Boolean groupdedup = false
        Boolean countReadLengths = false
        Boolean countsPerBarcode = false
        
        # Optional inputs for CountReadLengths
        File? python_script_readlengths
        String output_basename = "read_length_counts"
        Int workers = 4
        
        # Optional inputs for CountUMI
        File? process_barcodes_py
        String cb_tag = "CB"
        String umi_tag = "XM"
        
        # Runtime parameters
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
        Int? num_threads
        
        # Optional inputs for bcstats
        Int? bcstats_percentile
        String? bcstats_method
        
        # Optional inputs for groupdedup
        Boolean keep_non_real_cells = true
        String gcs_output_dir = ""
    }
    
    # Always perform the merge
    call MergeOnly {
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
    
    # Conditionally run CountReadLengths
    if (countReadLengths) {
        call CountReadLengthsTask {
            input:
                input_bam = MergeOnly.merged_sorted_bam,
                input_bam_index = MergeOnly.merged_sorted_bam_index,
                python_script = select_first([python_script_readlengths]),
                output_basename = output_basename,
                workers = workers
        }
    }
    
    # Conditionally run CountUMI
    if (countsPerBarcode) {
        call CountUMI {
            input:
                bam = MergeOnly.merged_sorted_bam,
                bai = MergeOnly.merged_sorted_bam_index,
                process_barcodes_py = select_first([process_barcodes_py]),
                sample_id = sample_id,
                cb_tag = cb_tag,
                umi_tag = umi_tag
        }
    }
    
    # Conditionally run isoseq bcstats
    if (bcstat) {
        call BCSTATS.isoseqBcstats {
            input:
                corrected_reads_bam = MergeOnly.merged_sorted_bam,
                sample_id = sample_id,
                bcstats_percentile = bcstats_percentile,
                bcstats_method = bcstats_method,
                mem_gb = mem_gb,
                disk_space_gb = disk_space_gb,
                cpu = cpu,
                preemptible_attempts = preemptible_attempts,
                boot_disk_size_gb = boot_disk_size_gb
        }
    }
    
    # Conditionally run pbGroupdedup
    if (groupdedup) {
        call PB.pbGroupdedup {
            input:
                input_bam = MergeOnly.merged_sorted_bam,
                sample_id = sample_id,
                keep_non_real_cells = keep_non_real_cells,
                num_threads = select_first([num_threads, 8]),
                gcs_output_dir = gcs_output_dir,
                mem_gb = mem_gb,
                preemptible_attempts = preemptible_attempts,
                disk_space_gb = disk_space_gb,
                boot_disk_size_gb = boot_disk_size_gb
        }
    }
    
    output {
        # Always available outputs
        File merged_sorted_bam = MergeOnly.merged_sorted_bam
        File merged_sorted_bam_index = MergeOnly.merged_sorted_bam_index
        
        # Conditional outputs
        File? read_length_counts = CountReadLengthsTask.read_length_counts
        File? umi_counts = CountUMI.counts
        File? bcstats_json = isoseqBcstats.bcstats_json
        File? bcstats_tsv = isoseqBcstats.bcstats_tsv
        String? dedup_out = pbGroupdedup.dedup_out
        File? deduped_bam = pbGroupdedup.deduped_bam
    }
}
