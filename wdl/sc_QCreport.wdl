version 1.0

task generateReport {
    meta {
        description: "Generates an HTML report consolidating flowcell level metrics and saturation plot for single cell runs"
    }
    
    input {
        # Required inputs
        String prefixes
        String barcodes
        String QC_plots_dir  # QC_plots from skera data_dir
        String sample_id
        File base_inserts_run
        File base_inserts_prefix
        Array[File] CCS_report  # Array of CCS report files
        File saturation_plot_png  # Saturation plot PNG file
        File reporting_script  # Path to generate_report.py script
        
        # Optional inputs
        String? output_filename
        String? docker_image
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }
    
    # Set defaults
    Int default_mem_gb = 4
    Int default_disk_space_gb = 100
    Int default_cpu = 1
    Int default_boot_disk_size_gb = 25
    String default_docker = "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/sc_reporting:latest"  
    
    # Determine output filename
    String final_output = select_first([output_filename, "${sample_id}_report.html"])
    
    command <<<
        set -euxo pipefail
        
        echo "Setting up script..."
        cp ~{reporting_script} ./stitchHtmlReport.py
        chmod +x ./stitchHtmlReport.py
        
        echo "Creating localized data directory..."
        mkdir -p localized_data
        
        echo "Localizing QC plots directory..."
        # Copy QC plots directory contents to localized directory
        gsutil -m cp -r "~{QC_plots_dir}"/* localized_data/
        
        echo "Localizing CCS report files..."
        # Copy all CCS report files to localized directory
        ~{sep=' ' CCS_report}
        for ccs_file in ~{sep=' ' CCS_report}; do
            cp "$ccs_file" localized_data/
        done
        
        echo "Localizing saturation plot..."
        # Copy saturation plot to localized directory
        cp ~{saturation_plot_png} localized_data/
        
        echo "Contents of localized_data directory:"
        ls -la localized_data/
        
        echo "Starting report generation..."
        echo "Sample ID: ~{sample_id}"
        echo "Prefixes: ~{prefixes}"
        echo "Barcodes: ~{barcodes}"
        echo "Localized data directory: localized_data"
        
        # Run the Python script with localized_data as the data_dir
        python stitchHtmlReport.py \
            --prefixes "~{prefixes}" \
            --barcodes "~{barcodes}" \
            --data-dir "localized_data" \
            --output "~{final_output}" \
            --sample-id "~{sample_id}" \
            --base-inserts-run "~{base_inserts_run}" \
            --base-inserts-prefix "~{base_inserts_prefix}"
        
        echo "Report generation completed!"
        echo "Output file: ~{final_output}"
        
        # Verify the output file was created
        if [[ -f "~{final_output}" ]]; then
            echo "Report file successfully created"
            ls -lh "~{final_output}"
        else
            echo "ERROR: Report file was not created"
            exit 1
        fi
    >>>
    
    output {
        File report = final_output
    }
    
    runtime {
        docker: select_first([docker_image, default_docker])
        memory: select_first([mem_gb, default_mem_gb]) + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        cpu: select_first([cpu, default_cpu])
        preemptible: select_first([preemptible_attempts, 0])
    }
}

workflow GenerateReportWorkflow {
    meta {
        description: "Workflow to generate HTML report consolidating QC reports across flowcells for a Kinnex single cell sample"
    }
    
    input {
        String prefixes
        String barcodes
        String QC_plots_dir 
        String sample_id
        File base_inserts_run
        File base_inserts_prefix
        Array[File] CCS_report 
        File saturation_plot_png  
        File reporting_script
        
        String? output_filename
        String? docker_image
        Int? mem_gb
        Int? disk_space_gb
        Int? cpu
        Int? preemptible_attempts
        Int? boot_disk_size_gb
    }
    
    call generateReport {
        input:
            prefixes = prefixes,
            barcodes = barcodes,
            QC_plots_dir = QC_plots_dir,  
            sample_id = sample_id,
            base_inserts_run = base_inserts_run,
            base_inserts_prefix = base_inserts_prefix,
            CCS_report = CCS_report,  
            saturation_plot_png = saturation_plot_png,  
            reporting_script = reporting_script,
            output_filename = output_filename,
            docker_image = docker_image,
            mem_gb = mem_gb,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            preemptible_attempts = preemptible_attempts,
            boot_disk_size_gb = boot_disk_size_gb
    }
    
    output {
        File report = generateReport.report
    }
}