version 1.0

task generateReport {
    meta {
        description: "Generates an HTML report using a Python script with pre-installed dependencies"
    }
    
    input {
        # Required inputs
        String prefixes
        String barcodes
        String data_dir
        String sample_id
        File base_inserts_run
        File base_inserts_prefix
        File reporting_script  # The generate_report.py script
        
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
    String default_docker = "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/sc_reporting:latest"  # Replace with your custom image
    
    # Determine output filename
    String final_output = select_first([output_filename, "${sample_id}_report.html"])
    
    command <<<
        set -euxo pipefail
        
        echo "Setting up script..."
        cp ~{reporting_script} ./stitchHtmlReport.py
        chmod +x ./stitchHtmlReport.py
        
        echo "Starting report generation..."
        echo "Sample ID: ~{sample_id}"
        echo "Prefixes: ~{prefixes}"
        echo "Barcodes: ~{barcodes}"
        echo "Data directory: ~{data_dir}"
        
        # Run the Python script
        python stitchHtmlReport.py \
            --prefixes "~{prefixes}" \
            --barcodes "~{barcodes}" \
            --data-dir "~{data_dir}" \
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
        description: "Workflow to generate HTML reports using Python script with custom Docker image"
    }
    
    input {
        String prefixes
        String barcodes
        String data_dir
        String sample_id
        File base_inserts_run
        File base_inserts_prefix
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
            data_dir = data_dir,
            sample_id = sample_id,
            base_inserts_run = base_inserts_run,
            base_inserts_prefix = base_inserts_prefix,
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