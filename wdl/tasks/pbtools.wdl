version 1.0

task pbSkerawQC {
    meta {
        description: "Given hifi reads, spilts MAS 8x or 16x array structure using provided adapters"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File hifi_bam
        String? sample_id
        File mas_adapters_fasta
        Int num_threads
        Int arraysize = 8
        String gcs_output_dir

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(hifi_bam, "GiB"))
    Int default_ram = 8
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25

    # Mem is in units of GB
    Int machine_mem = select_first([mem_gb,default_ram])
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    String skera_id = if defined(sample_id) then sample_id else sub(basename(hifi_bam,".bam"),".hifi_reads","")
    command <<<
        set -euxo pipefail
        
        echo "skera split initiated.."
        echo ~{skera_id}
        
        skera split -j ~{num_threads} ~{hifi_bam} ~{mas_adapters_fasta} ~{skera_id}.skera.bam
        echo "Skera split completed!"

        echo "Generating QC plots.."
        gsutil -m cp -r gs://mdl_terra_sandbox/tools/pb_plots/ .

        python ./pb_plots/plot_concat_hist.py \
        --csv ~{skera_id}.skera.read_lengths.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.concat_hist.png

        python ./pb_plots/plot_readlen_hist.py \
        --csv ~{skera_id}.skera.read_lengths.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.readlen_hist.png

        python ./pb_plots/plot_ligation_heatmap.py \
        --csv ~{skera_id}.skera.ligations.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.ligations_heatmap.png

        echo "Copying output to gcs path provided..."
        gsutil -m cp ~{skera_id}.skera.* ~{outdir}skera/
        echo "Copying skera files completed!"

        echo "Copying plots to gcs path QC_plots..."
        gsutil -m cp ~{skera_id}*.png ~{outdir}QC_plots/
        echo "Copying completed!"
    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        File skera_out        = "~{skera_id}.skera.bam"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }

}

task pbLimaBulk {
    meta {
        description: "Given deconcatenated S-reads for Bulk samples, de-multiplexes using provided primers fasta and trims PolyA tails using Isoseq Refine"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File skera_bam
        String? sample_id
        File bulk_barcodes_fasta
        Boolean trimPolyA = true
        Boolean clipAdapters = true
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 16
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    String isoseq_cmd = if trimPolyA then "isoseq refine --require-polya" else "isoseq refine"
    String lima_cmd = if clipAdapters then "lima --isoseq --log-level INFO" else "lima --isoseq --no-clip --log-level INFO"
    String skera_id = if defined(sample_id) then sample_id else basename(skera_bam,".skera.bam")
    command <<<
        set -euxo pipefail
        
        echo "Running lima demux.."
        ~{lima_cmd} -j ~{num_threads} ~{skera_bam} ~{bulk_barcodes_fasta} ~{skera_id}.lima.bam
        echo "Demuxing completed."

        echo "Copying output to gcs path provided..."
        gsutil -m cp ~{skera_id}*lima* ~{outdir}lima/
        echo "Copying lima files completed!"

        echo "Running Refine..."
        for i in `ls ./*_5p--3p.bam`;
        do
         echo `basename $i`
         a=`basename $i | awk -v FS='_5p--3p.bam' '{print $1}' | awk -v FS='.lima.' '{print $1"."$2}'`
        echo $a
         ~{isoseq_cmd} -j ~{num_threads} $i ~{bulk_barcodes_fasta} ./$a.refine.bam
        done
        echo "Refine completed."

        echo "Uploading refined bams..."
        gsutil -m cp ~{skera_id}*refine* ~{outdir}refine/
        echo "Copying extracted FLNC reads completed!"

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String refine_out  = "~{outdir}refine"
        String lima_out    = "~{outdir}lima"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag4"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 4])
    }
}

task bulkMerge {
    meta {
        description: "Given demuxed refined reads for bulk samples, collapse replicates if mergeBams set to true else plot counts"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        String refine_bampath
        String lima_dir
        String? datasetId = "Replicates_merged" 
        File barcode_to_sample
        File bulk_barcodes_fasta
        Boolean mergeBams = false 
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    #Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 16
    Int default_disk_space_gb = 500
    #Int default_disk_space_gb = ceil((2.5 * 1024 * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    String refinedir = sub(sub( refine_bampath + "/", "/+", "/"), "gs:/", "gs://")
    String limadir = sub(sub( lima_dir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        echo "Fetching refined bams to combine replicates..."
        gsutil -m cp ~{refinedir}*refine.bam .
        echo "Fetching lima counts files.."
        gsutil -m cp ~{limadir}*lima.counts .
 
        echo "plot counts and merge"

        gsutil -m cp -r gs://mdl_terra_sandbox/tools/mergeBam/mergeBams.py .
        python mergeBams.py \
            -idmap ~{barcode_to_sample} \
            -bampath . \
            -limacountsdir . \
            -outdir . \
            -mergeReplicates \
            -setTitleSamplePlot ~{datasetId} 

        echo "Uploading merged bams to merge dir..."  
        ls -lhrt  
        gsutil -m cp ./merge/* ~{outdir}merge/
        gsutil cp readcounts_by_sample.png ~{outdir}merge/
        gsutil cp aggregated_lima_counts_by_sample.tsv ~{outdir}merge/
        gsutil cp lima_counts_by_moviename.tsv ~{outdir}merge/

        echo "Completed copying merged outs and QC plot. All done"

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String merge_out        = "~{outdir}merge"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag5"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 4])
    }
}

task pbSingleCell {
    meta {
        description: "Given s-reads, performs all operations from lima and isoseq pipelines."
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File skera_bam
        String sample_id
        File primer_fasta
        File barcodes_list
        String read_design
        Boolean trimPolyA = true
        Boolean clipAdapters = true 
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }

     # Computing required disk size
    Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 32
    Int default_disk_space_gb = ceil((input_files_size_gb * 5) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    String isoseq_cmd = if trimPolyA then "isoseq refine --require-polya" else "isoseq refine"
    String lima_cmd = if clipAdapters then "lima --isoseq --log-level INFO" else "lima --isoseq --no-clip --log-level INFO"
    command <<<
        set -euxo pipefail

        gsutil -m cp ~{skera_bam} .
        echo "Running lima demux.."
        ~{lima_cmd} -j ~{num_threads} *.skera.bam ~{primer_fasta} ~{sample_id}.lima.bam
        echo "Demuxing completed."

        echo "Copying output to gcs path provided..."
        gsutil -m cp ~{sample_id}*lima* ~{outdir}lima/
        echo "Copying lima files completed!"

        echo "Running Tag, Refine and Correct..."
        for i in `ls ./*.5p--3p.bam`;
        do
           echo `basename $i`
           a=`basename $i | awk -v FS='.5p--3p.bam' '{print $1}' | awk -v FS='.lima.' '{print $1}'`
           echo $a
           echo "tagging.."
           isoseq tag --design ~{read_design} -j ~{num_threads} $i ./$a.tagged.bam 
           echo "refining.."
           ~{isoseq_cmd} -j ~{num_threads} ./$a.tagged.bam ~{primer_fasta} ./$a.refine.bam
           echo "correcting.."
           isoseq correct --barcodes ~{barcodes_list} -j ~{num_threads} ./$a.refine.bam ./$a.corrected.bam 
           echo "Done for this id!"
        done
        echo "All completed."

        echo "Uploading tagged bams..."
        gsutil -m cp *.tagged* ~{outdir}tag/
        echo "Copying extracted tagged reads completed!" 

        echo "Uploading refined bams..."
        gsutil -m cp *.refine* ~{outdir}refine/
        echo "Copying refined reads completed!"

        echo "Uploading corrected bams..."
        gsutil -m cp *.corrected* ~{outdir}correct/
        echo "Copying corrected reads completed!"  

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String corrected_reads  =  "~{outdir}correct"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag5"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }

}

task pbGroupdedup {
    meta {
        description: "Given merged, sorted bams, performs isoseq groupdededup pipelines."
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File input_bam
        String sample_id
        Boolean keep_non_real_cells = true
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }

     # Computing required disk size
    Float input_files_size_gb = 2.5*(size(input_bam, "GiB"))
    Int default_ram = 32
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = select_first([mem_gb,default_ram]) 
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    String isoseq_cmd = if keep_non_real_cells then "isoseq groupdedup --keep-non-real-cells" else "isoseq groupdedup"

    command <<<
        set -euxo pipefail

        echo "Running groupdedup.."
        ~{isoseq_cmd} -j ~{num_threads} ~{input_bam} ~{sample_id}.dedup.bam
        echo "Deduping completed."

        echo "Uploading deduped bams..."
        gsutil -m cp *.dedup* ~{outdir}groupdedup/
        echo "Copying extracted tagged reads completed!" 


    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String dedup_out = "~{outdir}groupdedup"        
        File deduped_bam  =  "~{sample_id}.dedup.bam"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag5"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }

}
