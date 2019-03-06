task bait_bias_filter_task_1 {
    
    File  MAF_file
    String ref_base
    String alt_base
    String name
    Float min_pval1
    

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions



    Float min_pval=select_first([min_pval1, 0.01])


    command {
        set -x
        python /opt/src/bait_bias_filter.py ${MAF_file} ${ref_base} ${alt_base} ${name} ${min_pval}
    }

    output {
        File FilteredMaf = "${name}.BaitBiasfilt_${ref_base}to${alt_base}.maf"
        String FilteredCount = read_string("${name}.BaitBiasfilt_${ref_base}to${alt_base}_ncut.txt")
    }

    runtime {
        docker : "chipstewart/bait_bias_filter_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '3'}"
    }

    meta {
        author : "Chip Stewart/ Benny Zhitomirsky"
        email : "stewart@broadinstitute.org"
    }
}

workflow bait_bias_filter {

    call bait_bias_filter_task_1 
}
