version 1.0



workflow mega_pipeline {
    input {
    File ref_fasta
    File DoC_interval_list
    File ref_fasta_index
    File ref_fasta_dict
    File inputCram
    File inputCramIndex
    File master_gene_list
    File GOF_gene_list
    String samplename
    String output_files_path
Int memoryGb
Int minBaseQuality
Int minMappingQuality
Int preemptible



    }

    call depthOfCov {
        input: refFasta=ref_fasta,
            intervalList=DoC_interval_list,
            refFastaIndex=ref_fasta_index,
            refFastaDict=ref_fasta_dict,
            inputCram=inputCram,
            inputCramIndex=inputCramIndex,
            sampleName=samplename,
            geneList=geneList,
            memoryGb=memoryGb,
            minBaseQuality=minBaseQuality,
            minMappingQuality=minMappingQuality,
            preemptible=preemptible
            
    }
    

  # Outputs that will be retained when execution is complete  
  output {
        #DOC:
        File sampleGeneSummary = depthOfCov.sampleGeneSummary
        File sampleSummary = depthOfCov.sampleSummary
        Float sampleMeanCoverage = depthOfCov.sampleMeanCoverage
        File sampleStatistics = depthOfCov.sampleStatistics
        File sampleIntervalSummary = depthOfCov.sampleIntervalSummary
        File sampleIntervalStatistics = depthOfCov.sampleIntervalStatistics
        File sampleCumulativeCoverageProportions = depthOfCov.sampleCumulativeCoverageProportions
        File sampleCumulativeCoverageCounts = depthOfCov.sampleCumulativeCoverageCounts
  }
    
}





task depthOfCov {
    input {
    File inputCram
    File inputCramIndex
    Int minBaseQuality
    Int minMappingQuality
    String sampleName
    File intervalList
    File geneList
    File refFasta
    File refFastaDict
    File refFastaIndex
    Int memoryGb
    Int preemptible
    #parameters
    Int addtional_disk_space_gb = 10
    Int disk_space_gb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    }
    command <<<
        ln -s ~{inputCram} Cramfile.cram
        ln -s ~{inputCramIndex} Cramfile.crai

        gatk --java-options "-Xmx16G" DepthOfCoverage -R ~{refFasta} -O ~{sampleName} --omit-depth-output-at-each-base true -pt sample -gene-list ~{geneList} -I Cramfile.cram -L ~{intervalList} --min-base-quality ~{minBaseQuality} --summary-coverage-threshold 10
        
        sed -n "2p" < "~{sampleName}.sample_summary" | cut -d',' -f3 > sample_mean_coverage.txt

        mv "~{sampleName}.sample_gene_summary" "~{sampleName}.sample_gene_summary.csv"
        mv "~{sampleName}.sample_summary" "~{sampleName}.sample_summary.csv"
        mv "~{sampleName}.sample_interval_summary" "~{sampleName}.sample_interval_summary.csv"
        mv "~{sampleName}.sample_statistics" "~{sampleName}.sample_statistics.csv"
        mv "~{sampleName}.sample_cumulative_coverage_proportions" "~{sampleName}.sample_cumulative_coverage_proportions.csv"
        mv "~{sampleName}.sample_interval_statistics" "~{sampleName}.sample_interval_statistics.csv"


    >>>

    output {
        File sampleGeneSummary = "~{sampleName}.sample_gene_summary.csv"
        File sampleSummary = "~{sampleName}.sample_summary.csv"
        Float sampleMeanCoverage = read_float("sample_mean_coverage.txt")
        File sampleStatistics = "~{sampleName}.sample_statistics.csv"
        File sampleIntervalSummary = "~{sampleName}.sample_interval_summary.csv"
        File sampleIntervalStatistics = "~{sampleName}.sample_interval_statistics.csv"
        File sampleCumulativeCoverageProportions = "~{sampleName}.sample_cumulative_coverage_proportions.csv"
        File sampleCumulativeCoverageCounts = "~{sampleName}.sample_cumulative_coverage_counts"
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: "~{memoryGb} GB"
        cpu: "1"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: "~{preemptible}"
    }
}
