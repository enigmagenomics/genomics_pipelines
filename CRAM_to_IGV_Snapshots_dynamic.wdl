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
    String TestType
File    geneList
Int runtime_disk
File    WGS_interval_list
File    WES_interval_list
Int diskSpace
Int resource_log_interval
Int runtime_cpus
String  runtime_docker
Int runtime_preemptible
Int memoryGb
Int minBaseQuality
Int minMappingQuality
Int preemptible
Int GATK_diskGb
Int GATK_memoryGb
Int memoryGb
Int RAM
Int HDD
Int boot_disk_gb
File    chain
File    Clinvar_genes
Int cpu_cores
File    GCD_genes
File    OMIM_genes
Int output_disk_gb
Int ram_gb
File    Rscript_file
Int bootDiskSizeGb_VEP
Int cpu_VEP
Int diskGb_VEP
Int fork
Int memoryGb_VEP
Int nearestGeneDistance


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
    
    
    
    call ChooseBed {
    input:
         TestType=TestType, 
         WGS_interval_list=WGS_interval_list,
        WES_interval_list=WES_interval_list

  }
    

    
    call interval_list_to_bed {
            input:
                interval_list=ChooseBed.model_interval_list
    }

  
    call deep_variant {
        input:
            sample=samplename,
            capture_bed=interval_list_to_bed.bed,
            Cram=inputCram,
            crai=inputCramIndex,
            reference_fasta=ref_fasta,
            reference_fasta_fai=ref_fasta_index,
            model_type=TestType,
            resource_log_interval=resource_log_interval,
            runtime_cpus=runtime_cpus,
            runtime_docker=runtime_docker,
            runtime_preemptible=runtime_preemptible
    }

    call bgzip {
        input:
            sample=samplename,
            uncompressed_vcf=deep_variant.vcf,
            runtime_disk=runtime_disk
    }

    call variantcount_vcf {
        input:
            vcf = bgzip.filtered_vcf,
            sampleId=samplename,
            HDD=HDD,
            RAM=RAM

    }   

    
    call UnZip { 
        input:
            vcfFileGz = bgzip.filtered_vcf,
            sampleId = samplename,
            RAM=RAM,
            HDD=HDD
    }
    
    call VTRecal {
        input:

            vcfFile=UnZip.vcfFile,
            refFasta=ref_fasta,
            refFastaDict=ref_fasta_dict,
            refFastaIdx=ref_fasta_index,
            sampleId=samplename, 
            RAM=RAM,
            HDD=HDD
    }
    
    call GATKVariantsToTable {
        input:
            normalizedvcfFileGz=VTRecal.normalizedVCF,
            refFasta=ref_fasta,
            refFastaFai=ref_fasta_index,
            refFastaDict=ref_fasta_dict,
            samplesetId=samplename,
            GATK_diskGb=GATK_diskGb,
            GATK_memoryGb=GATK_memoryGb
           
    }
    call vep_task {
        input:
       
            refFasta=ref_fasta,
            refFastaFai=ref_fasta_index,
            refFastaDict=ref_fasta_dict,
            samplesetId=samplename,
            normalizedvcfFileGz=VTRecal.normalizedVCF,
            bootDiskSizeGb_VEP=bootDiskSizeGb_VEP,
            cpu_VEP=cpu_VEP,
            diskGb_VEP=diskGb_VEP,
            fork=fork,
            memoryGb_VEP=memoryGb_VEP,
            nearestGeneDistance=nearestGeneDistance
    }
    
    call combineOutputFiles {
        input:
            samplesetId=samplename,
            vepOutputFile=vep_task.VEP_Output,
            gatkOutputFile=GATKVariantsToTable.GATK_output,
            diskSpace=diskSpace
            
    }

    call VariantFilter{
   input: 
     input_vcf=combineOutputFiles.vepannotated_vcf,
     master_gene_list=master_gene_list,
     GOF_gene_list=GOF_gene_list,
     samplename=samplename,
     boot_disk_gb=boot_disk_gb,
     chain=chain,
     Clinvar_genes=Clinvar_genes,
     cpu_cores=cpu_cores,
     GCD_genes=GCD_genes,
     OMIM_genes=OMIM_genes,
     output_disk_gb=output_disk_gb,
     ram_gb=ram_gb,
     Rscript_file=Rscript_file
  }

    call IGV_Snapshots{
    input: 
    inputCram=inputCram,
    inputCramIndex=inputCramIndex,
    path_var_HQ_IGV_bed=VariantFilter.path_var_HQ_IGV_bed,
    path_var_HQ_non_clinical_IGV_bed=VariantFilter.path_var_HQ_non_clinical_IGV_bed,
    path_var_LQ_IGV_bed=VariantFilter.path_var_LQ_IGV_bed,
    sampleID=samplename,
    memoryGb=memoryGb
  }

    call data_transfer{
            input:
    output_files_path=output_files_path,
    total_variant_count=VariantFilter.total_variant_count,
    DOC_sampleSummary=depthOfCov.sampleSummary,
    DOC_sampleCumulativeCoverageProportions=depthOfCov.sampleCumulativeCoverageProportions,
    DOC_sampleGeneSummary=depthOfCov.sampleGeneSummary,
    filtered_vcf=bgzip.filtered_vcf,
    DV_stats_report=deep_variant.stats_report,
    path_var_HQ=VariantFilter.path_var_HQ,
    path_var_HQ_non_clinical=VariantFilter.path_var_HQ_non_clinical,
    path_var_LQ=VariantFilter.path_var_LQ,
    variants_HQ_IGV_snapshots=IGV_Snapshots.variants_HQ_IGV_snapshots,
    variants_HQ_non_clinical_IGV_snapshots=IGV_Snapshots.variants_HQ_non_clinical_IGV_snapshots,
    variants_LQ_IGV_snapshots=IGV_Snapshots.variants_LQ_IGV_snapshots,
    all_variants_cancer_genes=VariantFilter.all_variants_cancer_genes, 
    all_variants_carrier_genes=VariantFilter.all_variants_carrier_genes,
    sampleID=samplename
    }

  # Outputs that will be retained when execution is complete  
  output {
        #DOC:
        File sampleGeneSummary = depthOfCov.sampleGeneSummary
        File sampleSummary = depthOfCov.sampleSummary
        Float sampleMeanCoverage = depthOfCov.sampleMeanCoverage
        Float sample10XCoverage = depthOfCov.sample10XCoverage
        File sampleStatistics = depthOfCov.sampleStatistics
        File sampleIntervalSummary = depthOfCov.sampleIntervalSummary
        File sampleIntervalStatistics = depthOfCov.sampleIntervalStatistics
        File sampleCumulativeCoverageProportions = depthOfCov.sampleCumulativeCoverageProportions
        File sampleCumulativeCoverageCounts = depthOfCov.sampleCumulativeCoverageCounts
        #DV:
        File gvcf = deep_variant.gvcf
        File resource_log = deep_variant.resource_log
        File stats_report = deep_variant.stats_report
        File filtered_vcf = bgzip.filtered_vcf
        File filtered_vcf_index = bgzip.filtered_vcf_index
        String variantcount = variantcount_vcf.variantcount
        File normalizedVCF= VTRecal.normalizedVCF
        File normalizedVCF_index= VTRecal.normalizedVCF_index
        #VEP:
        File vepannotated_vcf= combineOutputFiles.vepannotated_vcf
        #variant_filtering:
        File path_var_HQ= VariantFilter.path_var_HQ
        File path_var_HQ_non_clinical= VariantFilter.path_var_HQ_non_clinical        
        File path_var_LQ= VariantFilter.path_var_LQ
        #IGV:
        File variants_HQ_IGV_snapshots = IGV_Snapshots.variants_HQ_IGV_snapshots
        File variants_HQ_non_clinical_IGV_snapshots = IGV_Snapshots.variants_HQ_non_clinical_IGV_snapshots
        File variants_LQ_IGV_snapshots = IGV_Snapshots.variants_LQ_IGV_snapshots
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
        sed -n "2p" < "~{sampleName}.sample_summary" | cut -d',' -f7 > sample_10X_coverage.txt

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
        Float sample10XCoverage = read_float("sample_10X_coverage.txt")
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


task ChooseBed {
    input {
  File WES_interval_list
  File WGS_interval_list
  String TestType
}


   command <<<

# Check if the input variable is "WES"
if [ ~{TestType} = "WES" ]; then
    bed=~{WES_interval_list}
else
    bed=~{WGS_interval_list}
fi

    mv "$bed" "model_interval_list.bed"



    >>>


  output {
    File model_interval_list = "model_interval_list.bed"
  }
  
    runtime {
        memory: '2 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}


task interval_list_to_bed {
    input {
    File interval_list
    String bed_path = sub(basename(interval_list), 'interval_list', 'bed')
}
    command <<<
    set -xeuo pipefail

    # interval lists have headers that need to be removed and are 1-indexed
    # see also https://www.biostars.org/p/84686/
    grep -v '^@' ~{interval_list} \
    | awk -v OFS='\t' '{print $1, $2 - 1, $3}' \
    | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge \
    > ~{bed_path}
    >>>

    output {
        File bed = '~{bed_path}'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk 1 HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bedtools:2.28.0--hdf88d34_0'
    }
}




task deep_variant {
    input {
    String sample
    File Cram
    File crai
    String model_type 
    File reference_fasta
    File reference_fasta_fai
    File capture_bed

    Int runtime_cpus
    String runtime_docker
    Int runtime_memory = ceil(1.1 * runtime_cpus)
    Int addtional_disk_space_gb = 100
    Int disk_space_gb = ceil(size(Cram, "GB")  * 4 ) + addtional_disk_space_gb
    Int runtime_preemptible
    Int resource_log_interval
}
    command <<<
    # log resource usage for debugging purposes
    function runtimeInfo() {
        echo [$(date)]
        echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true;
        do runtimeInfo >> resource.log;
        sleep ~{resource_log_interval};
    done &
    lscpu

    set -xeuo pipefail

    # make symbolic links to ensure Cram and index are in expected structure even after localization
    ln -s ~{crai} reads.crai
    ln -s ~{Cram} reads.cram

    # make symbolic links to ensure reference and index are in expected structure even after localization
    ln -s ~{reference_fasta} reference.fa
    ln -s ~{reference_fasta_fai} reference.fa.fai

    mkdir deepvariant_tmp

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{model_type} \
        --ref=reference.fa \
        --reads=reads.cram \
        --regions=~{capture_bed} \
        --intermediate_results_dir=deepvariant_tmp \
        --output_vcf=~{sample}.vcf \
        --output_gvcf=~{sample}.gvcf.gz \
        --num_shards=~{runtime_cpus}
    >>>

    output {
        File vcf = '~{sample}.vcf'
        File gvcf = '~{sample}.gvcf.gz'
        File resource_log = 'resource.log'
        File stats_report = '~{sample}.visual_report.html'
    }
    
    runtime {
        memory: '~{runtime_memory} GB'
        cpu: '~{runtime_cpus}'
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: '~{runtime_preemptible}'
        docker: '~{runtime_docker}'
    }
}

task bgzip {
    input {
    String sample
    File uncompressed_vcf

    Int runtime_disk
}
    command <<<
    set -xeuo pipefail

    bcftools view -Oz -o ~{sample}.vcf.gz ~{uncompressed_vcf}
    bcftools index --tbi ~{sample}.vcf.gz

    # create version of VCF with only PASSing variants
    bcftools view -Oz -o ~{sample}.filtered_callset.vcf.gz -f PASS ~{uncompressed_vcf}
    bcftools index --tbi ~{sample}.filtered_callset.vcf.gz
    >>>

    output {
        File vcf = '~{sample}.vcf.gz'
        File vcf_index = '~{sample}.vcf.gz.tbi'
        File filtered_vcf = '~{sample}.filtered_callset.vcf.gz'
        File filtered_vcf_index = '~{sample}.filtered_callset.vcf.gz.tbi'
    }

    runtime {
        memory: '1 GB'
        disks: 'local-disk ~{runtime_disk} HDD'
        preemptible: 3
        docker: 'quay.io/biocontainers/bcftools:1.9--ha228f0b_3'
    }
}


task variantcount_vcf {
    input {
    File vcf
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    
        bcftools query -f 'pos=%POS\n' ~{vcf} -o temp.txt
        cat temp.txt | wc -l > ~{sampleId}.variantcount.txt
        
   >>>

    output {
        String variantcount = read_string("~{sampleId}.variantcount.txt")
        File variantcount_file= '~{sampleId}.variantcount.txt'
    }
    
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: RAM + " GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: 3
    }
}


task UnZip {
    input {
    File vcfFileGz
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    # Decompress bgzipped merged VCF file
    echo "bgzip -d -c ~{vcfFileGz} > vcfFile.vcf"
    bgzip -d -c ~{vcfFileGz} > ~{sampleId}.vcf

   >>>

    output {
        File vcfFile="~{sampleId}.vcf"
    }
    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + "GB"
        disks: "local-disk" + " " + HDD + " " + "HDD"
        preemptible: "3"
    }
}


task VTRecal {
    input {
    File vcfFile 
    File refFasta
    File refFastaIdx
    File refFastaDict
    String sampleId
    Int RAM
    Int HDD
}
    command <<<
    # VCF:
        echo "########### decompose VCF"
        /software/vt/./vt decompose -s \
        -o ~{sampleId}.vt1.vcf \
        ~{vcfFile}

        echo "########### normalize VCF using ch38 genome build"
        /software/vt/./vt normalize \
        -r ~{refFasta} \
        -o ~{sampleId}.vt2.vcf \
        ~{sampleId}.vt1.vcf
        
        echo "########### normalizing the spanning alleles (*):"
        sed 's/*/-/g' ~{sampleId}.vt2.vcf > ~{sampleId}.vt2_normalized_spanning_alleles.vcf
        bgzip ~{sampleId}.vt2_normalized_spanning_alleles.vcf
        
        echo "########### creating an index for vcf.gz:"
        tabix -p vcf ~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz 


   >>>

    output {
        File normalizedVCF="~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz"
        File normalizedVCF_index="~{sampleId}.vt2_normalized_spanning_alleles.vcf.gz.tbi"

    }

    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + " GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: "3"
    }
}


task vep_task {
    input {
    File normalizedvcfFileGz
    String samplesetId
    # Customizations
    Int nearestGeneDistance
    # Optimizations
    Int fork

    File refFasta
    File refFastaDict
    File refFastaFai
    Int bootDiskSizeGb_VEP
    Int cpu_VEP
    Int diskGb_VEP
    Int memoryGb_VEP
    
    # Cache files
    File speciesCacheTarGzFile="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/homo_sapiens_merged_vep_109_GRCh38.tar.gz"
    String speciesCacheLabel="homo_sapiens_merged"
    String speciesCacheParameter="--merged"
    
    #dbnSFP
    File dbNSFPData="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFPv4.2a_modified_header.gz"
    File dbNSFPDataTbi="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFPv4.2a_modified_header.gz.tbi"
    File dbNSFPPlugin="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFP.pm"
    File dbNSFPReadme="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/dbNSFPv4.2a/dbNSFP4.1a.readme.txt"
    
    #   MCAP
    #File mcap="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/mcap_v1.4_modified.vcf.gz"
    #File mcap_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/mcap_v1.4_modified.vcf.gz.tbi"

    #   S-CAP
    File scap="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/scap_COMBINED_v1.0_modified.vcf.gz"
    File scap_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/scap_COMBINED_v1.0_modified.vcf.gz.tbi"
    
    #gnomAD genome
    #File gnomAD_genome="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
    #File gnomAD_genome_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch37/annotation_files/VEP_files/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi"
    
   
    
     #   CLINVAR
    File clinvar="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/clinvar_20231217.vcf.gz"
    File clinvar_index="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/clinvar_20231217.vcf.gz.tbi"
}
    
    
    command <<<
        # Prepare directories for data
        echo "mkdir /opt/vep/.vep/~{speciesCacheLabel}/"
        mkdir /opt
        mkdir /opt/vep
        mkdir /opt/vep/.vep
        mkdir /opt/vep/.vep/~{speciesCacheLabel}/
        mkdir /opt/vep/.vep/Plugins
        mkdir /opt/vep/.vep/Plugins/data
        
        # Put Plugins in correct folder
        echo "symbolic links..."
       
      
        # dbNSFP
        ln -s ~{dbNSFPPlugin} /opt/vep/.vep/Plugins/dbNSFP.pm
        ln -s ~{dbNSFPData} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz
        ln -s ~{dbNSFPDataTbi} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz.tbi
        ln -s ~{dbNSFPReadme} /opt/vep/.vep/Plugins/data/dbNSFP4.1a.readme.txt
      
        
        # Uncompress the species cache to the .vep directory
        echo "# Uncompress the species cache to the .vep directory"
        echo "tar xzf ~{speciesCacheTarGzFile} -C ~"
        tar xzf ~{speciesCacheTarGzFile} -C .
        
        echo "ls -lh"
        ls -lh
        
        echo "mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}/"
        mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}
        echo "ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*"
        ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*
    
        # log progress to make sure that the VEP output is being generated
        set -xeuo pipefail
        function runtimeInfo() {            
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo;
            sleep 30;
        done &
   
        
        
        # Run VEP
        echo "running VEP..."
        vep -v -i ~{normalizedvcfFileGz} -o ~{samplesetId}.vep.txt \
        --tab \
        --offline --cache ~{speciesCacheParameter} --dir /opt/vep/.vep --fasta ~{refFasta} \
        --force_overwrite --stats_text --symbol --everything \
        --regulatory --distance ~{nearestGeneDistance}  \
        --total_length --numbers --domains --pick --variant_class --hgvs --hgvsg --ccds  --fork ~{fork} \
        --custom ~{scap},scap_v1.0,vcf,exact,0,Allele_region,rawscore,sensscore,rawscore_dom,sensscore_dom,rawscore_rec,senscore_rec \
        --plugin dbNSFP,/opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz,hg19_chr,hg19_pos,FATHMM_score,FATHMM_pred,PROVEAN_score,MetaSVM_score,MetaLR_score,MetaLR_pred,MetaRNN_score,MetaRNN_pred,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,MVP_score,Aloft_pred,LINSIGHT,CADD_raw,GenoCanyon_score,integrated_fitCons_score,Interpro_domain,gnomAD_genomes_MID_AC,gnomAD_genomes_MID_AN,gnomAD_genomes_MID_AF,gnomAD_genomes_MID_nhomalt \
        --custom ~{clinvar},ClinVar_updated_2023Feb,vcf,exact,0,ID,ALLELEID,CLNDN,CLNDISDB,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID 
    
        
        
        echo "ls -lh"
        ls -lh
        
        echo "ls -lh opt/vep/.vep/*"
        ls -lh /opt/vep/.vep/*
        
        echo "Number of VEP variants (grep -v # ~{samplesetId}.vep.txt | wc -l):"
        grep -v "#" ~{samplesetId}.vep.txt | wc -l
        
        gzip ~{samplesetId}.vep.txt
    >>>
    
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.2"    
        bootDiskSizeGb : "~{bootDiskSizeGb_VEP}"
        preemptible    : 0
        cpu            : "~{cpu_VEP}"
        disks          : "local-disk ~{diskGb_VEP} SSD"
        memory         : "~{memoryGb_VEP} GB"
    }

    output {        
        File VEP_Output="~{samplesetId}.vep.txt.gz"
        File VEP_Summary="~{samplesetId}.vep.txt_summary.txt"
    }
}


task GATKVariantsToTable {
input {   
    File normalizedvcfFileGz
    File refFasta
    File refFastaDict
    File refFastaFai
    String samplesetId
    Int GATK_diskGb
    Int GATK_memoryGb
 }   
    command <<<      
      echo "ls -lh"
      ls -lh
      ls ~{normalizedvcfFileGz}
      
      mv ~{normalizedvcfFileGz} vcfFile.vcf.gz
      
      echo "ls -lh"
      ls -lh
      
      echo "bgzip decompressing vt recal VCF file"
      bgzip --decompress vcfFile.vcf.gz
    
      echo "ls -lh"
      ls -lh 
      
      echo "ls -lh vcfFile.vcf"
      ls -lh vcfFile.vcf
      
      echo "########### Using GATK to extract variants into a table format (GRCh38)"
      java -jar /usr/GenomeAnalysisTK.jar -R ~{refFasta} -T VariantsToTable \
      -V vcfFile.vcf -o ~{samplesetId}.vt2_GATK_annotations.txt \
      -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -GF AD -GF DP  -GF GQ  -GF VAF -GF PL -GF GT --allowMissingData --showFiltered 

      echo '########### Done extracting relevant fields using GATK'

      # count the number of GATK variants:
      echo "########### number of GATK variants: "
      cat ~{samplesetId}.vt2_GATK_annotations.txt | wc -l
      
      gzip ~{samplesetId}.vt2_GATK_annotations.txt
      
      echo "ls -lh"
      ls -lh
    >>>
    
    output {
        File GATK_output="~{samplesetId}.vt2_GATK_annotations.txt.gz"
    }
    
    runtime {
        docker: "vanallenlab/gatk3.7_with_htslib1.9:1.0"
        memory: "~{GATK_memoryGb} GB"
        cpu: "1"
        disks: "local-disk ~{GATK_diskGb} HDD"
        preemptible: 0
    }
}


task combineOutputFiles {
    input {
    File vepOutputFile
    File gatkOutputFile
    String samplesetId
    Int diskSpace
}

    command <<<
      cp ~{vepOutputFile} vepOutputFile.txt.gz
      cp ~{gatkOutputFile} gatkOutputFile.txt.gz
    
      gunzip vepOutputFile.txt.gz
      gunzip gatkOutputFile.txt.gz
    
      # remove the '#' from the first line before parsing:
      echo "########### removing the # sign from the first line of the VEP output file"
      sed  "s/#Uploaded_variation/Uploaded_variation/g" vepOutputFile.txt > ~{samplesetId}_vt2_VEP_temp2.txt

      # remove the excess header part:
      grep "#" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_VEP_annotation_list.txt
      grep -v "##" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_vt2_VEP.txt
      #tail -n +124 ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_vt2_VEP.txt

      # count the number of VEP variants:
      echo "########### number of VEP variants"
      cat ~{samplesetId}_vt2_VEP.txt | wc -l
      
      echo "########### Combining VEP  output files"
      paste ~{samplesetId}_vt2_VEP.txt gatkOutputFile.txt | bgzip > ~{samplesetId}_vt2_VEP_Genotypes.txt.gz
    >>>
    
    output {
        File vepannotated_vcf="~{samplesetId}_vt2_VEP_Genotypes.txt.gz"
        File annotationsList="~{samplesetId}_VEP_annotation_list.txt"
    }
    
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.2"    
        preemptible    : 0
        cpu            : "1"
        disks          : "local-disk ~{diskSpace} HDD"
        memory         : "10 GB"
    }
}


task VariantFilter {
    input {
    File input_vcf
    File master_gene_list
    File GOF_gene_list
    File OMIM_genes
    File Clinvar_genes
    File GCD_genes
    File chain
    File Rscript_file
    Int cpu_cores
    Int ram_gb
    Int boot_disk_gb
    Int output_disk_gb
    String samplename
}
   command <<<
       Rscript ~{Rscript_file} ~{master_gene_list} ~{GOF_gene_list} ~{OMIM_genes} ~{Clinvar_genes} ~{GCD_genes} ~{input_vcf} ~{samplename} ~{chain}  
   >>>
   
   runtime { 
     docker : "lbwang/rocker-genome"
     bootDiskSizeGb: "~{boot_disk_gb}"
     preemtible : 0
     disks: "local-disk ~{output_disk_gb} HDD"
     cpu: "~{cpu_cores}"
     memory: "~{ram_gb}GB"
    }   
    output {
        File path_var_HQ="~{samplename}.pathogenic_variants_all_high_quality.csv"
        File path_var_HQ_non_clinical="~{samplename}.pathogenic_variants_all_high_quality_non_clinical.csv"        
        File path_var_LQ="~{samplename}.pathogenic_variants_all_low_quality.csv"
        File path_var_HQ_IGV_bed="~{samplename}_pathogenic_variants_all_high_quality_IGV.bed"
        File path_var_HQ_non_clinical_IGV_bed="~{samplename}_pathogenic_variants_all_high_quality_non_clinical_IGV.bed"
        File path_var_LQ_IGV_bed="~{samplename}_pathogenic_variants_all_low_quality_IGV.bed"
        File all_variants_cancer_genes="~{samplename}.all_variants_cancer_genes.csv.gz"
        File all_variants_carrier_genes="~{samplename}.all_variants_carrier_genes.csv.gz"
        File total_variant_count="~{samplename}.total_variant_count.csv"
        File variants_PGx="${samplename}.variants_PGx.csv"
       }  
}



task IGV_Snapshots {
    input {
    File inputCram
    File inputCramIndex
    File path_var_HQ_IGV_bed
    File path_var_HQ_non_clinical_IGV_bed
    File path_var_LQ_IGV_bed
    Int memoryGb
    Int addtional_disk_space_gb = 10
    Int disk_space_gb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    Int bootDiskSizeGb = ceil(size(inputCram, "GB")  * 2 ) + addtional_disk_space_gb
    String sampleID
    String genome = "hg38"
}
    command <<<
    ls

## HQ:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_HQ_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_HQ.tar.gz
    mv ~{sampleID}.IGV_Snapshots_HQ.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"




    ## HQ_non_clinical:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_HQ_non_clinical_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz
    mv ~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"



## LQ:
    echo "moving input files to correct directory..."
    
    mv ~{path_var_LQ_IGV_bed} /IGV-snapshot-automator/regions.bed

    cd /IGV-snapshot-automator
    
    mv ~{inputCram} /cromwell_root/~{sampleID}.cram
    mv ~{inputCramIndex} /cromwell_root/~{sampleID}.cram.crai

    /bin/sh -c "if [ ! -d 'IGV_Snapshots' ]; then mkdir IGV_Snapshots; fi"
    
    echo "creating batch file..."
    
    python3 make_IGV_snapshots.py -g ~{genome} /cromwell_root/~{sampleID}.cram
    
    xvfb-run --auto-servernum --server-num=1 igv.sh -g ~{genome} -b IGV_Snapshots/IGV_snapshots.bat

    echo "---- finished running igv snapshot automator ----"
    echo "compressing IGV_Snapshots directory for output..."
    
    tar -zcf IGV_Snapshots.tar.gz IGV_Snapshots
    
    echo "renaming & moving IGV_Snapshots directory to cromwell_root..."
    
    mv IGV_Snapshots.tar.gz ~{sampleID}.IGV_Snapshots_LQ.tar.gz
    mv ~{sampleID}.IGV_Snapshots_LQ.tar.gz /cromwell_root
    
    echo "---- completed running IGV_Snapshots ----"


    >>>

    output {
        File variants_HQ_IGV_snapshots = "~{sampleID}.IGV_Snapshots_HQ.tar.gz"
        File variants_HQ_non_clinical_IGV_snapshots = "~{sampleID}.IGV_Snapshots_HQ_non_clinical.tar.gz"
        File variants_LQ_IGV_snapshots = "~{sampleID}.IGV_Snapshots_LQ.tar.gz"
    }

    runtime {
        docker: "tylerchinskydfci/igv_snapshot:0.1"
        memory: "~{memoryGb} GB"
        cpu: "1"
        disks: "local-disk " + disk_space_gb + " HDD"
    }
}



task data_transfer {
    input {
    File DOC_sampleSummary
    File DOC_sampleCumulativeCoverageProportions
    File DOC_sampleGeneSummary
    File filtered_vcf
    File DV_stats_report
    File path_var_HQ
    File path_var_HQ_non_clinical
    File path_var_LQ
    File variants_HQ_IGV_snapshots
    File variants_HQ_non_clinical_IGV_snapshots
    File variants_LQ_IGV_snapshots
    File total_variant_count
    File all_variants_carrier_genes
    File all_variants_cancer_genes 
    String output_files_path
    String sampleID
}

    command <<<
    
    # copy the files to local disk:
    mv ~{DOC_sampleSummary}  ~/~{sampleID}.sample_summary.csv
    mv ~{DOC_sampleCumulativeCoverageProportions}  ~/~{sampleID}.sample_cumulative_coverage_proportions.csv
    mv ~{DOC_sampleGeneSummary}  ~/~{sampleID}.sample_gene_summary.csv
    mv ~{filtered_vcf}  ~/~{sampleID}.filtered_callset.vcf.gz
    mv ~{path_var_HQ}  ~/~{sampleID}.path_var_HQ.csv
    mv ~{path_var_HQ_non_clinical}  ~/~{sampleID}.path_var_HQ_non_clinical.csv
    mv ~{path_var_LQ}  ~/~{sampleID}.path_var_LQ.csv
    mv ~{variants_HQ_IGV_snapshots}  ~/~{sampleID}.variants_HQ_IGV_snapshots.tar.gz
    mv ~{variants_HQ_non_clinical_IGV_snapshots}  ~/~{sampleID}.variants_HQ_non_clinical_IGV_snapshots.tar.gz
    mv ~{variants_LQ_IGV_snapshots}  ~/~{sampleID}.variants_LQ_IGV_snapshots.tar.gz
    mv ~{total_variant_count}  ~/~{sampleID}.total_variant_count.csv
    mv ~{all_variants_carrier_genes}  ~/~{sampleID}.all_variants_carrier_genes.csv.gz
    mv ~{all_variants_cancer_genes}  ~/~{sampleID}.all_variants_cancer_genes.csv.gz
    mv ~{DV_stats_report}  ~/~{sampleID}.visual_report.html



    # Move file to the file repo:
    gsutil cp ~/~{sampleID}.sample_summary.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.sample_cumulative_coverage_proportions.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.sample_gene_summary.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.filtered_callset.vcf.gz    ~{output_files_path}
    gsutil cp ~/~{sampleID}.path_var_HQ.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.path_var_HQ_non_clinical.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.path_var_LQ.csv    ~{output_files_path}   
    gsutil cp ~/~{sampleID}.variants_HQ_IGV_snapshots.tar.gz    ~{output_files_path}
    gsutil cp ~/~{sampleID}.variants_HQ_non_clinical_IGV_snapshots.tar.gz    ~{output_files_path}
    gsutil cp ~/~{sampleID}.variants_LQ_IGV_snapshots.tar.gz    ~{output_files_path}
    gsutil cp ~/~{sampleID}.total_variant_count.csv    ~{output_files_path}
    gsutil cp ~/~{sampleID}.all_variants_carrier_genes.csv.gz   ~{output_files_path}
    gsutil cp ~/~{sampleID}.all_variants_cancer_genes.csv.gz   ~{output_files_path}
    gsutil cp ~/~{sampleID}.visual_report.html   ~{output_files_path}

   >>>

    output {

    }
    runtime {
        docker: "google/cloud-sdk"
        memory: "1GB"
        disks: 'local-disk 1 HDD' 
        preemptible: "0"
    }
}