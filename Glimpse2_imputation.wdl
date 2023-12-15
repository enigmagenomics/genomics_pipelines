version 1.0

workflow GlimpsePhaseLigate {
    input {
        String sampleid
        File bam_or_vcf
        File index

        String refid
        Array[Directory]+ refBins
        Array[String]+ chrs
    }

    scatter (i in range(length(chrs))) {
        call glimpsePhaseLigate {
            input:
              sampleid = sampleid,
              bam_or_vcf = bam_or_vcf,
              index = index,
              chr = chrs[i],

              refid = refid,
              refBins = refBins[i],
              chunks = "~{refBins[i]}/chunks.txt"
        }
    }
    
    call concat {
    	input:
        	sampleid = sampleid,
            vcfs = glimpsePhaseLigate.imputedVcf
    }

    output {
    	File imputedVcf = concat.combinedVcf
    }
}


task glimpsePhaseLigate {
    input {
        String sampleid
        File bam_or_vcf
        File index
        String inputFlag = '--bam-file' # or --input-gl
        String chr
        
        String refid
        Directory refBins 
        File chunks

        Int cpus = 4
        Int preemptible = 2
    }
    
    
    String base = sep("_", [sampleid, refid, chr])
    Int diskSpace = ceil(4 * size(bam_or_vcf, "GB")) + 12
    
    String here_target = basename(bam_or_vcf)
    String here_index = basename(index)
    
    command <<<
      mkdir sample_imputed
      ls
      
      # in case these were inappropriately separated
      mv ~{bam_or_vcf} ~{here_target}
      mv ~{index} ~{here_index}
      
      ### PHASE ###
      # loop logic taken from the GLIMPSE2 tutorial
      while IFS="" read -r LINE || [ -n "$LINE" ]
      do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        REGS=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f1)
        REGE=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f2)

        GLIMPSE2_phase \
        ~{inputFlag} ~{here_target} \
        --reference ~{refBins}/~{chr}_~{chr}_${REGS}_${REGE}.bin \
        --output sample_imputed/~{base}_${REGS}_${REGE}.bcf \
        --threads ~{cpus}
      done < ~{chunks}
      
      ls sample_imputed
      
      ### LIGATE ###
      ls -1v sample_imputed/~{base}_*.bcf > sample_chunks_list.txt
      GLIMPSE2_ligate --input sample_chunks_list.txt --output ~{base}_imputed.bcf

    >>>

    output {
    	File imputedVcf = "~{base}_imputed.bcf"
    }

    runtime {
        disks: "local-disk ${diskSpace} HDD"
        memory: "5 GB"
        cpus: cpus
        preemptible: preemptible
        docker: "simrub/glimpse:v2.0.0-27-g0919952_20221207"
    }
}

task concat {
    input {
      String sampleid
      Array[File] vcfs
      
      Int diskSpace
      Int memory = 2
      Int cpus = 6
      Int preemptible = 3
    }
    
    command <<<
      bcftools concat -Ov -o "~{sampleid}.vcf.gz" --write-index ~{sep=" " vcfs} --threads ~{cpus}
    >>>

    output {
        File combinedVcf = "~{sampleid}.vcf.gz"
    }

    runtime {
        disks: "local-disk ~{diskSpace} HDD"
        memory: "~{memory}GB"
        cpu: cpus
        preemptible: preemptible
        docker: "quay.io/biocontainers/bcftools:1.18--h8b25389_0"
    }
}