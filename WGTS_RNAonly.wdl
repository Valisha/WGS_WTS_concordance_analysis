version 1.0

## Workflow definition
workflow dragenSomatic {
    input {
      File ht_reference
      File gencode_v25_gtf
      File rna_csv
      String rna_id
    }

    call dragenRna {
      input:
        ht_reference = ht_reference,
        gencode_v25_gtf = gencode_v25_gtf,
        rna_csv = rna_csv,
        rna_id = rna_id,
    }


    output {
        File rna_bam = dragenRna.rna_bam
        File rna_bam_index = dragenRna.rna_bam_index
        File rna_chimeric_junction = dragenRna.rna_chimeric_junction
        File rna_fastqc = dragenRna.rna_fastqc
        File rna_fusion_features = dragenRna.rna_fusion_features
        File rna_fusion_filter = dragenRna.rna_fusion_filter
        File rna_fusion_final = dragenRna.rna_fusion_final
        File rna_fusion_preliminary = dragenRna.rna_fusion_preliminary
        File rna_insert_stats = dragenRna.rna_insert_stats
        File rna_mapping_metrics = dragenRna.rna_mapping_metrics
        File rna_quant_eq_classes = dragenRna.rna_quant_eq_classes
        File rna_quant_genes = dragenRna.rna_quant_genes
        File rna_quant_metrics = dragenRna.rna_quant_metrics
        File rna_quant_sf = dragenRna.rna_quant_sf
        File rna_quant_transcript_coverage = dragenRna.rna_quant_transcript_coverage
        File rna_quant_transcript_frag_len = dragenRna.rna_quant_transcript_frag_len
        File rna_replay_json = dragenRna.rna_replay_json
        File rna_sjdb_annotations_bin = dragenRna.rna_sjdb_annotations_bin
        File rna_sjdb_annotations_out = dragenRna.rna_sjdb_annotations_out
        File rna_sj_out = dragenRna.rna_sj_out
        File rna_sj_saturation = dragenRna.rna_sj_saturation
        File rna_time_metrics = dragenRna.rna_time_metrics
        File rna_trimmer_metrics = dragenRna.rna_trimmer_metrics
        File rna_unfiltered_sj = dragenRna.rna_unfiltered_sj
        File rna_wgs_contig_mean_cov = dragenRna.rna_wgs_contig_mean_cov
        File rna_wggs_cov_metrics = dragenRna.rna_wggs_cov_metrics
        File rna_wgs_fine_hist = dragenRna.rna_wgs_fine_hist
        File rna_wgs_hist = dragenRna.rna_wgs_hist
        File rna_wgs_overall_mean_cov = dragenRna.rna_wgs_overall_mean_cov
      }
}

## Task definitions
task dragenRna {
      input {
        File ht_reference
        File gencode_v25_gtf
        File rna_csv
        String rna_id
      }
    command <<<
        tar xvf ~{ht_reference}
        mv hash-hg38 /ephemeral/hash-hg38
        rm ~{ht_reference}
        mkdir -p /ephemeral/rna
        /opt/edico/bin/dragen \
          --lic-server=https://bdM2KvfvNyA=:uNaCli7EgMpjfMQdGjzkQT8R1lAcBlZR@license.edicogenome.com \
          --ref-dir=/ephemeral/hash-hg38 \
          --fastq-list=~{rna_csv} \
          --fastq-list-sample-id=~{rna_id} \
          --annotation-file=~{gencode_v25_gtf} \
          --intermediate-results-dir=/staging \
          --output-directory=/ephemeral/rna \
          --output-file-prefix=~{rna_id} \
          --enable-rna=true \
          --enable-rna-gene-fusion=true \
          --enable-rna-quantification=true
    >>>
    runtime {
        memory: "224 GB"
        cpu: 16
        disks: "/ephemeral 860 NVMe, /staging 2048 HDD"
        queueArn: "arn:aws:batch:us-east-1:876086242506:job-queue/somatic-dragen"
        docker: "876086242506.dkr.ecr.us-east-1.amazonaws.com/awsclinicaldragen:3.8.5"
    }
    output {
      File rna_bam = "/ephemeral/rna/~{rna_id}.bam"
      File rna_bam_index = "/ephemeral/rna/~{rna_id}.bam.bai"
      File rna_chimeric_junction = "/ephemeral/rna/~{rna_id}.Chimeric.out.junction"
      File rna_fastqc = "/ephemeral/rna/~{rna_id}.fastqc_metrics.csv"
      File rna_fusion_features = "/ephemeral/rna/~{rna_id}.fusion_candidates.features.csv"
      File rna_fusion_filter = "/ephemeral/rna/~{rna_id}.fusion_candidates.filter_info"
      File rna_fusion_final = "/ephemeral/rna/~{rna_id}.fusion_candidates.final"
      File rna_fusion_preliminary = "/ephemeral/rna/~{rna_id}.fusion_candidates.preliminary"
      File rna_insert_stats = "/ephemeral/rna/~{rna_id}.insert-stats.tab"
      File rna_mapping_metrics = "/ephemeral/rna/~{rna_id}.mapping_metrics.csv"
      File rna_quant_eq_classes = "/ephemeral/rna/~{rna_id}.quant.eq_classes.txt"
      File rna_quant_genes = "/ephemeral/rna/~{rna_id}.quant.genes.sf"
      File rna_quant_metrics = "/ephemeral/rna/~{rna_id}.quant.metrics.csv"
      File rna_quant_sf = "/ephemeral/rna/~{rna_id}.quant.sf"
      File rna_quant_transcript_coverage = "/ephemeral/rna/~{rna_id}.quant.transcript_coverage.txt"
      File rna_quant_transcript_frag_len = "/ephemeral/rna/~{rna_id}.quant.transcript_fragment_lengths.txt"
      File rna_replay_json = "/ephemeral/rna/~{rna_id}-replay.json"
      File rna_sjdb_annotations_bin = "/ephemeral/rna/~{rna_id}.sjdb.annotations.bin"
      File rna_sjdb_annotations_out = "/ephemeral/rna/~{rna_id}.sjdb.annotations.out.tab"
      File rna_sj_out = "/ephemeral/rna/~{rna_id}.SJ.out.tab"
      File rna_sj_saturation = "/ephemeral/rna/~{rna_id}.SJ.saturation.txt"
      File rna_time_metrics = "/ephemeral/rna/~{rna_id}.time_metrics.csv"
      File rna_trimmer_metrics = "/ephemeral/rna/~{rna_id}.trimmer_metrics.csv"
      File rna_unfiltered_sj = "/ephemeral/rna/~{rna_id}.unfiltered.SJ.out.tab"
      File rna_wgs_contig_mean_cov = "/ephemeral/rna/~{rna_id}.wgs_contig_mean_cov.csv"
      File rna_wggs_cov_metrics = "/ephemeral/rna/~{rna_id}.wgs_coverage_metrics.csv"
      File rna_wgs_fine_hist = "/ephemeral/rna/~{rna_id}.wgs_fine_hist.csv"
      File rna_wgs_hist = "/ephemeral/rna/~{rna_id}.wgs_hist.csv"
      File rna_wgs_overall_mean_cov = "/ephemeral/rna/~{rna_id}.wgs_overall_mean_cov.csv"
    }
}
