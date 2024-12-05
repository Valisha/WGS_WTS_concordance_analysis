version 1.0

## Workflow definition
workflow dragenSomatic {
    input {
      File ht_reference
      File nirvana
      File vc_systematic_noise
      File qc_coverage_region_1
      File qc_coverage_region_2
      File msi_microsatellites
      File gencode_v25_gtf
      File tumor_csv
      File normal_csv
      String normal_id
      String tumor_id
      String sex
    }


    call dragenGermline {
      input:
        ht_reference = ht_reference,
        normal_csv = normal_csv,
        normal_id = normal_id
    }

    call dragenTumorNormal {
      input:
        ht_reference = ht_reference,
        nirvana = nirvana,
        vc_systematic_noise = vc_systematic_noise,
        qc_coverage_region_1 = qc_coverage_region_1,
        qc_coverage_region_2 = qc_coverage_region_2,
        msi_microsatellites = msi_microsatellites,
        tumor_csv = tumor_csv,
        normal_csv = normal_csv,
        normal_id = normal_id,
        tumor_id = tumor_id,
        sex = sex,
        normal_vcf = dragenGermline.germline_filtered_vcf,
        normal_vcf_index = dragenGermline.germline_filtered_vcf_index
    }

    output {
        File germline_baf_bw = dragenGermline.germline_baf_bw
        File germline_bam = dragenGermline.germline_bam
        File germline_bam_index = dragenGermline.germline_bam_index
        File germline_fastqc = dragenGermline.germline_fastqc
        File germline_fragment_len_hist = dragenGermline.germline_fragment_len_hist
        File germline_filtered_vcf_baf_bw = dragenGermline.germline_filtered_vcf_baf_bw
        File germline_filtered_vcf = dragenGermline.germline_filtered_vcf
        File germline_filtered_vcf_index = dragenGermline.germline_filtered_vcf_index
        File germline_insert_stats = dragenGermline.germline_insert_stats
        File germline_gc_metrics = dragenGermline.germline_gc_metrics
        File germline_mapping_metrics = dragenGermline.germline_mapping_metrics
        File germline_pcr_model_0_log = dragenGermline.germline_pcr_model_0_log
        File germline_pcr_model_log = dragenGermline.germline_pcr_model_log
        File germline_ploidy_metrics = dragenGermline.germline_ploidy_metrics
        File germline_ploidy_vcf = dragenGermline.germline_ploidy_vcf
        File germline_ploidy_vcf_index = dragenGermline.germline_ploidy_vcf_index
        File germline_replay_json = dragenGermline.germline_replay_json
        File germline_roh_bed = dragenGermline.germline_roh_bed
        File germline_roh_metrics = dragenGermline.germline_roh_metrics
        File germline_time_metrics = dragenGermline.germline_time_metrics
        File germline_trimmer_metrics = dragenGermline.germline_trimmer_metrics
        File germline_vcf = dragenGermline.germline_vcf
        File germline_vcf_index = dragenGermline.germline_vcf_index
        File germline_vc_metrics = dragenGermline.germline_vc_metrics
        File germline_wgs_contig_mean_cov = dragenGermline.germline_wgs_contig_mean_cov
        File germline_wgs_coverage_metrics = dragenGermline.germline_wgs_coverage_metrics
        File germline_wgs_fine_hist = dragenGermline.germline_wgs_fine_hist
        File germline_wgs_hist = dragenGermline.germline_wgs_hist
        File germline_wgs_overall_mean_cov = dragenGermline.germline_wgs_overall_mean_cov
        File tn_baf_bw = dragenTumorNormal.tn_baf_bw
        File tn_baf_seg = dragenTumorNormal.tn_baf_seg
        File tn_baf_seg_bw = dragenTumorNormal.tn_baf_seg_bw
        File tn_normal_bam = dragenTumorNormal.tn_normal_bam
        File tn_normal_bam_index = dragenTumorNormal.tn_normal_bam_index
        File tn_cnv_excluded_intervals_bed = dragenTumorNormal.tn_cnv_excluded_intervals_bed
        File tn_cnv_gff3 = dragenTumorNormal.tn_cnv_gff3
        File tn_cnv_igv_session_xml = dragenTumorNormal.tn_cnv_igv_session_xml
        File tn_cnv_metrics = dragenTumorNormal.tn_cnv_metrics
        File tn_cnv_purity_coverage_models = dragenTumorNormal.tn_cnv_purity_coverage_models
        File tn_cnv_vcf_annotated_json = dragenTumorNormal.tn_cnv_vcf_annotated_json
        File tn_cnv_vcf_annotated_json_index = dragenTumorNormal.tn_cnv_vcf_annotated_json_index
        File tn_cnv_vcf = dragenTumorNormal.tn_cnv_vcf
        File tn_cnv_vcf_index = dragenTumorNormal.tn_cnv_vcf_index
        File tn_fastqc = dragenTumorNormal.tn_fastqc
        File tn_fragment_length_hist = dragenTumorNormal.tn_fragment_length_hist
        File tn_filtered_baf_bw = dragenTumorNormal.tn_filtered_baf_bw
        File tn_filtered_vcf_annotated_json = dragenTumorNormal.tn_filtered_vcf_annotated_json
        File tn_filtered_vcf_annotated_json_index = dragenTumorNormal.tn_filtered_vcf_annotated_json_index
        File tn_filtered_vcf = dragenTumorNormal.tn_filtered_vcf
        File tn_filtered_vcf_index = dragenTumorNormal.tn_filtered_vcf_index
        File tn_insert_stats = dragenTumorNormal.tn_insert_stats
        File tn_mapping_metrics = dragenTumorNormal.tn_mapping_metrics
        File tn_microsat_diffs = dragenTumorNormal.tn_microsat_diffs
        File tn_microsat_log = dragenTumorNormal.tn_microsat_log
        File tn_microsat_normal_dist = dragenTumorNormal.tn_microsat_normal_dist
        File tn_microsat_output_json = dragenTumorNormal.tn_microsat_output_json
        File tn_microsat_tumor_dist = dragenTumorNormal.tn_microsat_tumor_dist
        File tn_pcr_model_0_log = dragenTumorNormal.tn_pcr_model_0_log
        File tn_pcr_model_1_log = dragenTumorNormal.tn_pcr_model_1_log
        File tn_pcr_model_log = dragenTumorNormal.tn_pcr_model_log
        File tn_ploidy_estimation_metrics = dragenTumorNormal.tn_ploidy_estimation_metrics
        File tn_ploidy_vcf = dragenTumorNormal.tn_ploidy_vcf
        File tn_ploidy_vcf_index = dragenTumorNormal.tn_ploidy_vcf_index
        File tn_replay_json = dragenTumorNormal.tn_replay_json
        File tn_seg = dragenTumorNormal.tn_seg
        File tn_seg_bw = dragenTumorNormal.tn_seg_bw
        File tn_somatic_mpileup_normal = dragenTumorNormal.tn_somatic_mpileup_normal
        File tn_somatic_mpileup_tumor = dragenTumorNormal.tn_somatic_mpileup_tumor
        File tn_somatic_tallies_normal = dragenTumorNormal.tn_somatic_tallies_normal
        File tn_somatic_tallies_tumor = dragenTumorNormal.tn_somatic_tallies_tumor
        File tn_sv_metrics = dragenTumorNormal.tn_sv_metrics
        File tn_sv_vcf_annotated_json = dragenTumorNormal.tn_sv_vcf_annotated_json
        File tn_sv_vcf_annotated_json_index = dragenTumorNormal.tn_sv_vcf_annotated_json_index
        File tn_sv_vcf = dragenTumorNormal.tn_sv_vcf
        File tn_sv_vcf_index = dragenTumorNormal.tn_sv_vcf_index
        File tn_time_metrics = dragenTumorNormal.tn_time_metrics
        File tn_tmb_contig_mean_cov_normal = dragenTumorNormal.tn_tmb_contig_mean_cov_normal
        File tn_tmb_contig_mean_cov_tumor = dragenTumorNormal.tn_tmb_contig_mean_cov_tumor
        File tn_tmb_coverage_metrics_normal = dragenTumorNormal.tn_tmb_coverage_metrics_normal
        File tn_tmb_coverage_metrics_tumor = dragenTumorNormal.tn_tmb_coverage_metrics_tumor
        File tn_tmb_fine_hist_normal = dragenTumorNormal.tn_tmb_fine_hist_normal
        File tn_tmb_fine_hist_tumor = dragenTumorNormal.tn_tmb_fine_hist_tumor
        File tn_tmb_hist_normal = dragenTumorNormal.tn_tmb_hist_normal
        File tn_tmb_hist_tumor = dragenTumorNormal.tn_tmb_hist_tumor
        File tn_tmb_metrics = dragenTumorNormal.tn_tmb_metrics
        File tn_tmb_overall_mean_cov_normal = dragenTumorNormal.tn_tmb_overall_mean_cov_normal
        File tn_tmb_overall_mean_cov_tumor = dragenTumorNormal.tn_tmb_overall_mean_cov_tumor
        File tn_tmb_somatic_callable_regions = dragenTumorNormal.tn_tmb_somatic_callable_regions
        File tn_tmb_trace = dragenTumorNormal.tn_tmb_trace
        File tn_tn_bw = dragenTumorNormal.tn_tn_bw
        File tn_tn_tsv = dragenTumorNormal.tn_tn_tsv
        File tn_trimmer_metrics = dragenTumorNormal.tn_trimmer_metrics
        File tn_tumor_baf_bedgraph = dragenTumorNormal.tn_tumor_baf_bedgraph
        File tn_tumor_ballele_counts = dragenTumorNormal.tn_tumor_ballele_counts
        File tn_tumor_bam = dragenTumorNormal.tn_tumor_bam
        File tn_tumor_bam_index = dragenTumorNormal.tn_tumor_bam_index
        File tn_tumor_improper_pairs = dragenTumorNormal.tn_tumor_improper_pairs
        File tn_tumor_target_counts_bw = dragenTumorNormal.tn_tumor_target_counts_bw
        File tn_tumor_target_counts_gc_corrected = dragenTumorNormal.tn_tumor_target_counts_gc_corrected
        File tn_tumor_target_counts = dragenTumorNormal.tn_tumor_target_counts
        File tn_tn_vcf = dragenTumorNormal.tn_tn_vcf
        File tn_tn_vcf_index = dragenTumorNormal.tn_tn_vcf_index
        File tn_vc_metrics = dragenTumorNormal.tn_vc_metrics
        File tn_wgs_contig_mean_cov_normal = dragenTumorNormal.tn_wgs_contig_mean_cov_normal
        File tn_wgs_contig_mean_cov_tumor = dragenTumorNormal.tn_wgs_contig_mean_cov_tumor
        File tn_wgs_coverage_metrics_normal = dragenTumorNormal.tn_wgs_coverage_metrics_normal
        File tn_wgs_coverage_metrics_tumor = dragenTumorNormal.tn_wgs_coverage_metrics_tumor
        File tn_wgs_fine_hist_normal = dragenTumorNormal.tn_wgs_fine_hist_normal
        File tn_wgs_fine_hist_tumor = dragenTumorNormal.tn_wgs_fine_hist_tumor
        File tn_wgs_hist_normal = dragenTumorNormal.tn_wgs_hist_normal
        File tn_wgs_hist_tumor = dragenTumorNormal.tn_wgs_hist_tumor
        File tn_wgs_overall_mean_cov_normal = dragenTumorNormal.tn_wgs_overall_mean_cov_normal
        File tn_wgs_overall_mean_cov_tumor = dragenTumorNormal.tn_wgs_overall_mean_cov_tumor
        File tn_wgs_somatic_callable_regions = dragenTumorNormal.tn_wgs_somatic_callable_regions
        File n_gc_metrics = dragenTumorNormal.n_gc_metrics
        File t_gc_metrics = dragenTumorNormal.t_gc_metrics
      }
}

## Task definitions
task dragenGermline {
      input {
        File ht_reference
        File normal_csv
        String normal_id
      }
    command <<<
        /opt/edico/bin/dragen_reset -r
        rm -rf /ephemeral/hash-hg38
        rm -rf /ephemeral/germline
        tar xvf ~{ht_reference}
        mv hash-hg38 /ephemeral/hash-hg38
        rm ~{ht_reference}
        mkdir -p /ephemeral/germline
        /opt/edico/bin/dragen \
          --enable-sort=true \
          --enable-map-align-output=true \
          --enable-duplicate-marking=true \
          --intermediate-results-dir=/staging \
          --enable-map-align=true \
          --enable-variant-caller=true \
          --force \
          --gc-metrics-enable=true \
          --ref-dir=/ephemeral/hash-hg38 \
          --config-file=/opt/edico/config/dragen-user-defaults.cfg \
          --lic-server=https://bdM2KvfvNyA=:uNaCli7EgMpjfMQdGjzkQT8R1lAcBlZR@license.edicogenome.com \
          --fastq-list=~{normal_csv} \
          --fastq-list-sample-id=~{normal_id} \
          --output-directory=/ephemeral/germline \
          --output-file-prefix=~{normal_id}
    >>>
    runtime {
        memory: "224 GB"
        cpu: 16
        disks: "/ephemeral 860 NVMe, /staging 2048 HDD"
        queueArn: "arn:aws:batch:us-east-1:549456635097:job-queue/somatic-dragen"
        docker: "549456635097.dkr.ecr.us-east-1.amazonaws.com/awsclinicaldragen:3.8.5"
    }
    output {
      File germline_gc_metrics = "/ephemeral/germline/~{normal_id}.gc_metrics.csv"
      File germline_baf_bw = "/ephemeral/germline/~{normal_id}.baf.bw"
      File germline_bam = "/ephemeral/germline/~{normal_id}.bam"
      File germline_bam_index = "/ephemeral/germline/~{normal_id}.bam.bai"
      File germline_fastqc = "/ephemeral/germline/~{normal_id}.fastqc_metrics.csv"
      File germline_fragment_len_hist = "/ephemeral/germline/~{normal_id}.fragment_length_hist.csv"
      File germline_filtered_vcf_baf_bw = "/ephemeral/germline/~{normal_id}.hard-filtered.baf.bw"
      File germline_filtered_vcf = "/ephemeral/germline/~{normal_id}.hard-filtered.vcf.gz"
      File germline_filtered_vcf_index = "/ephemeral/germline/~{normal_id}.hard-filtered.vcf.gz.tbi"
      File germline_insert_stats = "/ephemeral/germline/~{normal_id}.insert-stats.tab"
      File germline_mapping_metrics = "/ephemeral/germline/~{normal_id}.mapping_metrics.csv"
      File germline_pcr_model_0_log = "/ephemeral/germline/~{normal_id}.pcr-model-0.log"
      File germline_pcr_model_log = "/ephemeral/germline/~{normal_id}.pcr-model.log"
      File germline_ploidy_metrics = "/ephemeral/germline/~{normal_id}.ploidy_estimation_metrics.csv"
      File germline_ploidy_vcf = "/ephemeral/germline/~{normal_id}.ploidy.vcf.gz"
      File germline_ploidy_vcf_index = "/ephemeral/germline/~{normal_id}.ploidy.vcf.gz.tbi"
      File germline_replay_json = "/ephemeral/germline/~{normal_id}-replay.json"
      File germline_roh_bed = "/ephemeral/germline/~{normal_id}.roh.bed"
      File germline_roh_metrics = "/ephemeral/germline/~{normal_id}.roh_metrics.csv"
      File germline_time_metrics = "/ephemeral/germline/~{normal_id}.time_metrics.csv"
      File germline_trimmer_metrics = "/ephemeral/germline/~{normal_id}.trimmer_metrics.csv"
      File germline_vcf = "/ephemeral/germline/~{normal_id}.vcf.gz"
      File germline_vcf_index = "/ephemeral/germline/~{normal_id}.vcf.gz.tbi"
      File germline_vc_metrics = "/ephemeral/germline/~{normal_id}.vc_metrics.csv"
      File germline_wgs_contig_mean_cov = "/ephemeral/germline/~{normal_id}.wgs_contig_mean_cov.csv"
      File germline_wgs_coverage_metrics = "/ephemeral/germline/~{normal_id}.wgs_coverage_metrics.csv"
      File germline_wgs_fine_hist = "/ephemeral/germline/~{normal_id}.wgs_fine_hist.csv"
      File germline_wgs_hist = "/ephemeral/germline/~{normal_id}.wgs_hist.csv"
      File germline_wgs_overall_mean_cov = "/ephemeral/germline/~{normal_id}.wgs_overall_mean_cov.csv"
    }
}

task dragenTumorNormal {
      input {
        File ht_reference
        File nirvana
        File vc_systematic_noise
        File qc_coverage_region_1
        File qc_coverage_region_2
        File msi_microsatellites
        File tumor_csv
        File normal_csv
        String normal_id
        String tumor_id
        String sex
        File normal_vcf
        File normal_vcf_index
      }
    command <<<
        /opt/edico/bin/dragen_reset -r
        rm -rf /ephemeral/hash-hg38
        rm -rf /ephemeral/somatic
        rm -rf /ephemeral/nirvana-grch38
        tar xvf ~{ht_reference}
        mv hash-hg38 /ephemeral/hash-hg38
        rm ~{ht_reference}
        tar xvf ~{nirvana}
        mv nirvana-grch38 /ephemeral/nirvana-grch38
        rm ~{nirvana}
        mkdir -p /ephemeral/somatic
        /opt/edico/bin/dragen \
          --vc-enable-vaf-ratio-filter=true \
          --msi-command=tumor-normal \
          --msi-coverage-threshold=60 \
          --output-format=bam \
          --enable-sort=true \
          --enable-duplicate-marking=true \
          --enable-bam-indexing=true \
          --enable-map-align-output=true \
          --intermediate-results-dir=/staging \
          --sample-sex=~{sex} \
          --tumor-fastq-list=~{tumor_csv} \
          --tumor-fastq-list-sample-id=~{tumor_id} \
          --enable-map-align=true \
          --enable-variant-caller=true \
          --enable-sv=true \
          --enable-cnv=true \
          --qc-coverage-region-2=~{qc_coverage_region_2} \
          --qc-coverage-filters-2='mapq<10,bq<30'\
          --qc-coverage-reports-2=full_res \
          --prepend-filename-to-rgid=1 \
          --qc-coverage-reports-1=callability \
          --cnv-normal-b-allele-vcf=~{normal_vcf} \
          --vc-systematic-noise=~{vc_systematic_noise} \
          --msi-microsatellites-file=~{msi_microsatellites} \
          --enable-tmb=true \
          --qc-coverage-region-1=~{qc_coverage_region_1} \
          --qc-coverage-tag-1=tmb \
          --vc-callability-tumor-thresh=40 \
          --enable-variant-annotation=true \
          --variant-annotation-assembly=GRCh38 \
          --variant-annotation-data=/ephemeral/nirvana-grch38 \
          --qc-somatic-contam-vcf=/opt/edico/config/somatic_sample_cross_contamination_resource_hg38.vcf.gz \
          --force \
          --gc-metrics-enable=true \
          --ref-dir=/ephemeral/hash-hg38 \
          --sv-reference=/ephemeral/hash-hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
          --lic-server=https://bdM2KvfvNyA=:uNaCli7EgMpjfMQdGjzkQT8R1lAcBlZR@license.edicogenome.com \
          --fastq-list=~{normal_csv} \
          --fastq-list-sample-id=~{normal_id} \
          --output-directory=/ephemeral/somatic \
          --output-file-prefix=~{tumor_id}--~{normal_id}
        ls -lh /ephemeral/somatic
    >>>
    runtime {
        memory: "224 GB"
        cpu: 16
        disks: "/ephemeral 860 NVMe, /staging 2048 HDD"
        queueArn: "arn:aws:batch:us-east-1:549456635097:job-queue/somatic-dragen" 
	docker: "549456635097.dkr.ecr.us-east-1.amazonaws.com/awsclinicaldragen:3.8.5" 
   }
    output {
      File n_gc_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.gc_metrics_normal.csv"
      File t_gc_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.gc_metrics_tumor.csv"
      File tn_baf_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.baf.bw"
      File tn_baf_seg = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.baf.seg"
      File tn_baf_seg_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.baf.seg.bw"
      File tn_normal_bam = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.bam"
      File tn_normal_bam_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.bam.bai"
      File tn_cnv_excluded_intervals_bed = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.excluded_intervals.bed.gz"
      File tn_cnv_gff3 = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.gff3"
      File tn_cnv_igv_session_xml = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.igv_session.xml"
      File tn_cnv_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv_metrics.csv"
      File tn_cnv_purity_coverage_models = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.purity.coverage.models.tsv"
      File tn_cnv_vcf_annotated_json = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.vcf.annotated.json.gz"
      File tn_cnv_vcf_annotated_json_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.vcf.annotated.json.gz.jsi"
      File tn_cnv_vcf = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.vcf.gz"
      File tn_cnv_vcf_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.cnv.vcf.gz.tbi"
      File tn_fastqc = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.fastqc_metrics.csv"
      File tn_fragment_length_hist = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.fragment_length_hist.csv"
      File tn_filtered_baf_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.hard-filtered.baf.bw"
      File tn_filtered_vcf_annotated_json = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.hard-filtered.vcf.annotated.json.gz"
      File tn_filtered_vcf_annotated_json_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.hard-filtered.vcf.annotated.json.gz.jsi"
      File tn_filtered_vcf = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.hard-filtered.vcf.gz"
      File tn_filtered_vcf_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.hard-filtered.vcf.gz.tbi"
      File tn_insert_stats = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.insert-stats.tab"
      File tn_mapping_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.mapping_metrics.csv"
      File tn_microsat_diffs = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.microsat_diffs.txt"
      File tn_microsat_log = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.microsat_log.txt"
      File tn_microsat_normal_dist = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.microsat_normal.dist"
      File tn_microsat_output_json = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.microsat_output.json"
      File tn_microsat_tumor_dist = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.microsat_tumor.dist"
      File tn_pcr_model_0_log = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.pcr-model-0.log"
      File tn_pcr_model_1_log = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.pcr-model-1.log"
      File tn_pcr_model_log = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.pcr-model.log"
      File tn_ploidy_estimation_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.ploidy_estimation_metrics.csv"
      File tn_ploidy_vcf = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.ploidy.vcf.gz"
      File tn_ploidy_vcf_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.ploidy.vcf.gz.tbi"
      File tn_replay_json = "/ephemeral/somatic/~{tumor_id}--~{normal_id}-replay.json"
      File tn_seg = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.seg"
      File tn_seg_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.seg.bw"
      File tn_somatic_mpileup_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.somatic-mpileup-normal.txt"
      File tn_somatic_mpileup_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.somatic-mpileup-tumor.txt"
      File tn_somatic_tallies_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.somatic-tallies-normal.tallies"
      File tn_somatic_tallies_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.somatic-tallies-tumor.tallies"
      File tn_sv_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.sv_metrics.csv"
      File tn_sv_vcf_annotated_json = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.sv.vcf.annotated.json.gz"
      File tn_sv_vcf_annotated_json_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.sv.vcf.annotated.json.gz.jsi"
      File tn_sv_vcf = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.sv.vcf.gz"
      File tn_sv_vcf_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.sv.vcf.gz.tbi"
      File tn_time_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.time_metrics.csv"
      File tn_tmb_contig_mean_cov_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_contig_mean_cov_normal.csv"
      File tn_tmb_contig_mean_cov_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_contig_mean_cov_tumor.csv"
      File tn_tmb_coverage_metrics_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_coverage_metrics_normal.csv"
      File tn_tmb_coverage_metrics_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_coverage_metrics_tumor.csv"
      File tn_tmb_fine_hist_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_fine_hist_normal.csv"
      File tn_tmb_fine_hist_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_fine_hist_tumor.csv"
      File tn_tmb_hist_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_hist_normal.csv"
      File tn_tmb_hist_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_hist_tumor.csv"
      File tn_tmb_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb.metrics.csv"
      File tn_tmb_overall_mean_cov_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_overall_mean_cov_normal.csv"
      File tn_tmb_overall_mean_cov_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_overall_mean_cov_tumor.csv"
      File tn_tmb_somatic_callable_regions = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb_somatic_callable_regions.bed"
      File tn_tmb_trace = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tmb.trace.tsv"
      File tn_tn_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tn.bw"
      File tn_tn_tsv = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tn.tsv.gz"
      File tn_trimmer_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.trimmer_metrics.csv"
      File tn_tumor_baf_bedgraph = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.baf.bedgraph.gz"
      File tn_tumor_ballele_counts = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.ballele.counts.gz"
      File tn_tumor_bam = "/ephemeral/somatic/~{tumor_id}--~{normal_id}_tumor.bam"
      File tn_tumor_bam_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}_tumor.bam.bai"
      File tn_tumor_improper_pairs = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.improper.pairs.bw"
      File tn_tumor_target_counts_bw = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.target.counts.bw"
      File tn_tumor_target_counts_gc_corrected = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.target.counts.gc-corrected.gz"
      File tn_tumor_target_counts = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.tumor.target.counts.gz"
      File tn_tn_vcf = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.vcf.gz"
      File tn_tn_vcf_index = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.vcf.gz.tbi"
      File tn_vc_metrics = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.vc_metrics.csv"
      File tn_wgs_contig_mean_cov_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_contig_mean_cov_normal.csv"
      File tn_wgs_contig_mean_cov_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_contig_mean_cov_tumor.csv"
      File tn_wgs_coverage_metrics_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_coverage_metrics_normal.csv"
      File tn_wgs_coverage_metrics_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_coverage_metrics_tumor.csv"
      File tn_wgs_fine_hist_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_fine_hist_normal.csv"
      File tn_wgs_fine_hist_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_fine_hist_tumor.csv"
      File tn_wgs_hist_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_hist_normal.csv"
      File tn_wgs_hist_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_hist_tumor.csv"
      File tn_wgs_overall_mean_cov_normal = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_overall_mean_cov_normal.csv"
      File tn_wgs_overall_mean_cov_tumor = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_overall_mean_cov_tumor.csv"
      File tn_wgs_somatic_callable_regions = "/ephemeral/somatic/~{tumor_id}--~{normal_id}.wgs_somatic_callable_regions.bed"
    }
}
