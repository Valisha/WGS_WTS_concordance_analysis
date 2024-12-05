version 1.0

workflow qc_metrics_parse {
    meta {
        description: "pulls QC metrics into json and uploads tsv"
    }
    input {
        String sample_id
        File mapping_metrics
        File wgs_cov_tumor
        File wgs_cov_normal
        File dragen_qc_metrics_py
        File cnv_metrics_file
        File gc_normal
        File gc_tumor
    }

call qc_json {
    input:
        sample_id = sample_id,
        mapping_metrics = mapping_metrics,
        wgs_cov_tumor = wgs_cov_tumor,
        wgs_cov_normal = wgs_cov_normal,
        dragen_qc_metrics_py = dragen_qc_metrics_py,
        cnv_metrics_file = cnv_metrics_file
        gc_normal = gc_normal
        gc_tumor = gc_tumor
        
    }
    output {
        File final_parsed_qc_metrics_tsv = qc_json.parsed_qc_metrics_tsv
    }
}

task qc_json {
    input {
        String sample_id
        File mapping_metrics
        File wgs_cov_tumor
        File wgs_cov_normal
        File dragen_qc_metrics_py
        File cnv_metrics_file
        File gc_normal
        File gc_tumor
    }
    command {
        python3 \
        ~{dragen_qc_metrics_py} \
        --mapping_metrics ~{mapping_metrics} \
        --wgs_cov_tumor ~{wgs_cov_tumor} \
        --wgs_cov_normal ~{wgs_cov_normal} \
        --cnv_metrics ~{cnv_metrics_file} \
        -gn ~{gc_normal} \
        -gt ~{gc_tumor} \
        -op ~{sample_id}.tsv
    }
    output {
        File parsed_qc_metrics_tsv = "~{sample_id}.tsv"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        queueArn: "arn:aws:batch:us-east-1:876086242506:job-queue/dragen-queue-noF1"
    }
}

