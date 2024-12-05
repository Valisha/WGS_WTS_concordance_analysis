#!/bin/bash
module load samtools

if [ $# -ne 1 ]; then
    echo "Usage: $0 OUT_DIR" >&2
    exit 1
fi

out_dir=$1

zcat $out_dir/chr1.out_small.seqz.gz > $out_dir/sample.seqz

zcat $out_dir/chr2.out_small.seqz.gz $out_dir/chr3.out_small.seqz.gz $out_dir/chr4.out_small.seqz.gz $out_dir/chr5.out_small.seqz.gz $out_dir/chr6.out_small.seqz.gz $out_dir/chr7.out_small.seqz.gz $out_dir/chr8.out_small.seqz.gz $out_dir/chr9.out_small.seqz.gz $out_dir/chr10.out_small.seqz.gz $out_dir/chr11.out_small.seqz.gz $out_dir/chr12.out_small.seqz.gz $out_dir/chr13.out_small.seqz.gz $out_dir/chr14.out_small.seqz.gz $out_dir/chr15.out_small.seqz.gz $out_dir/chr16.out_small.seqz.gz $out_dir/chr17.out_small.seqz.gz $out_dir/chr18.out_small.seqz.gz $out_dir/chr19.out_small.seqz.gz $out_dir/chr20.out_small.seqz.gz $out_dir/chr21.out_small.seqz.gz $out_dir/chr22.out_small.seqz.gz $out_dir/chrX.out_small.seqz.gz | gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' >> $out_dir/sample.seqz

bgzip $out_dir/sample.seqz

tabix -f -s 1 -b 2 -e 2 -S 1 $out_dir/sample.seqz.gz
