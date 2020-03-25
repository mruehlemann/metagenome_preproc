#!/bin/bash
#SBATCH -p ikmb_new
#SBATCH --qos=ikmb_a
#SBATCH -c 20
#SBATCH --mem=240gb
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/data_storage/raw_data/amato_et_al/log/%A_%a.out

cd /home/sukmb276/Isilon/Metagenomes/data_storage/raw_data/amato_et_al/per_sample_FASTQ/33269/

cpus=20
ADAPTERS="/ifs/data/nfs_share/ikmb_repository/software/bbmap/37.88/resources/truseq.fa.gz"
genome_index="/ifs/data/nfs_share/sukmb276/references/bbmask_references/hg19/"
READMINLEN=60

allsamples=($(ls ERUB*_R1_001.fastq.gz | cut -d '_' -f 1))

#sample="ERR2144929"
sample=${allsamples[$SLURM_ARRAY_TASK_ID]}

module load BBMap/38.57
module load Spades/3.9.0
module load Java/1.8.0

task_cpus=5

left=$(ls ${sample}_*_R1_001.fastq.gz)
right=$(ls ${sample}_*_R2_001.fastq.gz)

left_trimmed=tmp_${sample}.trimmed.R1.fastq.gz
right_trimmed=tmp_${sample}.trimmed.R2.fastq.gz
unpaired_trimmed=tmp_${sample}.trimmed.RU.fastq.gz
left_nophix=tmp_${sample}.nophix.R1.fastq.gz
right_nophix=tmp_${sample}.nophix.R2.fastq.gz
unpaired_nophix=tmp_${sample}.nophix.RU.fastq.gz
left_decon=tmp_${sample}.decon.R1.fastq.gz
right_decon=tmp_${sample}.decon.R2.fastq.gz
unpaired_decon=tmp_${sample}.decon.RU.fastq.gz

left_clean=${sample}.clean.R1.fastq.gz
right_clean=${sample}.clean.R2.fastq.gz
merged=tmp_${sample}.merged.RU.fastq.gz
unpaired_clean=${sample}.clean.RU.fastq.gz

###QC
### basic stats of input files (can be omitted)
reformat.sh threads=${task.cpus} in=${left} in2=${right}

### bbduk does adaptor trimming and artifact removal
bbduk.sh stats=$bbduk_adapter_stats threads=${task_cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=${unpaired_trimmed} ref=${ADAPTERS} ktrim=r k=23 mink=11 hdist=1 minlength=${READMINLEN} tpe tbo
bbduk.sh stats=$bbduk_artifact_stats threads=${task_cpus} in=${left_trimmed} in2=${right_trimmed} k=31 ref=artifacts,phix ordered cardinality out1=${left_nophix} out2=${right_nophix} minlength=${READMINLEN}
bbduk.sh threads=${task_cpus} in=${unpaired_trimmed}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_nophix} minlength=${READMINLEN}

### removal of host reads (human is used as reference), see for details: http://seqanswers.com/forums/archive/index.php/t-42552.html
bbwrap.sh -Xmx23g threads=${task_cpus} minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=20 minlength=${READMINLEN} in=${left_nophix},${unpaired_nophix} in2=${right_nophix},NULL path=${genome_index} outu1=${left_decon} outu2=${right_decon} outu=${unpaired_decon}

### merging of overlapping reads
bbmerge.sh threads=${task_cpus} in1=${left_decon} in2=${right_decon} out=${merged} outu1=${left_clean} outu2=${right_clean} mininsert=${READMINLEN}

### cleanup steps
zcat ${merged} ${unpaired_decon} | gzip -c > ${unpaired_clean}
reformat.sh threads=${task_cpus} in=${unpaired_clean}
reformat.sh threads=${task_cpus} in1=${left_clean} in2=${right_clean}

##SPADES ASSEMBLY
### assembly
spades.py --meta --pe1-1 $left_clean --pe1-2 $right_clean --pe1-s $unpaired_clean -k 21,33,55 -o ${sample}_spades_out -t ${task_cpus}
### renaming of contigs
outcontigs=${sample}.contigs.fasta
awk -v id=$id '/^>/{print ">Spades_"id"_contig_"++i; next}{print}' ${sample}_spades_out/scaffolds.fasta > $outcontigs

rm tmp_${sample}_*

exit








