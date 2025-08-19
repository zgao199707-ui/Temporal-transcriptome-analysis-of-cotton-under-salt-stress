# Prerequisites
The default conda environment has already been set up, and required software is installed and available in the current environment. Necessary directories have also been created.

# Directory structure
Reference genome files: /home/agis/huguanjing_group/gaozhan/time_RNAseq/raw/ref
Main loop script:       /home/agis/huguanjing_group/gaozhan/time_RNAseq/src/main.sh
FastQC script:          /home/agis/huguanjing_group/gaozhan/time_RNAseq/src/GetQuality.sh
fastp script:           /home/agis/huguanjing_group/gaozhan/time_RNAseq/src/TrimFqFile.sh
Kallisto script:        /home/agis/huguanjing_group/gaozhan/time_RNAseq/src/kallistoRun.sh
Working directory:      /home/agis/huguanjing_group/gaozhan/time_RNAseq/exp/20240913/run
Raw sequencing data: /home/agis/huguanjing_group/gaozhan/time_RNAseq/01.RawData/
SLURM command alias：alias mysrun='srun --nodes 1 --ntasks 1 --cpus-per-task 1 -p low,big --mem-per-cpu=8G'

# Step 1: Organize raw data
cp /home/agis/huguanjing_group/gaozhan/time_RNAseq/01.RawData/*/*.fq.gz dat/timeRNAseq/ 

# Step 2: Generate quality reports for raw data using FastQC
for i in `ls dat/timeRNAseq/*.fq.gz` ; do 
    sh ../../src/main.sh $(basename "${i%.fa.gz}") run/GetQuality.sh $i res/RawQC
done

## Run MultiQC to summarize raw FastQC reports
mysrun -o log/multiQCtrimedFq1.out -e log/multiQCtrimedFq1.err \
    multiqc -o res/RawQC/RawQC res/RawQC/

# Step 3: Perform quality trimming with fastp
for i in `ls dat/timeRNAseq/*.fq.gz | sed 's/_[12].fq.gz//' | sort | uniq` ; do 
    sh ../../src/main.sh $(basename "${i}") run/TrimFqFile.sh $i dat/TrimedFile
done

# Step 4: Summarize fastp trimming results
printf 'Raw_reads\tClean_reads\tQ30(%%)\tQ20(%%)\tGC\tClean base\tSample ID\n' > res/TrimSummary.txt
rm tmp/mytempfpsummary.txt   # temporary file for storing parsed stats

# Extract relevant metrics from fastp JSON output and compile into summary table
for i in `ls dat/TrimedFile/*.json`; do 
    grep "total_reads" $i | head -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    grep "total_reads" $i | head -n 2 | tail -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    grep "q20_rate" $i | head -n 2 | tail -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    grep "q30_rate" $i | head -n 2 | tail -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    grep "gc_content" $i | head -n 2 | tail -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    grep "total_bases" $i | head -n 2 | tail -n 1 | sed "s/^.*:\(.*\)/\1/" | sed "s/,$//" >> tmp/mytempfpsummary.txt
    echo "$(basename $i)" >> tmp/mytempfpsummary.txt
    paste -s -d '\t' tmp/mytempfpsummary.txt >> res/TrimSummary.txt
done

# Step 5: Generate quality reports for cleaned data using FastQC
for i in `ls dat/TrimedFile/*.clean.gz` ; do 
    sh ../../src/main.sh $(basename "${i%.clean.gz}") run/GetQuality.sh $i res/CleanQC
done

# Run MultiQC to summarize clean FastQC reports
mysrun -o log/multiQCtrimedFq2.out -e log/multiQCtrimedFq2.err \
    multiqc -o res/CleanQC/CleanQC res/CleanQC/


# Transcriptome Quantification with Kallisto
After filtering and quality control, proceed with transcript quantification, using `kallisto index` and `kallisto quant`. Before alignment, build the genome index file. Also prepare a chromosome-length tab-separated file if needed. Then perform paired-end quantification with `kallisto quant`.

# Step 1. Build kallisto index
srun --nodes 1 --ntasks 1 --cpus-per-task 8 -p low,big \
    -e log/index.err -o log/index.out \
    kallisto index -t 8 \
    -i /home/agis/huguanjing_group/gaozhan/time_RNAseq/raw/ref/UTX.idx \
    /home/agis/huguanjing_group/gaozhan/time_RNAseq/raw/ref/Ghirsutum_527_v2.0.fa

# Step 2. Create output directory
mkdir res/KallistoRes

# Step 3. Run kallisto quant on paired-end clean data
for i in `ls dat/TrimedFile/*.clean.gz | sed 's/_[12].clean.gz//' | sort | uniq`; do
    sh ../../src/main.sh $(basename "${i}") run/kallistoRun.sh $i res/KallistoRes/;
done

# Step 4. Summarize mapping rate
rm res/kallistoMapped.tsv
for i in `ls log/kallisto*.err`; do 
    echo -e "$(echo ${i%.err} | sed 's|log/kallisto_||')\t$(grep -Eo '[0-9]+.[0-9]+%' $i | tail -n1)" \
    >> res/kallistoMapped.tsv
done

# Step 5. Collect TPM results
cut -f1 res/KallistoRes/S9-3/abundance.tsv > res/kallisto_tpm.tsv
for i in `ls res/KallistoRes/*/abundance.tsv`; do
    cp res/kallisto_tpm.tsv tmp/test2.tsv
    myname="$(echo $i | sed 's|res/KallistoRes/||' | sed 's|/abundance.tsv||')"
    cut -f5 $i | sed "s/tpm/$myname/" > tmp/test.tsv
    paste tmp/test2.tsv tmp/test.tsv > res/kallisto_tpm.tsv
done

# Step 6. Collect estimated counts
cut -f1 res/KallistoRes/S9-3/abundance.tsv > res/kallisto_estCount.tsv
for i in `ls res/KallistoRes/*/abundance.tsv`; do
    cp res/kallisto_estCount.tsv tmp/test2.tsv
    myname="$(echo $i | sed 's|res/KallistoRes/||' | sed 's|/abundance.tsv||')"
    cut -f4 $i | sed "s/est_counts/$myname/" > tmp/test.tsv
    paste tmp/test2.tsv tmp/test.tsv > res/kallisto_estCount.tsv
done
