# Temporal-transcriptome-analysis-of-cotton-under-salt-stress
Analysis of salt stress time-series transcriptome data
1.Upstream RNA-seq analysis with kallisto
# Step 1: Organize raw data
cp /home/agis/huguanjing_group/gaozhan/time_RNAseq/01.RawData/*/*.fq.gz dat/timeRNAseq/ 

# Step 2: Generate quality reports for raw data using FastQC
# (edit input/output paths in GetQuality.sh if needed)
for i in `ls dat/timeRNAseq/*.fq.gz` ; do 
    sh ../../src/main.sh $(basename "${i%.fa.gz}") run/GetQuality.sh $i res/RawQC
done

# Run MultiQC to summarize raw FastQC reports
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
