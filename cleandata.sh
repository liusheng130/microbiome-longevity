for i in `cat  sample_List.txt`;do echo "#!/bin/bash
#PBS -N fastp
fastp -i fq/$i/$i\_1.fastq.gz -o cleandata/$i.1.clean.fastq.gz -I fq/$i/$i\_2.fastq.gz -O cleandata/$i.2.clean.fastq.gz -h cleandata/output/$i\.fastp.html" >>shell/$i.fastp.sh;
done

for i in `cat  sample_List.txt`;do qsub shell/$i.fastp.sh; done
