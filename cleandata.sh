for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N fastp
fastp -i fq/$i/$i\_1.fastq.gz -o cleandata/$i.1.clean.fastq.gz -I fq/$i/$i\_2.fastq.gz -O cleandata/$i.2.clean.fastq.gz \
-h cleandata/output/$i\.fastp.html" >>shell/$i.fastp.sh;
done

for i in `cat sample_List.txt`;do qsub shell/$i.fastp.sh; done

for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N bowtie2
bowtie2 -p 4 -x /public1/data/liusheng/data/GRCh38/hg38 -1 cleandata/$i.1.clean.fastq.gz -2 cleandata/$i.2.clean.fastq.gz \
--un-conc-gz rmhost/$i/$i.rmhost.fq.gz" >>shell/$i.rmhost.sh;
done

for i in `cat sample_List.txt`;do qsub shell/$i.rmhost.sh; done

