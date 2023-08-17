## quality control of raw data
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N fastp
fastp -i fq/$i/$i\_1.fastq.gz -o cleandata/$i.1.clean.fastq.gz -I fq/$i/$i\_2.fastq.gz -O cleandata/$i.2.clean.fastq.gz \
-h cleandata/output/$i\.fastp.html" >>shell/$i.fastp.sh;
done
for i in `cat sample_List.txt`;do qsub shell/$i.fastp.sh; done

## remove host sequences
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N bowtie2
bowtie2 -p 4 -x /public1/data/liusheng/data/GRCh38/hg38 -1 cleandata/$i.1.clean.fastq.gz -2 cleandata/$i.2.clean.fastq.gz \
--un-conc-gz rmhost/$i/$i.rmhost.fq.gz" >>shell/$i.rmhost.sh;
done
for i in `cat sample_List.txt`;do qsub shell/$i.rmhost.sh; done

## assembly using megahit
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#!/bin/bash
#PBS -N megahit
megahit -1 rmhost/$i/$i.rmhost_1.fq.gz -2 rmhost/$i/$i.rmhost_2.fastq.gz -o megahit/output/$i" >>shell/$i.megahit.sh;
done
for i in `cat sample_List.txt`;do qsub shell/$i.megahit.sh; done

## drop sequences with length shorter than 1000bp
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N seqtk
seqtk seq -L 1000 megahit/output/$i/final.contigs.fa > assembly_filter_results_megahit/$i.contigs.fna" >>shell/$i.seqtk.sh;
done
for i in `cat sample_List.txt`;do qsub shell/$i.seqtk.sh; done

## antismash to analyse BGCs
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N antismash
antismash --taxon bacteria --genefinding-tool prodigal assembly_filter_results_megahit/$i.contigs.fna \
-c 8 --output-dir antismash/output/$i" >>shell/$i.antismash.sh;
done
for i in `cat sample_List.txt`;do qsub shell/$i.antismash.sh; done

## extract the BGCs classes and number 
for i in `cat sample_List.txt`;do 
grep ' Region ' antismash/output/$i/index.html |awk '{print $8}' |sort |uniq -c >result/$i.BGC.txt; done

## calculate the abundance of contigs
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N salmon
salmon index -t assembly_filter_results_megahit/$i.contigs.fna -i salmon/output/$i.contig;
salmon quant -p 4 -l A -i salmon/output/$i.contig -1 rmhost/$i/$i.rmhost_1.fq.gz -2 rmhost/$i/$i.rmhost_2.fq.gz \
--meta -o salmon/output/$i.quant" >>shell/$i.salmon.sh; done
for i in `cat sample_List.txt`;do qsub shell/$i.salmon.sh; done

## extract the BGCs-derived contigs and calcualte the TPMs of BGCs
for i in `cat sample_List.txt`;do 
ls antismash/output/$i/k*.region0*.gbk |cut -d '/' -f 4 |cut -d '.' -f 1 |sort |uniq >contig.bgc/$i.contig;done
for i in `cat sample_List.txt`;do
        for j in `cat contig.bgc/$i.contig`;do
        grep "product=" antismash/output/$i/$j.region0*.gbk |cut -d '=' -f 2 |sed 's/"//g' |sort |uniq >contig.bgc.result/$i.$j.bgc;done;
done
for i in `cat sample_List.txt`;do
        for j in `cat contig.bgc/$i.contig`;do
        echo "awk '{print \$1\"\\\t\"\$2\"\\\t$j\"}' contig.bgc.result/$i.$j.bgc >>contig.bgc.result/$i.bgc.txt" >>total.sh;done
done
sh total.sh;
for i in `cat sample_List.txt`;do
        echo "perl -e 'open I,\"../salmon/output/$i.quant/quant.sf\";<I>;while(<I>){@arr=split;\$hash{\$arr[0]}=\$arr[3]};open F,\"contig.bgc.result/$i.bgc.txt\"; \
        while(<F>){@b=split;print \"\$b[0]\\\t\$b[1]\\\t\$hash{\$b[1]}\\\n\"}' \
        >contig.bgc.result/$i.bgc.TPM.txt" >>getTPM.sh;done
sh getTPM.sh;
for i in `cat sample_List.txt`;do
        awk '{sum[$1]+=$3}END{for(c in sum){print c,sum[c]}}' contig.bgc.result/$i.bgc.TPM.txt >contig.bgc.result/$i.bgc.TPM.sum.txt;
done

## mergy the BGCs results
cat contig.bgc.result/*.bgc.TPM.sum.txt |awk '{print $1}' |sort |uniq >all.uniq.BGCs;
perl mergy.TMP.pl allsample.BGCs.TPM.xls

## calculate the species abundance using MetaPhlAn
for i in `cat sample_List.txt`;do echo "#!/bin/bash
#PBS -N metaphlan
metaphlan rmhost/$i/$i.rmhost_1.fq.gz,rmhost/$i/$i.rmhost_2.fq.gz --bowtie2db /data/data1/liusheng/database/metaphlan/mpa_vJan21/ \
-x mpa_vJan21_CHOCOPhlAnSGB_202103 --nproc 4 --input_type fastq --bowtie2out $i.bowtie2.bz2 \
-o metaphlan/output/$i.metagenome.txt" >>shell/$i.metaphlan.sh

## mergy the metaphlan results of multiple samples
merge_metaphlan_tables.py metaphlan/output/*metagenome.txt >merge.metagenome.txt
sed '1'd merge.metagenome.txt |sed '1s/.metagenome//g' >merge.metagenome.xls
awk '$0~/clade_name/||$0~/s__/{print $0}' merge.metagenome.xls |cut -d '|' -f 7- >merge.metagenome.species.xls
sed 's/^s__//g' merge.metagenome.species.xls | sed '/|t__/'d |sed 's/clade_name//' >merge.metagenome.species.v1.xls




