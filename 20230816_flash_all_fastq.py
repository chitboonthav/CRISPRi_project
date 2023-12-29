#!/usr/bin/python3
#run in /Users/chutikarnchitboonthavisuk/Desktop/Nextseq-data_062723
import glob
import os

filenames = glob.glob("*_R1_001.fastq.gz")
for raw_name in filenames:
    sample = raw_name.split("_")
    sample_name = sample[0]
    sample_number = sample[1]
    print(str(sample_name))
  
    os.system("mkdir merged_sequence/{name}_{num}".format(name=sample_name, num=sample_number))
    os.system("/Users/chutikarnchitboonthavisuk/Desktop/Work/NGS-Program/FLASH-1.2.11/flash \
    {name}_{num}_L002_R1_001.fastq.gz \
    {name}_{num}_L002_R2_001.fastq.gz \
    --output-prefix {name}_{num} \
    --output-directory merged_sequence/{name}_{num} | tee {name}_{num}_flash.log".format(name=sample_name,num=sample_number))