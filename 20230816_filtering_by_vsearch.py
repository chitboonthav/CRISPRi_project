#!/usr/bin/python3
#run in /Users/chutikarnchitboonthavisuk/Desktop/Nextseq-data_062723/

import glob
import os

filenames = glob.glob("*_R1_001.fastq.gz")
for raw_name in filenames:
    sample = raw_name.split("_")
    sample_name = sample[0]
    sample_number = sample[1]
    print(str(sample_name))
  
    os.system("vsearch -fastq_filter /Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/Nextseq-data_062723/merged_sequence/{name}_{num}/{name}_{num}.extendedFrags.fastq \
    -fastq_maxee_rate 0.1 -fastqout {name}_{num}.filtered_e0.1.fastq".format(name=sample_name, num =sample_number))