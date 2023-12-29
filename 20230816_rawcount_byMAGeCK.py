#!/usr/bin/python3
#run in /Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/Nextseq-data_062723

import glob
import os

filenames = glob.glob("*_R1_001.fastq.gz")
for raw_name in filenames:
    sample = raw_name.split("_")
    sample_name = sample[0]
    sample_number = sample[1]
    print(str(sample_name))
  
  
    os.system("mkdir rawcount/{name}_{num}".format(name=sample_name,num=sample_number))
    os.chdir("/Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/Nextseq-data_062723/rawcount/{name}_{num}".format(name=sample_name,num=sample_number))
    os.system("mageck count -l /Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/20220816_CRISPRi_paper_data/final_combine_sgRNA.csv \
    --fastq /Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/Nextseq-data_062723/merged_sequence/{name}_{num}/{name}_{num}.extendedFrags.fastq \
    --sample-label {name}_{num} \
    -n {name}_{num}".format(name=sample_name,num=sample_number))
    os.chdir("/Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/Nextseq-data_062723")  