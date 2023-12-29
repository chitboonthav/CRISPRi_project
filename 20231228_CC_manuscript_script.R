#########################################################################################
## Project: Discovering bacterial factors governing bacteriophage host range systematically and comprehensively
## Script purpose: For the analysis of CRISPRi experiment in K12, B, and O121 strains
## Date: 2023-12-28
## Author: Chutikarn Chitboonthavisuk
#########################################################################################

# Library needed ------
library(dplyr)
library(ggplot2)
library(tidyr)
library(genefilter)
library(RColorBrewer)
library(GGally)
library(ggforce)
library(edgeR)
library(limma)
library(multtest)
library(metap)
library(knitr)
library(ggrepel)
library(gridExtra)
library(readxl)
library(seqinr)
library(ggpubr)
set.seed(1)

## Load countable and subset for matched sgRNA -----
# this is to check the distribution and coverage of sgRNA transformed in 10G and targeted host without any induction of CRISPRi
tbl1 = read.csv(file = "sgRNAdistribution_original.csv", header = TRUE) 
sgRNAall = read.csv(file = "splitsgRNA.csv", header = TRUE) # this is to get the sgRNA list and name
colnames(sgRNAall) <- c("loc_1", "loc_2", "loc_3", "Gene")
# because the original library contains matched, mismatched, and randomized sgRNA 
sgRNAall$matched = ifelse(is.na(sgRNAall$loc_3) & sgRNAall$Gene != "randomized", "mismatched", "matched")
sgRNAmatched = filter(sgRNAall, matched == "matched") # to select only matched sgRNA for this study 
sgRNAmatched$sgRNA = ifelse(is.na(sgRNAmatched$loc_3), (paste0(sgRNAmatched$loc_1,"-", sgRNAmatched$loc_2)),
                            paste0(sgRNAmatched$loc_1,"-",sgRNAmatched$loc_2,"-",sgRNAmatched$loc_3))
sgRNAmatched = sgRNAmatched[, grep(("sgRNA|Gene"), names(sgRNAmatched))] # to subset for only matched data
tbl2 = merge(sgRNAmatched, tbl1, by = c("sgRNA","Gene")) # this is the countable for only matched sgRNA

## Check distribution of no phage for the distribution of the whole library ----
# this code was made based on the lab strains involved in the experiment: 10G,BW25113,MG1655,BL21
function_distribution <- function(tbl, a)
{
  if (a == "10G") {
    x = tbl1[, grep(("sgRNA|Gene|10G"), names(tbl1))]
  } else if(a == "BW25113") {
    x = tbl1[, grep(("sgRNA|Gene|BW25113"), names(tbl1))]
  } else if (a == "MG1655") {
    x = tbl1[, grep(("sgRNA|Gene|MG1655"), names(tbl1))]
  } else {
    x = tbl1[, grep(("sgRNA|Gene|BL21"), names(tbl1))]
  }
  if (a == "10G") {
    x_mod = gather(x, "sample", "rawReads", 3:4)
  } else {
    x_mod = gather(x, "sample", "rawReads", 3:5)
  }
  sum = aggregate(x_mod$rawReads, list(x_mod$sample), FUN=sum)
  colnames(sum) = c("sample", "total")
  x_mod1 = merge(x_mod, sum, by.x = "sample", all = TRUE)
  col_order = c("sgRNA", "Gene", "sample", "rawReads", "total")
  x_mod1 <- x_mod1[,col_order]
  x_mod1$cutoff3 = ifelse(x_mod1$rawReads < 3, "reads < 3", "reads >= 3")
  x_filtered = filter(x_mod1, cutoff3 == "reads >= 3")
  x_filtered$normreads = x_filtered$rawReads/x_filtered$total
  x_wide = x_filtered[, which(colnames(x_filtered) %in% c("sgRNA","Gene","sample","normreads"))]
  x_wide =pivot_wider(x_wide, names_from = sample, values_from = normreads)
  
  allsgRNA = tbl1[,1:2]
  x_wide_merge = merge(allsgRNA, x_wide, by = c("sgRNA", "Gene"), all = TRUE)
  x_cal = x_wide_merge
  if (a == "BW25113"| a== "MG1655"| a=="BL21") {
    x_cal$avg = rowMeans(x_cal[,3:5], na.rm = TRUE) 
    x_cal$sd =rowSds(x_cal[,3:5], na.rm = TRUE)
  } else {
    x_cal$avg = rowMeans(x_cal[,3:4], na.rm = TRUE) 
    x_cal$sd =rowSds(x_cal[,3:4], na.rm = TRUE)
  }

  if ( a == "10G") {
    p1 = ggplot(x_cal, aes(x=avg)) +
      geom_histogram(bins = 100) +
      labs(title = "all sgRNA distribution in 10G WT",
           x = "average norm reads",
           y = "number of sgRNA count") +
      theme_classic()
    return(p1)
    
  } else if ( a == "BW25113") {
    p1 = ggplot(x_cal, aes(x=avg)) +
      geom_histogram(bins = 100) +
      labs(title = "all sgRNA distribution in BW25113",
           x = "average norm reads",
           y = "number of sgRNA count") +
      theme_classic()
    return(p1)
    
  } else if (a == "MG1655") {
    p1 = ggplot(x_cal, aes(x=avg)) +
      geom_histogram(bins = 100) +
      labs(title = "all sgRNA distribution in MG1655",
           x = "average norm reads",
           y = "number of sgRNA count") +
      theme_classic()
    return(p1)
    
  } else {
    p1 = ggplot(x_cal, aes(x=avg)) +
      geom_histogram(bins = 100) +
      labs(title = "sgRNA distribution in BL21",
           x = "all average norm reads",
           y = "number of sgRNA count") +
      theme_classic()
    return(p1)
  }
}

### run "function_distribution" to check the distribution ------
distribution = function_distribution(tbl = tbl1, a = "10G")
WT10G = distribution$data
hist_10G = print(distribution)

distribution = function_distribution(tbl = tbl1, a = "BW25113")
BWdist = distribution$data
hist_BWdist = print(distribution)

distribution = function_distribution(tbl = tbl1, a = "MG1655")
MGdist = distribution$data
hist_MGdist = print(distribution)

distribution = function_distribution(tbl = tbl1, a = "BL21")
BLdist = distribution$data
hist_BLdist = print(distribution)

### gather table and combine distribution -----
WT10G = WT10G[,grep(("sgRNA|Gene|avg"), names(WT10G))]
BW = BWdist[,grep(("sgRNA|Gene|avg"), names(BWdist))]
MG = MGdist[,grep(("sgRNA|Gene|avg"), names(MGdist))]
BL = BLdist[,grep(("sgRNA|Gene|avg"), names(BLdist))]

colnames(WT10G) = c("sgRNA", "Gene", "avg_10G")
colnames(BW) = c("sgRNA", "Gene", "avg_BW")
colnames(MG) = c("sgRNA", "Gene", "avg_MG")
colnames(BL) = c("sgRNA", "Gene", "avg_BL")

# at this point, we narrowed down to BW25113 and BL21 in comparison to 10G 
# There is no cutoff or filtering applied at this step

mergetbl = BW
mergetbl = merge(mergetbl, BL, by = c("sgRNA", "Gene"))
mergetbl = merge(mergetbl, WT10G, by = c("sgRNA", "Gene"))
colnames(mergetbl) = c("sgRNA", "Gene", "BW25113 with dCas9", "BL21 with dCas9", "WT_10G without dCas9" )
mergetbl_long = gather(mergetbl, "sample", "avg", 3:5)

p1 = ggplot(mergetbl_long, aes(x=avg, fill = sample)) +
  geom_histogram(bins = 500, position = "identity", alpha = 0.35) +
  labs(title = "All sgRNA distribution overnight culture, filtered at rawread cutoff of 3, no induction",
       x = "average norm reads",
       y = "number of sgRNA count") +
  facet_zoom(xlim = c(0,0.0001)) +
  theme_classic()
p1

# to save file in pdf
pdf(filename ="20221126_allsgRNAdist_noinduction.pdf",height = 6.5, width = 8.5)
p1
dev.off()

### check distribution of randomized sgRNA in 10G, BW25113, and BL21 -------
# use mergetbl
mergetbl1 = mergetbl
random = filter(mergetbl1, Gene == "randomized")
randomz_long = gather(random, "sample", "avg", 3:5)
p1 = ggplot(randomz_long, aes(x=avg, fill = sample)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.45) +
  labs(title = "Randomized sgRNA distribution no phage, no induction",
       x = "average norm reads",
       y = "number of sgRNA count") +
  facet_zoom(xlim = c(0,0.0002)) +
  theme_classic()
p1
pdf(filename ="20220816_combine_random_plot.pdf", height = 6.5, width = 8.5)
p1
dev.off()

### check the distribution of matched sgRNA only ---------
#use mergetbl
mergetbl2 = mergetbl
mergematched = merge(sgRNAmatched, mergetbl2, by = "sgRNA", all.x = TRUE)
mergematched = mergematched[,-2]
matchedlong = gather(mergematched, "sample", "avg", 3:5)
p1 = ggplot(matchedlong,aes(x=avg, fill = sample)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.45) +
  labs(title = "Matched sgRNA distribution no phage, no induction",
       x = "average norm reads",
       y = "number of sgRNA count") +
  facet_zoom(xlim = c(0,0.0002)) +
  theme_classic()
p1
pdf(filename ="20220816_combine_matched_plot.pdf",height = 6.5, width = 8.5)
p1
dev.off()

# Phage selection ---------
# I want to look at the changes in randomized sgRNA first to make sure that the data is valid
# for further determination of significant changes in targeted sgRNA related to phage infection
# To do this we need to load a new table made for phage selection experimental data

tbl3 = read.csv(file = "phageselection_labstrain.csv", header = TRUE)
tbl4 = merge(sgRNAmatched, tbl3, by = c("sgRNA","Gene")) # I will only look at matched sgRNA effect

# subsetting to make countable and run edgeR (this includes both before and after phage selection)
function_countablemaking = function(tbl, a){
  if (a == "BW25113") {
    x = tbl[, grep(("sgRNA|Gene|BW"), names(tbl))]
  } else if (a == "MG1655") {
    x = tbl[, grep(("sgRNA|Gene|MG"), names(tbl))]
  } else {
    x = tbl[, grep(("sgRNA|Gene|BL"), names(tbl))]
  }
  mergecount = merge(sgRNAmatched, x, by = c("sgRNA", "Gene"), all.x = TRUE)
}

countable_data = function_countablemaking(tbl3, "BW25113")
count_BW = countable_data
row.names(count_BW) = count_BW$sgRNA
count_BW = count_BW[,3:8]

countable_data = function_countablemaking(tbl3, "BL21")
count_BL = countable_data
row.names(count_BL) = count_BL$sgRNA
count_BL = count_BL[,3:8]

# running edgeR 
geneid = sgRNAmatched 
row.names(geneid) = geneid$sgRNA
geneid = geneid[-2]

function_edgeR_GLM = function(counttable, cutoffCPM) {
  group = c(rep("control", 3), rep("treatment", 3)) #this is the setup that considered unpaired controls and treatments 
  DGE_list = DGEList(count = counttable, genes = geneid, group = group)
  keep <- rowSums(cpm(DGE_list$count[,1:3])> cutoffCPM ) >=3 # this step filters the sgRNA *before* selection to pass CPM that we indicate
  # we only filter the samples before the selection because we want a decent amount of reads to present before the selection
  # we don't filter the samples after the selection because dropping out sgRNAs are expected
  DGE_list <- DGE_list[keep, , keep.lib.size = FALSE]
  DGE_list <- estimateDisp(DGE_list)
  treatment = as.factor(DGE_list$samples$group)
  design = model.matrix(~treatment)
  colnames(design)[2] <- "treatment"
  fit <- glmQLFit(DGE_list, design)
  lrt <- glmQLFTest(fit, coef = 2)
} # this function processes 

test_lrt = function_edgeR_GLM(counttable = count_BW, cutoffCPM = 15)
res_lrt_BW = test_lrt$table

test_lrt = function_edgeR_GLM(counttable = count_BL, cutoffCPM = 15)
res_lrt_BL = test_lrt$table
# result is the logFC generated by edgeR for each sgRNA between control and treatment

# this section merge sgRNA with the lrt table so that we know what genes those sgRNA targets
res_lrt_BW$sgRNA = row.names(res_lrt_BW)
res_lrt_BW_merge = merge(sgRNAmatched, res_lrt_BW, by = "sgRNA", all.x = TRUE)

res_lrt_BL$sgRNA = row.names(res_lrt_BL)
res_lrt_BL_merge = merge(sgRNAmatched, res_lrt_BL, by = "sgRNA", all.x = TRUE)

lrt_rand_BW = filter(res_lrt_BW_merge, Gene == "randomized") # this is to subset the randomized dataset 
lrt_rand_BL = filter(res_lrt_BL_merge, Gene == "randomized")

write.csv(x = res_lrt_BW_merge, file = "res_lrt_BW_merge_CPM15.csv",row.names = FALSE)
write.csv(x = res_lrt_BL_merge, file = "res_lrt_BL_merge_CPM15.csv",row.names = FALSE)

# Now, we need to create the distribution of randomized sgRNA 
# this distribution is for assessing the effect of randomized sgRNAs and to use as a reference for determining significance

# Import data for this randomized filter -----------
BW_indi_sgRNA = read.csv(file = "res_lrt_BW_merge_CPM15.csv", header = TRUE)
BL_indi_sgRNA = read.csv(file = "res_lrt_BL_merge_CPM15.csv", header = TRUE)

col_BW = c("sgRNA", "Gene","logFC_BW","logCPM_BW", "F_BW", "PValue_BW")
col_BL = c("sgRNA", "Gene","logFC_BL","logCPM_BL", "F_BL", "PValue_BL")
colnames(BW_indi_sgRNA) =col_BW
colnames(BL_indi_sgRNA) =col_BL

BW_rand <- filter(BW_indi_sgRNA, Gene == "randomized")
# this set makes a scatter plot of all randomized sgRNAs from their logFC 
p2 = ggplot(data = BW_rand , aes(x = sgRNA, y = logFC_BW)) +
  geom_point(shape = 20, alpha = 1) +
  labs(title = "Fold change of all randomized sgRNA in BW25113",
       x = "Randomized sgRNA",
       y = "Log fold change") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p2
pdf(file = "scattteredplot_randsgRNA_BW_CPM15.pdf",width = 8, height = 6)
p2
dev.off()

BL_rand <- filter(BL_indi_sgRNA, Gene == "randomized")
p2 = ggplot(data = BL_rand , aes(x = sgRNA, y = logFC_BL)) +
  geom_point(shape = 20, alpha = 1) +
  labs(title = "Fold change of all randomized sgRNA in BL21",
       x = "Randomized sgRNA",
       y = "Log fold change") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p2
pdf(file = "scatteredplot_randsgRNA_BL_CPM15.pdf",width = 8, height = 6)
p2
dev.off()

# this is to calculate the log of mean FC and log of SD FC 
# then we are plotting the histogram distribution and draw the line for mean and SD on the plotted histogram

BW_rand$FC_BW <- 2^(BW_rand$logFC_BW) # this is to transform logFC to FC
BL_rand$FC_BL <- 2^(BL_rand$logFC_BL)

# this section is to calculate meanFC of all 1000 sgRNA and log-transform the FC to log2(meanFC)
# We calculate mean of FC because it is a real effect not a log-transformed data
# the mean and SD will be log-transformed as a reference for determining significance of the gene

mean_BW = mean(x = BW_rand$FC_BW, na.rm = TRUE) 
log_mean_BW = round(log(mean_BW, 2), digits = 2)
sd_BW = sd(x = BW_rand$FC_BW, na.rm=TRUE)
log_sd_BW = round(log(sd_BW,2), digits = 2)
  
mean_BL = mean(x = BL_rand$FC_BL, na.rm = TRUE)
log_mean_BL = round(log(mean_BL, 2), digits = 2)
sd_BL = sd(x = BL_rand$FC_BL, na.rm=TRUE)
log_sd_BL = round(log(sd_BL,2), digits = 2)

# this section is to make histogram for randomized distribution 
# I decided to plot 
#BW25113
p1 = ggplot(data = BW_rand, aes(x = FC_BW)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity", color = "black")+
  labs(title = "Distribution of FC of randomized sgRNA in BW25113 selection (CPM15)",
       x = "Fold Change",
       y = "Number of sgRNA count") +
  geom_vline(xintercept = c(mean_BW-sd_BW, mean_BW+sd_BW), linetype = "dashed", color = "#488286")+
  geom_vline(xintercept = mean_BW, linetype = "solid", color = "pink")+
  geom_vline(xintercept = 1)+
  annotate("text", x = 3, y=40, label = paste0("Mean FC = ",round(mean_BW, digits = 2))) +
  annotate("text", x = 3, y=35, label = paste0("SD FC = ", round(sd_BW, digits = 2)))+
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p1
pdf(file = "20231228_dist_randomized_BW_CPM15_unpaired_FC.pdf",width = 8, height = 6)
p1
dev.off()

#BL21
p5 = ggplot(data = BL_rand, aes(x = FC_BL)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity", color = "black")+
  labs(title = "Distribution of FC of randomized sgRNA in BL21 selection (CPM15)",
       x = "Fold Change",
       y = "Number of sgRNA count") +
  geom_vline(xintercept = c(mean_BL), linetype = "solid", color = "pink")+
  geom_vline(xintercept = c(mean_BL-sd_BL, mean_BL+sd_BL), linetype = "dashed", color = "#488286")+
  geom_vline(xintercept = 1)+
  annotate("text", x = 3, y=60, label = paste0("Mean FC = ",round(mean_BL, digits = 2))) +
  annotate("text", x = 3, y=55, label = paste0("SD FC = ",round(sd_BL, digits = 2)))+
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p5
pdf(file = "20231228_unpaired_dist_randomized_BL_FC.pdf",width = 8, height = 6)
p5
dev.off()

# this section is to take the final combined p-value and mean FC of each gene and then filtering by mean and SD of randomized data
function_finaltbl = function(lrtcountable, hostfull, hostshort, adj_method, threshold, sd_rand, mean_rand){
  list_gene = unique(sgRNAmatched['Gene']) # to get the list of all genes uniquely
  main_table = matrix(c('Gene','PValue'), ncol=2) # to make new table
  j =0
  for( i in 1:nrow(list_gene)){
    print(list_gene[i,]) # to print each gene out when it's processing
    filtertbl = filter(lrtcountable, Gene == list_gene[i,])
    filtertbl$adjPValue = p.adjust(filtertbl$PValue, method = adj_method) #to adjust p-value by FDR method
    combinedPvaluetbl = sumz(filtertbl$adjPValue, na.action = na.omit) #to combine p-value using Stouffer's method
    main_table <- rbind(main_table, list(list_gene[i,], combinedPvaluetbl$p)) # to write data on the table
  }
  combinetbl = main_table
  colnames(combinetbl) = combinetbl[1,]
  combinetbl = combinetbl[-1,]
  filenametbl = paste0(hostshort,"_combinedpvalue.csv")
  write.csv(combinetbl, file = filenametbl) # to save the file
  tblnew = read.csv(file = filenametbl, header = TRUE) # read the saved file
  
  lrtprocess = lrtcountable
  lrtprocess$FC = 2^(lrtprocess$logFC) # to transform logFC of each sgRNA to Fold Change (FC)
  sumFC <- lrtprocess %>%  group_by(Gene) %>%  summarise_at(vars(FC), list(avg_FC = mean, SD_FC = sd), na.rm = TRUE) # to calculate mean of FC and SD of FC
  sumFC$logFC = log(sumFC$avg_FC,2) #after we get mean of FC --> we log-transformed the mean FC to log of meanFC
  
  tblnew = tblnew[,-1]
  tbl_merge = merge(tblnew, sumFC, by = "Gene") # to merge the 2 tables by Gene 
  tbl_merge$Expression = ifelse(tbl_merge$logFC >= (mean_rand+sd_rand) & tbl_merge$PValue <= threshold, "Up/More Phage",
                                ifelse(tbl_merge$logFC  <= (mean_rand-sd_rand) & tbl_merge$PValue <= threshold, "Down/Less Phage", "Unchanged"))
  # assigning the category to each gene by referring to log of meanFC and log of SDFC, dn threshold (p-value)
  tbl_merge$neglogPValue = -log(tbl_merge$PValue, 10) # confident calculation
  finalresname = paste0(hostshort,"_combinedpvaluelogFC_",threshold,".csv")
  write.csv(tbl_merge, file = finalresname)
  plottitlename = paste0(hostfull," phage selection p-value threshold = ", threshold)
  
  p1 = ggplot(tbl_merge, aes(x=logFC, y=neglogPValue)) +
    geom_point(aes(color = Expression), size = 1) + # color by "Expression" column
    xlab(expression("log"[2]*"FC")) +
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("#E8175D", "#BEBEBE", "#2F9599")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_hline(yintercept = -log(threshold,10), linetype = "dashed") +
    geom_vline(xintercept = c(mean_rand+sd_rand, mean_rand-sd_rand), linetype = "dashed",color = "#488286") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey" ) +
    geom_vline(xintercept = c(mean_rand), linetype = "solid", color = "pink" ) +
    ggtitle(label = plottitlename) +
    geom_text_repel(data = tbl_merge[which(tbl_merge$Expression %in% c("Up/More Phage","Down/Less Phage" )),],
                    aes(x=logFC, y=neglogPValue, label = Gene), size = 3,
                    min.segment.length = 0.1,
                    na.rm = TRUE,
                    show.legend = FALSE)+
    theme_bw()+
    theme(panel.grid.major=element_blank())+
    theme(panel.grid.minor=element_blank())+
    theme(legend.key=element_blank())+
    theme(legend.position="none")
  # ylim(0,10) #add this because I want to expand the graph data for BW25113
  # ylim(0,40) #add this for y-axis limit for BL21
  return(p1)
}
#lrtcountable: the lrt table of logFC generated from edgeR for each sgRNA in the library (in csv)
#hostfull: this is the fullname of the host (ex. "BW25113")
#hostshort: this is the shortname of the host (ex. "BW")
#adj_method: enter the adjustment method for the p-value
#threshold: the p-value you want to use
#sd_rand: the log of SD of fold change of randomized sgRNA data
#mean_rand: the log of mean of fold change of randomized sgRNA data

testrun_BW = function_finaltbl(
  lrtcountable = res_lrt_BW_merge, 
  hostfull = "BW25113",
  hostshort = "BW",
  adj_method = "fdr",
  threshold = 0.05,
  sd_rand = abs(log(sd_BW,2)), #it has to be an absolute number of SD to be used as a reference
  mean_rand = log(mean_BW,2)) # this is log-transformed

pdf(file = "20231228_BW_CPM15_filterRandom_attheEnd_redo.pdf", width = 9, height = 6.5)
print(testrun_BW)
dev.off()

testrun_BL = function_finaltbl(
  lrtcountable = res_lrt_BL_merge,
  hostfull = "BL21",
  hostshort = "BL",
  adj_method = "fdr",
  threshold = 0.05,
  sd_rand = abs(log(sd_BL,2)),
  mean_rand = log(mean_BL,2))

pdf(file = "20231228_BL_CPM15_filterRandom_attheEnd_redo_new.pdf", width = 9, height = 6.5)
print(testrun_BL)
dev.off()

# Here is for O121
O121 = read_excel(path = "/Users/chutikarnchitboonthavisuk/Desktop/Current/Sequencing_Work/20230509_clinicalstrain/20230511_O121_data.xlsx", sheet = 1)

# genelist --------
O121list = O121[,1:2]
tbl1 = O121
x = tbl1[, grep(("sgRNA|Gene|O121"), names(tbl1))]
colnames(x) = c("sgRNA","Gene","O121_noPhage_45m_1","O121_noPhage_2hr_1","O121_PhageMOI1_45m_1", 
                "O121_PhageMOI1_2hr_1", "O121_PhageMOI2_45m_1", "O121_PhageMOI2_2hr_1","O121_noPhage_45m_2","O121_noPhage_45m_3",
                "O121_PhageMOI2_45m_2", "O121_PhageMOI2_45m_3", "O121_noPhage_2hr_2", "O121_noPhage_2hr_3", "O121_PhageMOI2_2hr_2",
                "O121_PhageMOI2_2hr_3")
## for sanity check ----------
x_mod = gather(x, "sample", "rawReads", 3:16)
sum = aggregate(x_mod$rawReads, list(x_mod$sample), FUN=sum)
colnames(sum) = c("sample", "total")
x_mod1 = merge(x_mod, sum, by.x = "sample", all = TRUE)
col_order = c("sample", "sgRNA", "Gene", "rawReads", "total")
x_mod1 <- x_mod1[,col_order]
x_mod1$cutoff3 = ifelse(x_mod1$rawReads < 3, "reads < 3", "reads >= 3")
x_filtered = x_mod1
x_filtered$normreads = x_filtered$rawReads/x_filtered$total
x_wide = x_filtered[, which(colnames(x_filtered) %in% c("sgRNA","Gene","sample","normreads"))]
x_wide =pivot_wider(x_wide, names_from = sample, values_from = normreads)

x_cal = x_wide
x_cal$type = ifelse(x_cal$Gene == "randomized", "Randomized","Targeted")

# EdgeR analysis for 45 minutes 
tbl1 = O121
x = tbl1[, grep(("sgRNA|45m"), names(tbl1))]
colnames(x) = c("sgRNA", "O121_noPhage_45m_1", "O121_PhageMOI1_45m_1", "O121_PhageMOI2_45m_1",
                "O121_noPhage_45m_2", "O121_noPhage_45m_3","O121_PhageMOI2_45m_2", "O121_PhageMOI2_45m_3")
col_order <- c("sgRNA","O121_noPhage_45m_1", "O121_noPhage_45m_2", "O121_noPhage_45m_3",
               "O121_PhageMOI1_45m_1", "O121_PhageMOI2_45m_1", "O121_PhageMOI2_45m_2", "O121_PhageMOI2_45m_3")
count_x = x[, col_order]
count_x = count_x[, grep(("sgRNA|noPhage|PhageMOI2"), names(count_x))]
  
row.names(count_x) = count_x$sgRNA
count_x = as.matrix(count_x)
count_x = count_x[,2:7]
# group = c(rep("control",3 ), rep("treatment", 3))
count_x = as.data.frame(count_x)
count_x$O121_noPhage_45m_1 <- as.numeric(count_x$O121_noPhage_45m_1)
count_x$O121_noPhage_45m_2 <- as.numeric(count_x$O121_noPhage_45m_2)
count_x$O121_noPhage_45m_3 <- as.numeric(count_x$O121_noPhage_45m_3)
count_x$O121_PhageMOI2_45m_1 <- as.numeric(count_x$O121_PhageMOI2_45m_1)
count_x$O121_PhageMOI2_45m_2 <- as.numeric(count_x$O121_PhageMOI2_45m_2)
count_x$O121_PhageMOI2_45m_3 <- as.numeric(count_x$O121_PhageMOI2_45m_3)

group = c(rep("control", 3), rep("treatment", 3))
DGE_list = DGEList(count = count_x, group = group)
keep <- rowSums(cpm(DGE_list$count[,1:3])>3 ) >=3 # this is for CPM = 3
DGE_list <- DGE_list[keep, , keep.lib.size = FALSE]

### making design matrix 
designtbl <- data.frame()
sample <- colnames(count_x)
subject <- factor(c(1,2,3,1,2,3)) # this is a pairwise 
treatment <- factor(c("C","C","C","T","T","T"))
designtbl <- data.frame(sample, subject, treatment)

Subject <- factor(designtbl$subject)
Treat <- factor(designtbl$treatment, levels = c("C","T"))
design <- model.matrix(~Subject+Treat)

DGE_list <- estimateDisp(DGE_list, design)
fit <- glmQLFit(DGE_list, design)
lrt <- glmQLFTest(fit)

test_lrt = lrt$table
test_lrt$sgRNA = row.names(test_lrt)
test_lrt_merge = merge(O121list, test_lrt,by = "sgRNA", all.x = TRUE)
write.csv(x = test_lrt_merge, file = "45min_lrt_CPM3_paired.csv", row.names = FALSE)
lrt_45min = test_lrt_merge # this is for lrt 45 minutes

tbl1 = O121
x = tbl1[, grep(("sgRNA|2h"), names(tbl1))]
colnames(x) = c("sgRNA", "O121_noPhage_2hr_1", "O121_PhageMOI1_2hr_1", "O121_PhageMOI2_2hr_1",
                "O121_noPhage_2hr_2", "O121_noPhage_2hr_3","O121_PhageMOI2_2hr_2", "O121_PhageMOI2_2hr_3")
col_order <- c("sgRNA","O121_noPhage_2hr_1", "O121_noPhage_2hr_2", "O121_noPhage_2hr_3",
               "O121_PhageMOI1_2hr_1", "O121_PhageMOI2_2hr_1", "O121_PhageMOI2_2hr_2", "O121_PhageMOI2_2hr_3")
count_x = x[, col_order]
count_x = count_x[, grep(("sgRNA|noPhage|PhageMOI2"), names(count_x))]

row.names(count_x) = count_x$sgRNA
count_x = as.matrix(count_x)
count_x = count_x[,2:7]
# group = c(rep("control",3 ), rep("treatment", 3))
count_x = as.data.frame(count_x)
count_x$O121_noPhage_2hr_1 <- as.numeric(count_x$O121_noPhage_2hr_1)
count_x$O121_noPhage_2hr_2 <- as.numeric(count_x$O121_noPhage_2hr_2)
count_x$O121_noPhage_2hr_3 <- as.numeric(count_x$O121_noPhage_2hr_3)
count_x$O121_PhageMOI2_2hr_1 <- as.numeric(count_x$O121_PhageMOI2_2hr_1)
count_x$O121_PhageMOI2_2hr_2 <- as.numeric(count_x$O121_PhageMOI2_2hr_2)
count_x$O121_PhageMOI2_2hr_3 <- as.numeric(count_x$O121_PhageMOI2_2hr_3)

group = c(rep("control", 3), rep("treatment", 3))
DGE_list = DGEList(count = count_x, group = group)
keep <- rowSums(cpm(DGE_list$count[,1:3])>3 ) >=3 # I will try lower CPM; this step should be done before anything else
DGE_list <- DGE_list[keep, , keep.lib.size = FALSE]

### making design matrix 
designtbl <- data.frame()
sample <- colnames(count_x)
subject <- factor(c(1,2,3,1,2,3))
treatment <- factor(c("C","C","C","T","T","T"))
designtbl <- data.frame(sample, subject, treatment)

Subject <- factor(designtbl$subject)
Treat <- factor(designtbl$treatment, levels = c("C","T"))
design <- model.matrix(~Subject+Treat)

DGE_list <- estimateDisp(DGE_list, design) # at this step, you can adjust pseudocount by using "prior.count = ___", the default pseudocount is 0.5
fit <- glmQLFit(DGE_list, design)
lrt <- glmQLFTest(fit)

test_lrt = lrt$table
test_lrt$sgRNA = row.names(test_lrt)
test_lrt_merge = merge(O121list, test_lrt,by = "sgRNA", all.x = TRUE)
write.csv(x = test_lrt_merge, file = "lrt_2hr_CPM3_paired.csv", row.names = FALSE)

#import lrt 
lrt_45min <- read.csv(file = "45min_lrt_CPM3_paired.csv", header= TRUE)
lrt_2hr <- read.csv(file = "lrt_2hr_CPM3_paired.csv", header = TRUE)
col_45min = c("sgRNA", "Gene","logFC_45min","logCPM_45min", "F_45min", "PValue_45min")
col_2hr = c("sgRNA", "Gene","logFC_2hr","logCPM_2hr", "F_2hr", "PValue_2hr")
colnames(lrt_45min) =col_45min
colnames(lrt_2hr) =col_2hr

## 45-min ------

lrt_45min_mod = lrt_45min
lrt_45min_mod$FC_45min = 2^(lrt_45min_mod$logFC_45min)
lrt_45min_rand = filter(lrt_45min_mod, Gene == "randomized")

mean_45min = round(mean(lrt_45min_rand$FC_45min, na.rm = TRUE), digits = 2)
log_mean45min = log(mean_45min,2)
SD_45min = round(sd(lrt_45min_rand$FC_45min, na.rm = TRUE), digits = 2)
log_SD45min = log(SD_45min, 2)


p1 = ggplot(data = lrt_45min_rand, aes(x = FC_45min)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  labs(title = "Distribution of FC of randomized sgRNA in 45-minutes selection",
       x = "Fold Change",
       y = "Number of sgRNA count") +
  
  geom_vline(xintercept = c(mean_45min), linetype = "solid", color = "pink")+
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(mean_45min+SD_45min, mean_45min-SD_45min), linetype = "dashed",color = "#488286") +
  annotate("text", x = 3, y=200, label = paste0("Mean = ",mean_45min)) +
  annotate("text", x = 3, y=150, label = paste0("SD = ", SD_45min)) +
  annotate("text", x = 3, y=100, label = paste0("Depleted sgRNA = ",349)) +
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p1

pdf(file = "20231228_FC_randomdist_45min.pdf",width = 8, height = 6)
p1
dev.off()

# this is to plot the log2FC distribution 
logmean_45 = log(mean_45min,2)
logSD_45 = log(SD_45min,2)
p1 = ggplot(data = lrt_45min_rand, aes(x = logFC_45min)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  labs(title = "Distribution of log2FC of randomized sgRNA in 45-minutes selection",
       x = "Log2 Fold Change",
       y = "Number of sgRNA count") +
  geom_vline(xintercept = logmean_45, linetype = "solid", color = "pink")+
  # geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(logmean_45+logSD_45, logmean_45-logSD_45), linetype = "dashed",color = "#488286") +
  annotate("text", x = 0, y=10, label = paste0("Mean = ",round(logmean_45, digits = 2))) +
  annotate("text", x = 0, y=9, label = paste0("SD = ",round(logSD_45, digits = 2))) +
  annotate("text", x = 0, y=8, label = paste0("Depleted sgRNA = ",349)) +
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p1

pdf(file = "log_meanFC_dist_45min.pdf",width = 8, height = 6)
p1
dev.off()

## 2-hr ------------
lrt_2hr_mod = lrt_2hr
lrt_2hr_rand = filter(lrt_2hr_mod, Gene == "randomized")
lrt_2hr_rand$FC_2hr = 2^lrt_2hr_rand$logFC_2hr

mean_2hr = round(mean(lrt_2hr_rand$FC_2hr, na.rm = TRUE), digits = 2)
SD_2hr = round(sd(lrt_2hr_rand$FC_2hr, na.rm = TRUE), digits = 2)

p2 = ggplot(data = lrt_2hr_rand , aes(x=FC_2hr)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  labs(title = "Distribution of FC of randomized sgRNA in 2-hours selection",
       x = "Log2 Fold Change",
       y = "Number of sgRNA count") +
  geom_vline(xintercept = c(mean_2hr), linetype = "solid", color = "pink")+
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(mean_2hr+SD_2hr, mean_2hr-SD_2hr), linetype = "dashed",color = "#488286") +
  annotate("text", x = 3, y=200, label = paste0("Mean = ", mean_2hr)) +
  annotate("text", x = 3, y=150, label = paste0("SD = ",SD_2hr)) +
  annotate("text", x = 3, y=100, label = paste0("Depleted sgRNA = ",357)) +
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p2
pdf(file = "20231228_FCdist_randomsgRNA_2hr.pdf",width = 8, height = 6)
p2
dev.off()

logmean_2 = log(mean_2hr,2)
logSD_2 = log(SD_2hr,2)

p2 = ggplot(data = lrt_2hr_rand , aes(x=logFC_2hr)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  labs(title = "Distribution of log2FC of randomized sgRNA in 2-hours selection",
       x = "Log2 Fold Change",
       y = "Number of sgRNA count") +
  geom_vline(xintercept = logmean_2, linetype = "solid", color = "pink")+
  geom_vline(xintercept = c(logmean_2+logSD_2, logmean_2-logSD_2), linetype = "dashed",color = "#488286") +
  annotate("text", x = 0, y=10, label = paste0("Mean = ",round(logmean_2, digits = 2))) +
  annotate("text", x = 0, y=9, label = paste0("SD = ",round(logSD_2, digits = 2))) +
  annotate("text", x = 0, y=8, label = paste0("Depleted sgRNA = ",357)) +
  theme_bw() +
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(strip.background=element_blank())
p2
pdf(file = "20231228_logmeanFCdist_randomsgRNA_2hr.pdf",width = 8, height = 6)
p2
dev.off()

# making volcano plot 

function_finaltbl = function(lrtcountable, hostfull, hostshort, adj_method, threshold, mean_rand, sd_rand){
  list_gene = unique(O121list['Gene'])
  list_gene = as.data.frame(list_gene)
  main_table = matrix(c('Gene','PValue'), ncol=2)
  j =0
  for( i in 1:nrow(list_gene)){
    print(list_gene[i,])
    filtertbl = filter(lrtcountable, Gene == list_gene[i,])
    filtertbl$adjPValue = p.adjust(filtertbl$PValue, method = adj_method)
    combinedPvaluetbl = sumz(filtertbl$adjPValue, na.action = na.omit)
    main_table <- rbind(main_table, list(list_gene[i,], combinedPvaluetbl$p))
  }
  combinetbl = main_table
  colnames(combinetbl) = combinetbl[1,]
  combinetbl = combinetbl[-1,]
  filenametbl = paste0(hostshort,"_combinedpvalue.csv")
  write.csv(combinetbl, file = filenametbl)
  tblnew = read.csv(file = filenametbl, header = TRUE)
  
  lrtprocess = lrtcountable
  lrtprocess$FC = 2^(lrtprocess$logFC)
  sumFC <- lrtprocess %>%  group_by(Gene) %>%  summarise_at(vars(FC), list(avg_FC = mean, SD_FC = sd), na.rm = TRUE)
  sumFC$logFC = log(sumFC$avg_FC,2)
  
  tblnew = tblnew[,-1]
  tbl_merge = merge(tblnew, sumFC, by = "Gene")
  tbl_merge$Expression = ifelse(tbl_merge$logFC > (mean_rand + abs(sd_rand)) & tbl_merge$PValue <= threshold, "Up/More Phage",
                                ifelse(tbl_merge$logFC <  (mean_rand - abs(sd_rand)) & tbl_merge$PValue <= threshold, "Down/Less Phage", "Unchanged"))
  tbl_merge$neglogPValue = -log(tbl_merge$PValue, 10)
  
  finalresname = paste0(hostshort,"_combinedpvaluelogFC_",threshold,".csv")
  write.csv(tbl_merge, file = finalresname)
  plottitlename = paste0(hostfull," phage selection p-value threshold = ", threshold)
  
  p1 = ggplot(tbl_merge, aes(x=logFC, y=neglogPValue)) +
    geom_point(aes(color = Expression), size = 1) +
    xlab(expression("log"[2]*"FC")) +
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("#BEBEBE", "#BEBEBE", "#2F9599")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_hline(yintercept = -log(threshold,10), linetype = "dashed") +
    geom_vline(xintercept = c(mean_rand, mean_rand + sd_rand, mean_rand - sd_rand), linetype = "dashed") +
    ggtitle(label = plottitlename) +
    geom_text_repel(data = tbl_merge[which(tbl_merge$Expression %in% c("Up/More Phage")),],
                    aes(x=logFC, y=neglogPValue, label = Gene), size = 3,
                    min.segment.length = 0.1,
                    na.rm = TRUE,
                    show.legend = FALSE)+
    scale_y_continuous(trans='log10')+
    theme_bw()+
    theme(panel.grid.major=element_blank())+
    theme(panel.grid.minor=element_blank())+
    theme(legend.key=element_blank())+
    theme(legend.position="none")
  return(p1)
}
O121_run = function_finaltbl(lrtcountable = lrt_45min, hostfull = "O121",
                             hostshort = "O121", adj_method = "fdr", threshold = 0.05,
                             mean_rand = logmean_45,sd_rand = logSD_45)

O121_run
pdf(file = "20231212_test_big_logscale.pdf", width = 20, height = 14.5)
print(O121_run)
dev.off()

# noted that when you run this, you can remove scale transformation log
O121_2hr = function_finaltbl(lrtcountable = lrt_2hr, hostfull = "O121",
                             hostshort = "O121", adj_method = "fdr", threshold = 0.05,
                             mean_rand = logmean_2,sd_rand = logSD_2)
O121_2hr


