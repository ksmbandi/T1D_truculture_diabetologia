#trying to get all comparisons. this might not work. refer otherwise to original.
#original analysis is T1D_TruCulture_analysis_for_Kameron2 filtered2 combat corrected for quant batch.r

## scripts to asses variation in T1D TruCulture nanostring data

#install packages missing
# install.packages("ordinal")
# install.packages("ggthemes")
# install.packages("devtools")
# source("https://bioconductor.org/biocLite.R")
# biocLite("scran")
# biocLite("ggrepel")
# library(devtools)
# install_github("mjdufort/countSubsetNorm")
# install_github("mjdufort/RNAseQC")
# install_github("mjdufort/limmaTools")
# install_github("mjdufort/miscHelpers")
# install_github("mjdufort/geneSetTools")

#does removing HC16 LPS change the combat batch normalization a lot?
#does removing HC16 LPS change the significant DE genes a lot?
library(ggrepel)

##### set up environment: load packages #####

# load general packages
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1)))
library(gplots)
library(ordinal)
library(ggthemes)

# load analysis_specific packages
library(edgeR)
library(limma)

# load my relevant functions (at github.com/mjdufort)
library(countSubsetNorm)
library(RNAseQC)
library(limmaTools)
library(miscHelpers)
library(geneSetTools)

# load magrittr
library(magrittr)


##### read in annotation and TruCulture nanostring data #####

setwd("../Data shared with BRI TruCulture_T1D_innate_stim_project")
# setwd("./Data shared with BRI TruCulture_T1D_innate_stim_project")

# read in raw data

#clinical_nanostring_data.T1D <-
#  readxl::read_xlsx(
#    "truculture_nanostring_tarbell_data_T1D_filtered_final2.xlsx",
#    sheet=5) %>%
#colnames(clinical_nanostring_data.T1D)[
#  1:which(colnames(clinical_nanostring_data.T1D)=="stimulus")] <-
#  colnames(clinical_nanostring_data.T1D)[
#    1:which(colnames(clinical_nanostring_data.T1D)=="stimulus")] %>%
#  standardize_names()



# read in normalized nanostring data
clinical_nanostring_data.T1D <-
  readxl::read_xlsx(
    "truculture_nanostring_tarbell_data_T1D_filtered_final2.xlsx",
    sheet=1) %>%
  plyr::rename(
    c("PatientType"="group",
      "StimulusId"="stimulus",
      "Gender"="sex",
      "extraction batch"="extract_batch",
      "quant machine"="quant_batch"))
colnames(clinical_nanostring_data.T1D)[
  1:which(colnames(clinical_nanostring_data.T1D)=="stimulus")] <-
  colnames(clinical_nanostring_data.T1D)[
    1:which(colnames(clinical_nanostring_data.T1D)=="stimulus")] %>%
  standardize_names()

# standardize annotation column names and add sample_id
clinical_nanostring_data.T1D <-
  clinical_nanostring_data.T1D %>%
  tibble::add_column(
    sample_id =
      with(clinical_nanostring_data.T1D,
           paste(patientid, stimulus, sep="_")),
    .after="stimulus") %>%
  as.data.frame()
glimpse(clinical_nanostring_data.T1D)


# extract count data
nanostring.norm_counts.T1D <-
  clinical_nanostring_data.T1D[
    , (which(colnames(clinical_nanostring_data.T1D)=="sample_id") + 1):
      ncol(clinical_nanostring_data.T1D)] %>%
  t() %>%
  set_colnames(clinical_nanostring_data.T1D$sample_id)


##### generate "master" annotation object for compatibility with downstream code #####

master <-
  clinical_nanostring_data.T1D[
    , 1:which(colnames(clinical_nanostring_data.T1D)=="sample_id")]
master$group <-
  factor(master$group, levels=c("HC", "T1D"))
master$stimulus <- 
  factor(master$stimulus, levels=c("Null", "IFNb", "IL1", "LPS", "PolyIC", "SEB"))
master$sex <-
  factor(master$sex, levels = c("male","female"))
master$extract_batch <-
  factor(master$extract_batch, levels = c("1","2","3"))
master$quant_batch <-
  factor(master$quant_batch, levels = c("qubit","tecan"))

rm(data)
#'data' object
load("nanostring data combat quant batch corrected nonlog2.rda")
saveRDS(object = data,file = "Nanostring data combat quant batch corrected nonlog2.rda")
saveRDS(object = master, file = "master.rda")

#data2<- data #this is everything
#data<- subset(x=data2, select = -(HC16_LPS))#this is without HC16 LPS sample (but it was included for ComBat batch effect removal) 

#calculate ranks for each column and make dataRanks dataframe



#'new_data' object
#load("nanostring_data_calcNormFactors quant batch effect corrected unlogged.rda")

# generate vwts object for general use in plotting
# rm(vwts.all)
if (!exists("vwts.all")) {
  vwts.all <-
    #nanostring.norm_counts.T1D
    data  %>%
    #new_data 
    calc_norm_counts(
      design=master, libID_col="sample_id",
      min_count = 1, min_libs_perc = 0.15,
      log2_transform=FALSE, normalize=FALSE) %>%
    voomWithQualityWeights(plot=TRUE, span=0.2) # reduce span to better model low-count genes
}

# output normalized, log-transformed counts for use outside R (sharing, etc.)
vwts.all$E %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  write.csv(
    "data_output/counts.log2_cpm_plus_0.5_filtered2_combat_quant_batch_corrected.csv",
    row.names=FALSE)


##### Set up and run limma for patient group and stimulus, with interaction #####

condition.tmp <- 
  with(master, !is.na(group) & !is.na(stimulus))
master.group.stimulus <-
  master %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.group.stimulus.tmp <-
  calc_norm_counts(
    counts= data, #nanostring.norm_counts.T1D,
    design=master.group.stimulus, libID_col="sample_id",
    min_cpm = 1, min_libs_perc = 0.15,
    #normalize=TRUE, log2_transform = TRUE)
    normalize=FALSE, log2_transform=FALSE)
master.group.stimulus <-
  master.group.stimulus[
    match(colnames(DGECounts.group.stimulus.tmp),
          master.group.stimulus[,"sample_id"]),]

# rm(vwts.group.stimulus.with_interaction)
if (!exists("vwts.group.stimulus.with_interaction")) {
  DesignMat.group.stimulus.with_interaction <-
    model.matrix(
      # ~ group + stimulus + group:stimulus + sex + quant_batch,
      # ~ group + stimulus + group:stimulus + sex + extract_batch + quant_batch,
      #~ group + stimulus + group:stimulus + quant_batch,
      
      ~ group + stimulus + group:stimulus, #original and correct
      #~0+group+stimulus+group:stimulus,
      data=master.group.stimulus)
  vwts.group.stimulus.with_interaction <-
    voomWithQualityWeights(
      DGECounts.group.stimulus.tmp,
      design=DesignMat.group.stimulus.with_interaction,
      plot=TRUE, span=0.2) # reduce span to better model low-count genes
  corfit.group.stimulus.with_interaction <-
    duplicateCorrelation(
      vwts.group.stimulus.with_interaction,
      design=DesignMat.group.stimulus.with_interaction,
      block=master.group.stimulus$patientid)
  vfit.group.stimulus.with_interaction <-
    lmFit(vwts.group.stimulus.with_interaction,
          block=master.group.stimulus$patientid,
          correlation=corfit.group.stimulus.with_interaction$consensus.correlation) %>%
    eBayes()
  topGenes.group.stimulus.with_interaction <-
    vfit.group.stimulus.with_interaction %>%
    topTable(
      coef = 2:ncol(.$coefficients),
      number=Inf, sort.by="F")
}
saveRDS(object = vfit.group.stimulus.with_interaction,file = "vfit_group_stimulus_with_interaction.rda")

# what are the most DE genes?
head(topGenes.group.stimulus.with_interaction[
  order(topGenes.group.stimulus.with_interaction$P.Value),], 20)
# head(topGenes.group.stimulus.with_interaction[
#   order(abs(topGenes.group.stimulus.with_interaction$logFC), decreasing=TRUE),], 20)


##### extract results for the difference in each stimulus effect between T1D and HC #####

## extract topGenes object for each stimulus T1D interaction effect
topGenes.group.stimulus.with_interaction.interaction_by_stimulus <- list()
for (i in (levels(master$stimulus)[-1])) {
  topGenes.group.stimulus.with_interaction.interaction_by_stimulus[[i]] <-
    topTable(
      vfit.group.stimulus.with_interaction,
      coef = paste0("groupT1D:stimulus", i), 
      number=Inf, sort.by="P") %>%
    tibble::rownames_to_column(var="gene") %>%
    tibble::add_column(stimulus=i, .before="gene") %>%
    set_rownames(.[,"gene"])
}

for (i in names(topGenes.group.stimulus.with_interaction.interaction_by_stimulus))
  write.csv(topGenes.group.stimulus.with_interaction.interaction_by_stimulus[[i]],
            file=paste0("data_output/topGenes.group.stimulus.combat_quant_batch_corrected_filtered2_interaction_by_", i, ".csv"),
            quote=FALSE)


nanostringTableForPublication <- data.frame(vwts.group.stimulus.with_interaction$E,stringsAsFactors = FALSE)
nanostringTableForPublication$IFNbSignatureGenes <- rownames(nanostringTableForPublication) %in% IFNbGenes
nanostringTableForPublication$IFNgSignatureGenes <- rownames(nanostringTableForPublication) %in% IFNgGenes
nanostringTableForPublication$IL1bSignatureGenes <- rownames(nanostringTableForPublication) %in% IL1Bgenes
nanostringTableForPublication$TNFaSignatureGenes <- rownames(nanostringTableForPublication) %in% TNFAgenes
write.table(x = nanostringTableForPublication, file = "./nanostring data for publication.csv",sep = ",",col.names = NA )


# "(Intercept)", 
# "groupT1D", 
# "stimulusIFNb", 
# "stimulusIL1", 
# "stimulusLPS", 
# "stimulusPolyIC",
# "stimulusSEB", 
# "groupT1D:stimulusIFNb", 
# "groupT1D:stimulusIL1",
# "groupT1D:stimulusLPS", 
# "groupT1D:stimulusPolyIC", 
# "groupT1D:stimulusSEB"

contrasts.T1D_vs_HC <-
  cbind( T1D_Null_vs_HC_Null =     c(0,1,0,0,0,0,0,0,0,0,0,0),
         T1D_IFNb_vs_HC_IFNb =     c(0,1,0,0,0,0,0,1,0,0,0,0),
         T1D_PolyIC_vs_HC_PolyIC = c(0,1,0,0,0,0,0,0,0,0,1,0),
         T1D_LPS_vs_HC_LPS =       c(0,1,0,0,0,0,0,0,0,1,0,0),
         T1D_IL1b_vs_HC_IL1b =     c(0,1,0,0,0,0,0,0,1,0,0,0),
         T1D_SEB_vs_HC_SEB =       c(0,1,0,0,0,0,0,0,0,0,0,1),
         IFNb_vs_Null = c(0,0,1,0,0,0,0,0.5,0,0,0,0),
         PolyIC_vs_Null = c(0,0,0,0,0,1,0,0,0,0,0.5,0),
         LPS_vs_Null = c(0,0,0,0,1,0,0,0,0,0.5,0,0),
         IL1b_vs_Null = c(0,0,0,1,0,0,0,0,0.5,0,0,0),
         SEB_vs_Null = c(0,0,0,0,0,0,1,0,0,0,0,0.5))
         
vfit.T1D_vs_HC.by_stimulation <-
  vfit.group.stimulus.with_interaction %>%
  contrasts.fit(contrasts= contrasts.T1D_vs_HC) %>%
  eBayes()


fc_cut.tmp <- 0
p_cut.tmp <- 0.0549

lapply(
  topGenes.group.stimulus.with_interaction.interaction_by_stimulus,
  function(x) x[x$adj.P.Val < p_cut.tmp,])



Null <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_Null_vs_HC_Null", resort.by = "P", n=Inf)
LPS <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_LPS_vs_HC_LPS", resort.by = "P", n=Inf)
PolyIC <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_PolyIC_vs_HC_PolyIC", resort.by = "P", n=Inf)
IFNb <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_IFNb_vs_HC_IFNb", resort.by = "P", n=Inf)
IL1b <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_IL1b_vs_HC_IL1b", resort.by = "P", n=Inf)
SEB <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="T1D_SEB_vs_HC_SEB", resort.by = "P", n=Inf)

IFNbstim <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="IFNb_vs_Null", resort.by = "P", n=Inf)
PolyICstim <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="PolyIC_vs_Null", resort.by = "P", n=Inf)
LPSstim <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="LPS_vs_Null", resort.by = "P", n=Inf)
IL1bstim <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="IL1b_vs_Null", resort.by = "P", n=Inf)
SEBstim <- topTable(vfit.T1D_vs_HC.by_stimulation, coef="SEB_vs_Null", resort.by = "P", n=Inf)



write.csv(Null,"comparisons after combat/noninteraction comparisons/T1D_Null - HC_Null.csv")
write.csv(LPS,"comparisons after combat/noninteraction comparisons/T1D_LPS - HC_LPS with HC16.csv")
write.csv(PolyIC,"comparisons after combat/noninteraction comparisons/T1D_PolyIC - HC_PolyIC.csv")
write.csv(IFNb,"comparisons after combat/noninteraction comparisons/T1D_IFNb - HC_IFNb.csv")
write.csv(IL1b,"comparisons after combat/noninteraction comparisons/T1D_IL1b - HC_IL1b.csv")
write.csv(SEB,"comparisons after combat/noninteraction comparisons/T1D_SEB - HC_SEB.csv")

write.csv(LPSstim,"comparisons after combat/noninteraction comparisons/LPS - Null.csv")
write.csv(PolyICstim,"comparisons after combat/noninteraction comparisons/PolyIC - Null.csv")
write.csv(IFNbstim,"comparisons after combat/noninteraction comparisons/IFNb - Null.csv")
write.csv(IL1bstim,"comparisons after combat/noninteraction comparisons/IL1b - Null.csv")
write.csv(SEBstim,"comparisons after combat/noninteraction comparisons/SEB - Null.csv")


## volcano plots for each stimulus

# standardize X and Y limits across the conditions
xy_lims.tmp <-
  get_xy_lims(
    topGenes.group.stimulus.with_interaction.interaction_by_stimulus,
    pairwise=TRUE)

for (i in names(topGenes.group.stimulus.with_interaction.interaction_by_stimulus)) {
  plot_volcano_2var(
    topGenes.group.stimulus.with_interaction.interaction_by_stimulus[[i]],
    file_prefix=paste0("plots/TruCulture/volcano.group.stimulus.with_interaction.interaction_combat_quant_batch_corrected_filtered_final_AdjPvalCutoff0.05_", i),
    plotdims=c(6,6),
    x_lim=xy_lims.tmp$x, y_lim=xy_lims.tmp$y,
    fc_cut=fc_cut.tmp, p_cut=p_cut.tmp,
    color_by_threshold=FALSE, my_cols="black",
  #  gene_labs="ellipse", x_cut=2.5, y_cut=1.4, gene_lab_size=5)
     gene_labs="threshold", x_cut=0, y_cut=1.26, gene_lab_size=2, x_cut_direction = "both")
}


plot_volcano_2var(
  Null,
  file_prefix=paste0("plots/TruCulture/volcano.group.stimulus.with_interaction.interaction_combat_quant_batch_corrected_filtered_final_AdjPvalCutoff0.05_", "Null_forFigure"),
  plotdims=c(6,6),
  x_lim=xy_lims.tmp$x, y_lim=xy_lims.tmp$y,
  fc_cut=fc_cut.tmp, p_cut=p_cut.tmp,
  color_by_threshold=FALSE, my_cols="black",
  #  gene_labs="ellipse", x_cut=2.5, y_cut=1.4, gene_lab_size=5)
  gene_labs="threshold", x_cut=0, y_cut=1.26, gene_lab_size=2, x_cut_direction = "both")



rm(plot_volcano_2var)

?geom_text_repel
library(ggrepel)
plot_volcano_2var <- function (topGenes, my_cols = c("darkcyan", "darkorange"), file_prefix = NULL, 
          plotdims = c(9, 9), color_by_threshold = TRUE, fc_cut = log2(1.5), 
          p_cut = 0.01, x_lim = "auto", y_lim = "auto", gene_labs = NULL, 
          x_cut = 0, y_cut = 0, x_cut_direction = "both", gene_labs_repel = TRUE, 
          gene_lab_size = 3, ...) {
  if (color_by_threshold & (is.null(fc_cut) | is.null(p_cut))) 
    stop("Cannot plot points by threshold with null values of fc_cut or p_cut.")
  if (identical(x_lim, "auto") | identical(y_lim, "auto")) {
    xy_lims <- get_xy_lims(topGenes, min_x_abs = fc_cut, 
                           min_y2 = -log10(p_cut))
    if (identical(x_lim, "auto")) 
      x_lim <- xy_lims[["x"]]
    if (identical(y_lim, "auto")) 
      y_lim <- xy_lims[["y"]]
  }
  topGenes$genes <- rownames(topGenes)
  if (color_by_threshold & is.null(topGenes$threshold)) 
    topGenes$threshold <- (abs(topGenes$logFC) > fc_cut) & 
    (topGenes$adj.P.Val < p_cut)
  if (color_by_threshold) {
    volcano <- ggplot(data = topGenes, aes(x = logFC, y = -log10(adj.P.Val), 
                                           colour = threshold)) + geom_point(alpha = 0.6, size = 3, 
                                                                             shape = 16) + scale_colour_manual(values = my_cols)
  }
  else {
    volcano <- ggplot(data = topGenes, aes(x = logFC, y = -log10(adj.P.Val))) + 
      geom_point(alpha = 0.6, size = 2, shape = 16, color = my_cols[1])
  }
  volcano <- volcano + theme(legend.position = "none") + xlab("log2 fold change") + 
    ylab("-log10 Adj P")
  if (!is.null(fc_cut)) {
    volcano <- volcano + geom_vline(xintercept = fc_cut, 
                                    linetype = "dotted", size = 0.5) + geom_vline(xintercept = -fc_cut, 
                                                                                linetype = "dotted", size = 0.5)
  }
  if (!is.null(p_cut)) {
    volcano <- volcano + geom_hline(yintercept = -log10(p_cut), 
                                    linetype = "dotted", size = 0.5)
  }
  if (!is.null(x_lim)) {
    volcano <- volcano + xlim(x_lim)
  }
  if (!is.null(y_lim)) {
    volcano <- volcano + ylim(y_lim)
  }
  if (!is.null(gene_labs)) {
    gene_labs <- match.arg(gene_labs, choices = c("ellipse", 
                                                  "threshold"))
    if (gene_labs == "ellipse") {
      topGenes.tmp <- topGenes[((topGenes$logFC^2)/(x_cut^2) + 
                                  (log10(topGenes$adj.P.Val)^2)/(y_cut^2)) > 1, 
                               ]
    }
    else if (gene_labs == "threshold") {
      x_cut_direction <- match.arg(x_cut_direction, choices = c("both", 
                                                                "lower", "upper"))
      if (x_cut_direction == "both") {
        topGenes.tmp <- topGenes[(abs(topGenes$logFC) > 
                                    x_cut) & (-log10(topGenes$adj.P.Val) > y_cut), 
                                 ]
      }
      else if (x_cut_direction == "lower") {
        topGenes.tmp <- topGenes[(topGenes$logFC < x_cut) & 
                                   (-log10(topGenes$adj.P.Val) > y_cut), ]
      }
      else if (x_cut_direction == "upper") {
        topGenes.tmp <- topGenes[(topGenes$logFC > x_cut) & 
                                   (-log10(topGenes$adj.P.Val) > y_cut), ]
      }
    }
    if (gene_labs_repel) {
      volcano <- volcano + geom_text_repel(data = topGenes.tmp, 
                                           aes(label = genes), color = "black", size = gene_lab_size)
    }
    else {
      volcano <- volcano + geom_text(data = topGenes.tmp, 
                                     aes(label = genes), color = "black", size = gene_lab_size, 
                                     vjust = 1, hjust = 0.5)
    }
  }
  if (!is.null(file_prefix)) {
    pdf(file = paste(file_prefix, "pdf", sep = "."), w = plotdims[1], 
        h = plotdims[2], ...)
    on.exit(dev.off(), add = TRUE)
  }
  else quartz(plotdims[1], plotdims[2])
  print(volcano)
}


plot_volcano_2var <-
  function(topGenes, my_cols=c("darkcyan", "darkorange"),
           file_prefix=NULL, plotdims=c(9,9),
           color_by_threshold=TRUE, fc_cut=log2(1.5), p_cut=0.01,
           x_lim="auto", y_lim="auto",
           gene_labs=NULL, x_cut=0, y_cut=0, x_cut_direction="both",
           gene_labs_repel=TRUE, gene_lab_size=3,
           ...) {
    if (color_by_threshold & (is.null(fc_cut) | is.null(p_cut)))
      stop("Cannot plot points by threshold with null values of fc_cut or p_cut.")
    if (identical(x_lim, "auto") | identical(y_lim, "auto")) {
      xy_lims <- get_xy_lims(topGenes, min_x_abs=fc_cut, min_y2=-log10(p_cut))
      if (identical(x_lim, "auto")) x_lim <- xy_lims[["x"]]
      if (identical(y_lim, "auto")) y_lim <- xy_lims[["y"]]
    }

    # add "genes" column to topGenes
    topGenes$genes <- rownames(topGenes)

    # if threshold not specified, calculate it
    if (color_by_threshold & is.null(topGenes$threshold))
      topGenes$threshold <- (abs(topGenes$logFC) > fc_cut) & (topGenes$adj.P.Val < p_cut)

    # generate volcano plot
    if (color_by_threshold) {
      volcano <-
        ggplot(data = topGenes,
               aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
        geom_point(alpha=0.6, size=3, shape=16) +
        scale_colour_manual(values=my_cols)
    } else {
      volcano <-
        ggplot(data = topGenes,
               aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(alpha=0.6, size=2, shape=16, color=my_cols[1])
    }
    volcano <- volcano +
      theme(legend.position = "none") +
      xlab("log2 fold change") + ylab("-log10 Adj P")

    if (!is.null(fc_cut)) {
      volcano <- volcano +
        geom_vline(xintercept = fc_cut, linetype="dotted", size=0.5) +
        geom_vline(xintercept = -fc_cut, linetype="dotted", size=0.5)
    }

    if (!is.null(p_cut)) {
      volcano <- volcano +
        geom_hline(yintercept = -log10(p_cut), linetype="dotted",size=0.5)
    }

    if (!is.null(x_lim)) {volcano <- volcano + xlim(x_lim)}
    if (!is.null(y_lim)) {volcano <- volcano + ylim(y_lim)}

    if (!is.null(gene_labs)) {
      gene_labs <- match.arg(gene_labs, choices=c("ellipse", "threshold"))
      if (gene_labs=="ellipse") {
        topGenes.tmp <-
          topGenes[((topGenes$logFC^2)/(x_cut^2) + (log10(topGenes$adj.P.Val)^2)/(y_cut^2)) > 1,]
      } else if (gene_labs=="threshold") {
        x_cut_direction <- match.arg(x_cut_direction, choices=c("both", "lower", "upper"))
        if (x_cut_direction=="both") {
          topGenes.tmp <-
            topGenes[
              (abs(topGenes$logFC) > x_cut) & (-log10(topGenes$adj.P.Val) > y_cut),]
        } else if (x_cut_direction=="lower") {
          topGenes.tmp <-
            topGenes[
              (topGenes$logFC < x_cut) & (-log10(topGenes$adj.P.Val) > y_cut),]
        } else if (x_cut_direction=="upper") {
          topGenes.tmp <-
            topGenes[
              (topGenes$logFC > x_cut) & (-log10(topGenes$adj.P.Val) > y_cut),]
        }
      }

      if (gene_labs_repel) {
        volcano <- volcano +
          geom_text_repel(
            data=topGenes.tmp,
            aes(label=genes),
            color="black", size=gene_lab_size)
      } else {
        volcano <- volcano +
          geom_text(
            data=topGenes.tmp,
            aes(label=genes),
            color="black", size=gene_lab_size, vjust=1, hjust=0.5)
      }
    }

    # output volcano plot to file or plot window
    if (!is.null(file_prefix)) {
      pdf(file=paste(file_prefix, "pdf", sep="."), w=plotdims[1], h=plotdims[2], ...)
      on.exit(dev.off(), add=TRUE) # close plotting device on exit
    } else quartz(plotdims[1],plotdims[2])
    print(volcano)
  }

for (i in names(topGenes.group.stimulus.with_interaction.interaction_by_stimulus)) {
  plot_volcano_2var(
    topGenes.group.stimulus.with_interaction.interaction_by_stimulus[[i]],
    file_prefix=paste0("plots/TruCulture/volcano.group.stimulus.with_interaction.interaction_combat_quant_batch_corrected_filtered2_AdjPvalCutoff0.05_", i),
    plotdims=c(6,6),
    x_lim=xy_lims.tmp$x, y_lim=xy_lims.tmp$y,
    fc_cut=fc_cut.tmp, p_cut=p_cut.tmp,
    color_by_threshold=FALSE, my_cols="black",
    gene_labs="ellipse", x_cut=2.5, y_cut=1.4, gene_lab_size=5,gene_labs_repel= TRUE)
  # gene_labs="threshold", x_cut=0, y_cut=1.4, gene_lab_size=6)
}



#alternative volcano plots using code from Daniel Nachun
# volcanoPlot <- function(volcano.df, filename, plot.title, pp.column, pp.threshold, log.column = "logFC", xlabel = "Log Fold Change", rotate.x = FALSE) {
#   sig.df <- filter_(volcano.df, str_c("adj p val <= 0.05 & ", pp.column, " > ", pp.threshold))
#   notsig.df <- filter_(volcano.df, str_c("adj p val > 0.05 | ", pp.column, " < ", pp.threshold))
#   plot.title <- str_c(plot.title, " (", nrow(sig.df), "/", nrow(volcano.df), " Transcripts)")
#   
#   p <- ggplot() +
#     geom_point(aes_string(x = log.column, y = "logBF"), color = "gray", notsig.df) + 
#     geom_point(aes_string(x = log.column, y = "logBF"), color = "mediumblue", sig.df) + 
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.border = element_rect(size = 1, color = "black"),
#           legend.position = "none",
#           plot.background = element_blank(),
#           plot.title = element_text(hjust = 0.5)) + 
#     xlab(xlabel) + ylab("pvalue") + 
#     xlim(c(-0.4, 0.4)) +
#     ggtitle(plot.title) 
#   
#   CairoPDF(filename, width = 5, height = 5, bg = "transparent")
#   print(p)
#   dev.off()
# }







