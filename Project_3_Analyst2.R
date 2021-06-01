#analyst project 3
#part 6
#computing concordance 

#read limma files 
limma_flu <- read.csv("/projectnb/bf528/users/frazzled/project_3/analyst/limma_results_flu.csv")
limma_ifo <- read.csv("/projectnb/bf528/users/frazzled/project_3/analyst/limma_results_ifo.csv")
limma_lef <- read.csv("/projectnb/bf528/users/frazzled/project_3/analyst/limma_results_lef.csv")

#subset limma results by pval and logFC
limma_flu <- subset(limma_flu, P.Value < 0.05)
limma_ifo <- subset(limma_ifo, P.Value < 0.05)
limma_lef <- subset(limma_lef, P.Value < 0.05)

limma_flu <- subset(limma_flu, abs(logFC) >1.5)
limma_ifo <- subset(limma_ifo, abs(logFC) >1.5)
limma_lef <- subset(limma_lef, abs(logFC) >1.5)

#open deseq results
deseq_flu <- read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/CAR_PXR_DESeq_results_padj.csv")
deseq_ifo <- read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/DNA_DAMAGE_DESeq_results_padj.csv")
deseq_lef <- read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/AhR_DESeq_results_padj.csv")

#subset deseq results by pval and logFC
deseq_flu <- subset(deseq_flu, pvalue < 0.05)
deseq_ifo <- subset(deseq_ifo, pvalue < 0.05)
deseq_lef <- subset(deseq_lef, pvalue < 0.05)

deseq_flu <- subset(deseq_flu, abs(log2FoldChange) >1.5)
deseq_ifo <- subset(deseq_ifo, abs(log2FoldChange) >1.5)
deseq_lef <- subset(deseq_lef, abs(log2FoldChange) >1.5)

#read affy map 
affy_map <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")


#computing concordance 

#match ids 
limma_flu_match <- affy_map$PROBEID %in% limma_flu$X
deseq_flu_match <- affy_map$REFSEQ %in% deseq_flu$X

limma_index_flu <- limma_flu$logFC[match(affy_map$PROBEID,limma_flu$X)]
deseq_index_flu <- deseq_flu$log2FoldChange[match(affy_map$REFSEQ,deseq_flu$X)]
samedir_logFC_flu <- sign(limma_index_flu) == sign(deseq_index_flu)

# ----- FLUCONAZOLE 
intersect <- limma_flu_match & deseq_flu_match & samedir_logFC_flu
n_0 <- sum(intersect)
total <- sum(limma_flu_match | deseq_flu_match)
n_1 <- dim(deseq_flu)[1]
n_2 <- dim(limma_flu)[1]
#N <- n_1 + n_2
#which N to use? varies greatly.... 
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_flu <- 2 * n_x / total
# 0.12135449... 



# ----- IFOSFAMIDE
limma_ifo_match <- affy_map$PROBEID %in% limma_ifo$X
deseq_ifo_match <- affy_map$REFSEQ %in% deseq_ifo$X

limma_index_ifo <- limma_ifo$logFC[match(affy_map$PROBEID,limma_ifo$X)]
deseq_index_ifo <- deseq_ifo$log2FoldChange[match(affy_map$REFSEQ,deseq_ifo$X)]
samedir_logFC_ifo <- sign(limma_index_ifo) == sign(deseq_index_ifo)

intersect <- limma_ifo_match & deseq_ifo_match & samedir_logFC_ifo
n_0 <- sum(intersect)
total <- sum(limma_ifo_match | deseq_ifo_match)
n_1 <- dim(deseq_ifo)[1]
n_2 <- dim(limma_ifo)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_ifo <- 2 * n_x / total
# 5.08 *E-05


# ------ LEFLUNOMIDE 
limma_lef_match <- affy_map$PROBEID %in% limma_lef$X
deseq_lef_match <- affy_map$REFSEQ %in% deseq_lef$X

limma_index_lef <- limma_lef$logFC[match(affy_map$PROBEID,limma_lef$X)]
deseq_index_lef <- deseq_lef$log2FoldChange[match(affy_map$REFSEQ,deseq_lef$X)]
samedir_logFC_lef <- sign(limma_index_lef) == sign(deseq_index_lef)

intersect <- limma_lef_match & deseq_lef_match & samedir_logFC_lef
n_0 <- sum(intersect)
total <- sum(limma_lef_match | deseq_lef_match)
n_1 <- dim(deseq_lef)[1]
n_2 <- dim(limma_lef)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_lef <- 2 * n_x / total
# 0.0448827... 



library(ggplot2)
library(ggpubr)
#remaking Fig 2A 

#plot concordance vs number of DE genes for microarray analysis
scatter_limma <- data.frame("Concordance"= c(C_flu, C_ifo, C_lef), "Treatment"= c(nrow(limma_flu), nrow(limma_ifo), nrow(limma_lef)))
labels <- c("CAR", "IFO", "LEF")
plot1 <- ggplot(scatter_limma, aes(x=Treatment, y=Concordance)) +
  geom_point() +
  labs(title="DEGs from Microarray Analysis", x="Treatment Effect", y="Concordance of DEGs") +
  geom_text(label=labels, nudge_x = c(6, 6, 6), nudge_y = c(0, 0, 0)) 

#plot concordance vs number of DE genes for RNA-seq analysis  
scatter_deseq <- data.frame("Concordance"= c(C_flu, C_ifo, C_lef), "Treatment" = c(nrow(deseq_flu), nrow(deseq_ifo), nrow(deseq_lef)))
plot2 <- ggplot(scatter_deseq, aes(x=Treatment, y=Concordance)) +
  geom_point() +
  labs(title="DEGs from RNA-Seq", x="Treatment Effect", y="Concordance of DEGs") +
  geom_text(label=labels, nudge_x = c(10, 10, 10), nudge_y = c(0.005, 0.005, 0.005))

#arrange plots into one figure
ggarrange(plot1, plot2,
          labels = c('A', 'B'),
          ncol = 2, nrow = 1)






# Fig 3 B & C?

#separating above / below median 
med_limma_flu <- median(limma_flu$AveExpr)
above_limma_flu <- subset(limma_flu, AveExpr > med_limma_flu)
below_limma_flu <- subset(limma_flu, AveExpr < med_limma_flu)

med_deseq_flu <- median(deseq_flu$baseMean)
above_deseq_flu <- subset(deseq_flu, baseMean > med_deseq_flu) 
below_deseq_flu <- subset(deseq_flu, baseMean < med_deseq_flu)





med_limma_ifo <- median(limma_ifo$AveExpr)
above_limma_ifo <- subset(limma_ifo, AveExpr > med_limma_ifo)
below_limma_ifo <- subset(limma_ifo, AveExpr < med_limma_ifo)

med_deseq_ifo <- median(deseq_ifo$baseMean)
above_deseq_ifo <- subset(deseq_ifo, baseMean > med_deseq_ifo) 
below_deseq_ifo <- subset(deseq_ifo, baseMean < med_deseq_ifo)




med_limma_lef <- median(limma_lef$AveExpr)
above_limma_lef <- subset(limma_lef, AveExpr > med_limma_lef)
below_limma_lef <- subset(limma_lef, AveExpr < med_limma_lef)

med_deseq_lef <- median(deseq_lef$baseMean)
above_deseq_lef <- subset(deseq_lef, baseMean > med_deseq_lef) 
below_deseq_lef <- subset(deseq_lef, baseMean < med_deseq_lef)


#compute new concordance values

# ----- FLUCONAZOLE ----- ABOVE MEDIAN 
above_limma_flu_match <- affy_map$PROBEID %in% above_limma_flu$X
above_deseq_flu_match <- affy_map$REFSEQ %in% above_deseq_flu$X

above_limma_index_flu <- above_limma_flu$logFC[match(affy_map$PROBEID,above_limma_flu$X)]
above_deseq_index_flu <- above_deseq_flu$log2FoldChange[match(affy_map$REFSEQ,above_deseq_flu$X)]
above_samedir_logFC_flu <- sign(above_limma_index_flu) == sign(above_deseq_index_flu)

intersect <- above_limma_flu_match & above_deseq_flu_match & above_samedir_logFC_flu
n_0 <- sum(intersect)
total <- sum(above_limma_flu_match | above_deseq_flu_match)
n_1 <- dim(above_deseq_flu)[1]
n_2 <- dim(above_limma_flu)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
ABOVE_C_flu <- 2 * n_x / total
#0.12604

# ----- IFO ----- ABOVE MEDIAN 
above_limma_ifo_match <- affy_map$PROBEID %in% above_limma_ifo$X
above_deseq_ifo_match <- affy_map$REFSEQ %in% above_deseq_ifo$X

above_limma_index_ifo <- above_limma_ifo$logFC[match(affy_map$PROBEID,above_limma_ifo$X)]
above_deseq_index_ifo <- above_deseq_ifo$log2FoldChange[match(affy_map$REFSEQ,above_deseq_ifo$X)]
above_samedir_logFC_ifo <- sign(above_limma_index_ifo) == sign(above_deseq_index_ifo)

intersect <- above_limma_ifo_match & above_deseq_ifo_match & above_samedir_logFC_ifo
n_0 <- sum(intersect)
total <- sum(above_limma_ifo_match | above_deseq_ifo_match)
n_1 <- dim(above_deseq_ifo)[1]
n_2 <- dim(above_limma_ifo)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
ABOVE_C_ifo <- 2 * n_x / total
#0


# ----- LEF ----- ABOVE MEDIAN 
above_limma_lef_match <- affy_map$PROBEID %in% above_limma_lef$X
above_deseq_lef_match <- affy_map$REFSEQ %in% above_deseq_lef$X

above_limma_index_lef <- above_limma_lef$logFC[match(affy_map$PROBEID,above_limma_lef$X)]
above_deseq_index_lef <- above_deseq_lef$log2FoldChange[match(affy_map$REFSEQ,above_deseq_lef$X)]
above_samedir_logFC_lef <- sign(above_limma_index_lef) == sign(above_deseq_index_lef)

intersect <- above_limma_lef_match & above_deseq_lef_match & above_samedir_logFC_lef
n_0 <- sum(intersect)
total <- sum(above_limma_lef_match | above_deseq_lef_match)
n_1 <- dim(above_deseq_lef)[1]
n_2 <- dim(above_limma_lef)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
ABOVE_C_lef <- 2 * n_x / total
#0.038794...



# ----- FLUCONAZOLE ----- BELOW MEDIAN 
below_limma_flu_match <- affy_map$PROBEID %in% below_limma_flu$X
below_deseq_flu_match <- affy_map$REFSEQ %in% below_deseq_flu$X

below_limma_index_flu <- below_limma_flu$logFC[match(affy_map$PROBEID,below_limma_flu$X)]
below_deseq_index_flu <- below_deseq_flu$log2FoldChange[match(affy_map$REFSEQ,below_deseq_flu$X)]
below_samedir_logFC_flu <- sign(below_limma_index_flu) == sign(below_deseq_index_flu)

intersect <- below_limma_flu_match & below_deseq_flu_match & below_samedir_logFC_flu
n_0 <- sum(intersect)
total <- sum(below_limma_flu_match | below_deseq_flu_match)
n_1 <- dim(below_deseq_flu)[1]
n_2 <- dim(below_limma_flu)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
BELOW_C_flu <- 2 * n_x / total
#0.009909....


# ----- IFO ----- BELOW MEDIAN 
below_limma_ifo_match <- affy_map$PROBEID %in% below_limma_ifo$X
below_deseq_ifo_match <- affy_map$REFSEQ %in% below_deseq_ifo$X

below_limma_index_ifo <- below_limma_ifo$logFC[match(affy_map$PROBEID,below_limma_ifo$X)]
below_deseq_index_ifo <- below_deseq_ifo$log2FoldChange[match(affy_map$REFSEQ,below_deseq_ifo$X)]
below_samedir_logFC_ifo <- sign(below_limma_index_ifo) == sign(below_deseq_index_ifo)

intersect <- below_limma_ifo_match & below_deseq_ifo_match & below_samedir_logFC_ifo
n_0 <- sum(intersect)
total <- sum(below_limma_ifo_match | below_deseq_ifo_match)
n_1 <- dim(below_deseq_ifo)[1]
n_2 <- dim(below_limma_ifo)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
BELOW_C_ifo <- 2 * n_x / total
#0



# ----- LEF ----- BELOW MEDIAN 
below_limma_lef_match <- affy_map$PROBEID %in% below_limma_lef$X
below_deseq_lef_match <- affy_map$REFSEQ %in% below_deseq_lef$X

below_limma_index_lef <- below_limma_lef$logFC[match(affy_map$PROBEID,below_limma_lef$X)]
below_deseq_index_lef <- below_deseq_lef$log2FoldChange[match(affy_map$REFSEQ,below_deseq_lef$X)]
below_samedir_logFC_lef <- sign(below_limma_index_lef) == sign(below_deseq_index_lef)

intersect <- below_limma_lef_match & below_deseq_lef_match & below_samedir_logFC_lef
n_0 <- sum(intersect)
total <- sum(below_limma_lef_match | below_deseq_lef_match)
n_1 <- dim(below_deseq_lef)[1]
n_2 <- dim(below_limma_lef)[1]
N <- dim(affy_map)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
BELOW_C_lef <- 2 * n_x / total
#0.00811029...





#bar plot combining overall concordance measures 
bar_plot <- data.frame("Concordance"=c(ABOVE_C_flu, BELOW_C_flu, C_flu, ABOVE_C_ifo, BELOW_C_ifo, C_ifo, ABOVE_C_lef, BELOW_C_lef, C_lef),
                  "Analysis"=c("Above", "Below", "Overall", "Above", "Below", "Overall", "Above", "Below", "Overall"),
                  "Chemical"=c("Fluconazole","Fluconazole","Fluconazole",
                               "Ifosfamide", "Ifosfamide", "Ifosfamide",
                               "Leflunomide", "Leflunomide", "Leflunomide"))

# Plot data with a side-by-side bar chart
ggplot(bar_plot, aes(fill=Analysis, x=Chemical, y=Concordance)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(title="Above and Below Median Concordance")