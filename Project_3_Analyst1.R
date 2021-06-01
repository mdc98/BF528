#Analyst Project 3
#part 5
#run limma & plotting histograms, scatterplots 

library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_3_mic_info.csv', header = TRUE, as.is = TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

# subset the full expression matrix to just those in this comparison
rma.subset.lef <- rma[paste0('X',samples$array_id[samples$chemical=='LEFLUNOMIDE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]
rma.subset.flu <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]
rma.subset.ifo <- rma[paste0('X',samples$array_id[samples$chemical=='IFOSFAMIDE' | (samples$chemical=='Control' & samples$vehicle=='SALINE_100_%')])]

# construct a design matrix modeling treatment vs control for use by limma
design.lef <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','LEFLUNOMIDE')
  )
)

colnames(design.lef) <- c('Intercept','LEFLUNOMIDE')

design.flu <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)

colnames(design.flu) <- c('Intercept','FLUCONAZOLE')

design.ifo <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='SALINE_100_%'],
    levels=c('Control','IFOSFAMIDE')
  )
)

colnames(design.ifo) <- c('Intercept','IFOSFAMIDE')


# run limma 
fit_lef <- lmFit(rma.subset.lef, design=design.lef)
fit_lef <- eBayes(fit_lef)
t_lef <- topTable(fit_lef, coef=2, n=nrow(rma.subset.lef), adjust='BH')

dim(rma.subset.flu)
dim(design.flu)


# write out the results to file
write.csv(t_lef,'limma_results_lef.csv')

# run limma 
fit_flu <- lmFit(rma.subset.flu, design=design.flu)
fit_flu <- eBayes(fit_flu)
t_flu <- topTable(fit_flu, coef=2, n=nrow(rma.subset.flu), adjust='BH')
#write results into file 
write.csv(t_flu, 'limma_results_flu.csv')

# run limma 
fit_ifo <- lmFit(rma.subset.ifo, design=design.ifo)
fit_ifo <- eBayes(fit_ifo)
t_ifo <- topTable(fit_ifo, coef=2, n=nrow(rma.subset.ifo), adjust='BH')
#write results into file 
write.csv(t_ifo, 'limma_results_ifo.csv')



# get significant results

#based off adjusted p value 
sig.lef <- subset(t_lef, adj.P.Val < 0.05)
#note: 466 sig genes

sig.flu <- subset(t_flu, adj.P.Val < 0.05)
#note: 1997 sig genes

sig.ifo <- subset(t_ifo, adj.P.Val < 0.05)
#note: 0 significant genes 

#sig genes based off nominal p value 
lef_significant <- t_lef[t_lef$P.Value < 0.05,]
flu_significant <- t_flu[t_flu$P.Value < 0.05,]
ifo_significant <- t_ifo[t_ifo$P.Value < 0.05,]

#note:
  #need table of top 10 DE genes

library(ggplot2)
library('ggpubr')
#plotting histograms of fold changes
p1 <- ggplot(lef_significant, aes(x=logFC)) +
  geom_histogram(binwidth = 0.1) +
  labs(title = 'Leflunomide', x = 'Log Fold Change')

p2 <- ggplot(flu_significant, aes(x=logFC)) +
  geom_histogram(binwidth = 0.1) +
  labs(title = 'Fluconazole', x = 'Log Fold Change')

p3 <- ggplot(ifo_significant, aes(x=logFC)) +
  geom_histogram(binwidth = 0.1) +
  labs(title = 'Ifosfamide', x = 'Log Fold Change')

#arrange plots into one figure
ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

#plotting scatterplots 
# fold change vs nominal p value 
s1 <- ggplot(lef_significant, aes(x=logFC, y=P.Value)) +
  geom_point(size = 0.25) +
  labs(title = 'Leflunomide Plot', x = 'Log Fold Change')

s2 <- ggplot(flu_significant, aes(x=logFC, y=P.Value)) +
  geom_point(size = 0.25) +
  labs(title = 'Fluconazole Plot')

s3 <- ggplot(ifo_significant, aes(x=logFC, y=P.Value)) +
  geom_point(size = 0.25) +
  labs(title = 'Ifosfamide Plot')

#arrange plots into one figure 
ggarrange(s1, s2, s3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)