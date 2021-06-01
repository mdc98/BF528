#calling to use tidyverse packages 
library(tidyverse)
library(dplyr)

#reading fpkm_tracking tables 
Ad_1 <- read.table("/projectnb/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking", header = TRUE)
Ad_2 <- read.table("/projectnb/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking", header = TRUE)
P0_1 <- read.table("/projectnb/bf528/users/frazzled/project_2/ledia/P0_1_cufflinks/cuffdiff_out/genes.fpkm_tracking", header = TRUE) %>% 
  rename(FPKM = P0_FPKM)
P0_2 <- read.table("/projectnb/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking", header = TRUE)
P4_1 <- read.table("/projectnb/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking", header = TRUE)
P4_2 <- read.table("/projectnb/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking", header = TRUE)
P7_1 <- read.table("/projectnb/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking", header = TRUE)
P7_2 <- read.table("/projectnb/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking", header = TRUE)

#assigning gene names from Fig. 1D 
genes_sar <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
genes_mit <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
genes_cc <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")

#function to get fpkm value
subset.gene <- function(x1,x2,x3,x4,target.gene){
 x1 <- x1 %>% 
   filter(gene_short_name %in% target.gene) %>% 
   select(gene_short_name, FPKM) %>% 
   mutate(sample = 1)
 x2 <- x2 %>% 
   filter(gene_short_name %in% target.gene) %>% 
   select(gene_short_name, FPKM) %>% 
   mutate(sample = 2)
 x3 <- x3 %>% 
   filter(gene_short_name %in% target.gene) %>% 
   select(gene_short_name, FPKM) %>% 
   mutate(sample = 3)
 x4 <- x4 %>% 
   filter(gene_short_name %in% target.gene) %>% 
   select(gene_short_name, FPKM) %>% 
   mutate(sample = 4)
 return(rbind(x1, x2, x3, x4))
}

#running function 
sar_gene1.fpkm <- subset.gene(P0_1,P4_1,P7_1,Ad_1,genes_sar)
sar_gene1.fpkm
mit_gene1.fpkm <- subset.gene(P0_1,P4_1,P7_1,Ad_1,genes_mit)
mit_gene1.fpkm
cc_gene1.fpkm <- subset.gene(P0_1,P4_1,P7_1,Ad_1,genes_cc)
cc_gene1.fpkm
sar_gene2.fpkm <- subset.gene(P0_2,P4_2,P7_2,Ad_2,genes_sar)
sar_gene2.fpkm
mit_gene2.fpkm <- subset.gene(P0_2,P4_2,P7_2,Ad_2,genes_mit)
mit_gene2.fpkm
cc_gene2.fpkm <- subset.gene(P0_2,P4_2,P7_2,Ad_2,genes_cc)
cc_gene2.fpkm


#plotting function 

p <- ggplot(data = sar_gene1.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample 2", y = "FPKM", title = "Sarcomere (Sample 1)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")

p <- ggplot(data = sar_gene2.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample 2", y = "FPKM", title = "Sarcomere (Sample 2)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")

p <- ggplot(data = mit_gene1.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line()+
  labs(x = "Sample 1", y = "FPKM", title = "Mitochondria (Sample 1)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")

p <- ggplot(data = mit_gene2.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line()+
  labs(x = "Sample 2", y = "FPKM", title = "Mitochondria (Sample 2)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")

p <- ggplot(data = cc_gene1.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line()+
  labs(x = "Sample 1", y = "FPKM", title = "Cell Cycle (Sample 1)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")

p <- ggplot(data = cc_gene2.fpkm, aes(x = sample, y = FPKM, color = gene_short_name)) +
  geom_point() +
  geom_line()+
  labs(x = "Sample 2", y = "FPKM", title = "Cell Cycle (Sample 2)", color = "Gene Name") +
  theme(plot.title = element_text(hjust = 0.5))
p
p + scale_x_discrete("P0                  P4                  P7                  Ad")






#7.2 

#open DAVID files, and paper files to compare 
downreg <- read.csv("/projectnb/bf528/users/frazzled/project_2/downregulated_genes_DAVID.csv",header=T)
downreg_ref <- read.csv("/projectnb/bf528/users/frazzled/project_2/paper_down_reg.csv",header = T)
#match up terms in common with an asterick 
match_downreg <- match(downreg[,2], downreg_ref[,2])
match_downreg[!is.na(match_downreg)] <- "*"
downreg <- cbind(downreg, match_downreg)
#making down regulated table 
down_table <- downreg[,c(2,11,12,13,14)]
down_table <- down_table[1:50,]
write.csv(down_table,"DAVID_GOterms_down.csv")


#open DAVID file / paper file 
upreg <- read.csv("/projectnb/bf528/users/frazzled/project_2/upregulated_genes_DAVID.csv",header=T)
upreg_ref <- read.csv("/projectnb/bf528/users/frazzled/project_2/paper_up_reg.csv",header = T)
match_upreg <- match(upreg[,2], upreg_ref[,2])
match_upreg[!is.na(match_upreg)] <- "*"
upreg <- cbind(upreg, match_upreg)
#making up regulated table
up_table <- upreg[,c(2,11,12,13,14)]
up_table <- up_table[1:50,]
write.csv(up_table,"DAVID_GOterms_up.csv")







#Part 3, heatmap:  
#get FPKM values avg for each sample/gene name

P0_1 <- P0_1[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
P0_2 <- P0_2[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
P4_1 <- P4_1[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
P4_2 <- P4_2[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
P7_1 <- P7_1[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
P7_2 <- P7_2[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
Ad_1 <- Ad_1[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)
Ad_2 <- Ad_2[c("gene_short_name","FPKM")] %>% 
  group_by(gene_short_name) %>% 
  summarise_all(sum)

#gene expression diff
diff_genes <- read.table("/projectnb/bf528/users/frazzled/project_2/project-2-frazzled/cuffdiff_out/gene_exp.diff", header = TRUE)
#get top 500 diff expressed genes by qval
diff_genes <- diff_genes %>% 
  arrange(q_value) %>% 
  head(500)

df_all <- P0_1 %>% 
  left_join(P0_2,"gene_short_name") %>% 
  left_join(P4_1,"gene_short_name") %>% 
  left_join(P4_2,"gene_short_name") %>% 
  left_join(P7_1,"gene_short_name") %>% 
  left_join(P7_2,"gene_short_name") %>% 
  left_join(Ad_1,"gene_short_name") %>% 
  left_join(Ad_2,"gene_short_name")
df_all <- df_all %>% 
  filter(gene_short_name %in% diff_genes$gene) %>% 
  filter(complete.cases(.))
df_all <- df_all %>%
  rename(P0_1 = FPKM.x) %>%
  rename(P0_2 = FPKM.y) %>%
  rename(P4_1 = FPKM.x.x) %>%
  rename(P4_2 = FPKM.y.y) %>%
  rename(P7_1 = FPKM.x.x.x) %>%
  rename(P7_2 = FPKM.y.y.y) %>%
  rename(Ad_1 = FPKM.x.x.x.x) %>%
  rename(Ad_2 = FPKM.y.y.y.y)

#create heatmap with genes along rows, and samples along columns
  #add dendrograms and labels 
heatmap(as.matrix(df_all[2:9]))