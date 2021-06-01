#individual project

#project 2
#analyst 

data <- read.table('/projectnb/bf528/users/frazzled/project_2/project-2-frazzled/cuffdiff_out/gene_exp.diff', header=T)

data <- data[order(data$q_value),]

diff_data <- data[1:10,]
top_10 <- diff_data[ -c(1:2, 4:7, 11,14) ]
write.table(top_10,'/projectnb/bf528/users/frazzled/project_5/mdc98/top_10_genes.txt',sep = '\t')


hist <- hist(data$log2.fold_change., breaks = 20, labels = TRUE,
             main = "Histogram of log2 fold change",
             xlab = "log2 fold change values")  

sig_data <- data[which(data$significant=="yes"),]

hist2 <- hist(sig_data$log2.fold_change., breaks=20, labels = TRUE,
     main = "Histogram of significant log2 fold change values",
     xlab = "log2 fold change values")




up_reg <- subset(sig_data,log2.fold_change.>0)
down_reg <- subset(sig_data,log2.fold_change.<0)


write.table(up_reg[3], "/projectnb/bf528/users/frazzled/project_5/mdc98/up_reg.txt", row.names = F, col.name = F, quote = F)
write.table(down_reg[3], "/projectnb/bf528/users/frazzled/project_5/mdc98/down_reg.txt", row.names = F, col.name = F, quote = F)