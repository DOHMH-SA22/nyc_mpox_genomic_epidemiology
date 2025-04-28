##### NYC Mpox genomic epidemiology

library(tidyverse)
library(rstatix)
library(ggpubr)
library(incidence2)
library(data.table)
library(lubridate)


############ Figure 1A ################

lineage_ct <- read_tsv("lineage_counts.txt", show_col_types = FALSE)

lineage_ct %>% ggplot(aes(Lineage)) + geom_bar(aes(y=counts), stat="identity")  +
  theme_bw() + labs(x = "Lineage", y = "Sequence counts") + facet_wrap(~Geography, ncol = 1) + 
  theme(axis.text.x = element_text(size=14, angle = 300, hjust = 0.05), axis.text.y = element_text(size=14), axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), strip.text.x = element_text(size=14))

############ Figure 3 ################

apo_muts_ct <- read_tsv("apobec_muts_counts.txt", show_col_types = FALSE)

apo_muts_ct$group <- factor(apo_muts_ct$group, levels = c("Pre-Outbreak Global", "Outbreak Global", "Outbreak NYC"))

stat_test <- apo_muts_ct %>% 
  group_by(type) %>% 
  t_test(mut_ct ~ group) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance()

stat_test2 <-stat_test %>% filter(type=="APOBEC3")  %>%
  mutate(y.position = c(80, 90, 100))

apo_muts_ct %>%
  ggplot(aes(x=group, y=mut_ct, color=type)) + geom_boxplot(outlier.colour = "#525252") +scale_color_manual(values = c("#a63603", "#0570b0")) +
  theme_bw() + 
  theme(axis.text = element_text(size=10), legend.text = element_text(size=10),
        legend.title = element_blank(), legend.position = "top", axis.title= element_text(size=10), legend.background = element_rect(color="#969696")) + labs(x="Group", y="Number of Mutations per Sequence") + stat_pvalue_manual(stat_test2, label = "p.adj.signif", size = 7, )


############ Figure 5A and 5B ################

intra_host <- read_tsv("intra_host.txt", show_col_types = FALSE)

intra_host %>% filter(apobec_prop >= 0.5, apobec_prop < 1)

p1 <- intra_host %>% dplyr::mutate(category = if_else(apobec_prop == 0, "Without-APOBEC", "With-APOBEC")) %>% ggplot(aes(as.character(snpdists))) + geom_bar(aes(fill=category)) + scale_fill_manual(values = c("#a63603", "#0570b0"))  + theme_bw() + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12),
        legend.title = element_blank(), legend.position = "bottom", legend.background = element_rect(color="#969696")) + labs(x="Maximum SNP Distance Between Sequences within a Patient", y="Patient Count")

p2 <- intra_host  %>% ggplot(aes(apobec_prop)) + geom_histogram() + theme_bw() + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12)) + labs(x="Proportion of SNPs Between Sequences due to APOBEC", y="Patient Count")

ggarrange(p1, p2, nrow=2, common.legend = TRUE, labels = c("A", "B"), font.label = list(size = 14))

############ Supp Figure 1 ################

fig1a <- read_tsv("data_sup_fig_1A.txt", show_col_types = FALSE)

p1 <- fig1a %>%
  ggplot(aes(x=as.factor(collection_month),y=Count, fill=Lesion)) + geom_bar(stat="identity",  position = position_dodge()) + 
  scale_fill_manual(values = c("#f46d43", "#4393c3")) + theme_bw() + 
  theme(axis.text.x = element_text(size=12, angle=90, vjust=0.6, hjust=0.1), axis.text.y = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = "top", axis.title = element_text(size=12), legend.background = element_rect(color="#969696") ) + labs( x= "Month of Collection", y = "# of Patients")

fig1b <- read_tsv("data_sup_fig_1B.txt", show_col_types = FALSE)

p2 <- fig1b %>%
  ggplot(aes(x=Borough,y=Count, fill=Lesion)) + geom_bar(stat="identity",  position = position_dodge()) + 
  scale_fill_manual(values = c("#f46d43", "#4393c3")) + theme_bw()+ 
  theme(axis.text.x = element_text(size=12, angle=90, vjust=0.6, hjust=0.1), axis.text.y = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = "top", axis.title = element_text(size=12), legend.background = element_rect(color="#969696") ) + labs( x= "NYC Borough", y = "# of Patients")


fig1c <- read_tsv("data_sup_fig_1C.txt", show_col_types = FALSE)

p3 <- fig1c %>%
  ggplot(aes(x=Site,y=Count, fill=Lesion)) + geom_bar(stat="identity",  position = position_dodge()) + 
  scale_fill_manual(values = c("#f46d43", "#4393c3")) + theme_bw()+ 
  theme(axis.text.x = element_text(size=12, angle=90, vjust=0.6, hjust=0.1), axis.text.y = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = "top", axis.title = element_text(size=12), legend.background = element_rect(color="#969696") ) + labs( x= "Body Site", y = "Sequenced Specimens")

fig1d <- read_tsv("data_sup_fig_1D.txt", show_col_types = FALSE)

p4 <- fig1d %>%
  ggplot(aes(x=lesion_per_site,y=specimen_ct, fill=Lesion)) + geom_bar(stat="identity") + 
  scale_fill_manual(values = c("#f46d43")) + theme_bw()+ 
  theme(axis.text = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = "top", axis.title = element_text(size=12), legend.background = element_rect(color="#969696") ) + labs( x= "Lesion(s)/Site/Patient", y = "Sequenced Specimens")

ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, common.legend = TRUE, labels = c("A", "B", "C", "D"), font.label = list(size = 14))


############ Supp Figure 3 ################

muts_freq <- read_tsv("Muts_freq.txt", show_col_types = FALSE )

muts_freq %>% ggplot(aes(Pos, freq)) + geom_point(aes(color=Effect), shape=21, size=2.5) + 
  geom_segment( aes(x=Pos, xend=Pos, y=0, yend=freq), alpha=0.5) + scale_color_manual(values = c("#4d4d4d", "#e41a1c", "#377eb8")) + 
  facet_wrap(~Cluster_Id, ncol=2) + theme_bw() + theme(axis.text.x = element_text(face="bold", size=14), strip.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14), legend.text = element_text(face="bold", size=14), legend.title = element_blank(), legend.position = "top", axis.title.x = element_text(face="bold", size=14),   axis.title.y = element_text(face="bold", size=14)) + labs(x = "Position", y = "Frequency")


############ Supp Figure 8 ################

mpxv_biweekly_cts <- read_tsv("mpox_temporal_distribution.txt", show_col_types = FALSE)

mpxv_biweekly_cts %>% ggplot(aes(x=as.factor(date_index),y=count, fill=Category)) + geom_bar(stat="identity") + 
  scale_fill_manual(values = c("#f46d43", "#4393c3", "grey")) + theme_bw() + 
  theme(axis.text.x = element_text(face="bold", size=10, angle=90, vjust=0.6, hjust=0.1), axis.text.y = element_text(face="bold", size=10), legend.text = element_text(face="bold", size=10),
        legend.title = element_blank(), legend.position = c(0.87,0.9), axis.title.x = element_text(face="bold", size=10), 
        axis.title.y = element_text(face="bold", size=10), legend.background = element_rect(color="#969696") ) + 
  labs( x= "Collection dates", y = "# of Sequences")
