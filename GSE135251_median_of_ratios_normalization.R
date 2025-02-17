#Have a look at the dataset GSE135251 (MASHLD), to have a look at differentially expressed genes
#DOI: 10.1126/scitranslmed.aba4448 
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(ggpubr)
library(rstatix)
library(org.Hs.eg.db)

#import the raw count data???
GSE135251 <- getRNASeqData("GSE135251")

#assemble DESeqDataSet
dds <- DESeqDataSet(GSE135251, design = ~ group.in.paper.ch1)
#get more recognizable gene names
rownames(dds) <- rowData(dds)$Symbol

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)



#which genes to compare
genes <- c("RBCK1", "SHARPIN", "RNF31", "CASP3", "BAX", "BAK1", "RIPK1")

#get the counts, normalized (through median of ratios)
counts <-  normalized_counts[genes, dds$geo_accession[dds$group.in.paper.ch1!="control"]]

#order the rows in counts according to the gene vector provided
counts <- counts[genes, ]



#turn into dataframe with only one column for the data and the other columns as descriptors
counts <- data.frame(gene = rep(row.names(counts), times =length(colnames(counts))), 
                     label = rep(colData(dds)[colnames(counts), "group.in.paper.ch1"], each = length(genes)),
                     count = unlist(as.list(counts)))
counts$label <- factor(counts$label, levels = c("NAFL", "NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4"))
counts$gene <- factor(counts$gene, levels = genes)


#with geom_pwc, ggadjust_pvalue adjusts the values across panels
#Is this too conservative????
ggadjust_pvalue(ggplot(counts, aes(x = label, y = log10(count), fill = label)) +
                  geom_boxplot(outliers = FALSE,
                               staplewidth = 0.5) +
                  scale_fill_manual(values = c("#D3D3D3", "#FECCCB", "#FF9899", "#FC6667", "#FE0000")) +
                  facet_wrap(~ gene, ncol = 2, scales = "free") +
                  geom_point(position = position_jitter(height = 0), color = "black", size = 1) +
                  geom_pwc(ref.group = "NAFL",
                           label = "p", 
                           method = "wilcox_test", 
                           method.args = list(alternative = "two.sided")
                  ) +
                  ylab("normalized counts") +
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        line = element_line(color = "black", linewidth = 0.5),
                        legend.position = "none",
                        strip.background = element_blank(),
                        axis.line.x.bottom = element_line(),
                        axis.line.y.left = element_line(),
                        axis.line.x.top = element_blank(),
                        axis.line.y.right = element_blank(),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_blank()
                  ) +
                  ylab("normalized counts") +
                  xlab(NULL) +
                  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
                  ggtitle("conservative testing"),
                p.adjust.method = "BH",
                label = "p.adj",
                hide.ns = TRUE
)

#ggsave("figure_output/GSE135251_cell_death_comparison.svg", device = "svg", scale = 3)







#alternative, less conservative testing under the assumption, that the different genes are independent of each other
'
#let padj get calculate in each panel on its own
ggplot(counts, aes(x = label, y = count, fill = label)) +
  geom_boxplot(outliers = FALSE,
               staplewidth = 0.5) +
  scale_fill_manual(values = c("#D3D3D3", "#FECCCB", "#FF9899", "#FC6667", "#FE0000")) +
  facet_wrap(~ gene, ncol = 2, scales = "free") +
  geom_point(position = position_jitter(height = 0), 
             color = "black", size = 1) +
  geom_pwc(ref.group = "NAFL",
           label = "p.adj", 
           method = "wilcox_test",
           p.adjust.method = "bonferroni",
           hide.ns = FALSE,
           remove.bracket = FALSE,
           method.args = list(alternative = "two.sided")
           ) +
  ylab("normalized counts") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(color = "black", linewidth = 0.5),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()
  ) +
  ylab("normalized counts") +
  xlab(NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  ggtitle("unconservative testing")



#make custom statistics first, then add them to the plot (this does the same as the above)
stat.test <- counts |> 
  group_by(gene) |> 
  wilcox_test(formula = count ~ label,
              ref.group = "NAFL",
              p.adjust.method = "bonferroni"
              )
stat.test <- stat.test |> 
  add_y_position()

stat.test$y.position <- stat.test$y.position*0.8

ggplot(counts, aes(x = label, y = count)) +
  geom_boxplot(aes(fill = label),
               outliers = FALSE,
               staplewidth = 0.5) +
  scale_fill_manual(values = c("#D3D3D3", "#FECCCB", "#FF9899", "#FC6667", "#FE0000")) +
  facet_wrap(~ gene, ncol = 2, scales = "free") +
  geom_point(position = position_jitter(height = 0),
             color = "black",
             size = 1) +
  ylab("normalized counts") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(color = "black", linewidth = 0.5),
        legend.position = "none",
        strip.background = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()
        ) +
  ylab("normalized counts") +
  xlab(NULL) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     label = "p.adj"
                     ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
'
