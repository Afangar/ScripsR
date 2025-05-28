library(tidyr)
library(dplyr)
library(ggplot2)
library(topGO)
library(clusterProfiler)
library(readxl)

setwd("Your_Folder")

#Go terms related to genes or proteins of your species, normally downloaded in Uniprot
termsgo <- read.csv("Termsgo.csv") 

#load proteins of interest
proteinasdeinteres <- read_excel("Your_File.xlsx",col_names = FALSE) #the proteins you want to perform enrichment analysis on
colnames(proteinasdeinteres)[1] <- "ID"

# Separate proteins into individual rows if multiple proteins are listed in the names, this is common in proteomics
myInterestingGenes <- proteinasdeinteres %>%
  mutate(ID = strsplit(as.character(ID), ";")) %>%
  unnest(ID)

# Make sure the ID and V2 columns of termsgo are vectors
colnames(termsgo)[1] <- "ID"
colnames(termsgo)[2] <- "V2"
termsgo$ID <- as.vector(termsgo$ID)
termsgo$V2 <- as.vector(termsgo$V2)

# Create the gene2GO list from clean termsgo
gene2GO <- strsplit(as.character(termsgo$V2), "; ")
names(gene2GO) <- termsgo$ID

# Make sure that todosgenes is a data frame with a column named ID
todosgenes <- termsgo[1]
colnames(todosgenes)[1] <- "ID"

if (!is.data.frame(todosgenes) || !("ID" %in% colnames(todosgenes))) {
  stop("todosgenes is not a valid data frame or does not have the ID column.")
}

# Convert todosgenes to a named numeric vector
geneList <- as.factor(as.integer(todosgenes$ID %in% myInterestingGenes$ID))
names(geneList) <- todosgenes$ID

# Create the topGOdata object with the cleaned data, MF: Molecular function, BP: biological process, CP: cellular
GOdata <- new("topGOdata",
              ontology = "MF"
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

#Fisher's statistical analysis,
allGO <- usedGO(GOdata)
Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# replace the values "< 1e-30" with a very small number ( < 1e-30 cannot be corrected with BH)
table$Classic[table$Classic == "< 1e-30"] <- 1e-30 
table$Classic <- as.numeric(table$Classic)  
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]

write.table(all_res_final[1:100,], file = "Your_Table.csv", quote=FALSE, row.names=FALSE, sep = ",")

# get only 20 most signigicat restults
ntop <- 20

ggdata <- all_res_final[1:ntop,]
ggdata <- ggdata[complete.cases(ggdata), ]
aux <- go2term(ggdata$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")
ggdata$Classic <- as.numeric(ggdata$Classic)

ggdata <- ggdata[order(ggdata$Classic),]
# fixes order
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm))
#ggdata$Lterm <- stringr::str_wrap(ggdata$Lterm, width = 50)  #Text on multiple lines if the name does not fit
#ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm))

#Minimalist graphics
gg1 <- ggplot(ggdata, aes(x = Lterm, y = -log10(Classic) ))+
  geom_point(size = 4, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  xlab('GO Term') +
  ylab('-log(p)') +
  labs(title = 'GO Biological Process')+
  theme_bw(base_size = 24) +
  coord_flip()
ggsave(paste0("Your_Folder", "Your_Grafic.png"), device = "png", width = 40, height = 30, dpi = 300, units = "cm")


#graph with more information
#library(stringr) If the name of the go term is too long, modify the # de Lterm

ggdata <- ggdata %>%
  mutate(
    log_p = -log10(Classic),
    gene_ratio = Significant / Annotated,
    #Lterm = str_wrap(Lterm, width = 60)
  )

gg1 <- ggplot(ggdata, aes(
  x = reorder(Lterm, log_p),  # Sort terms by significance
  y = log_p,
  color = log_p,              # Color by significance level
  size = Significant          # Size by number of genes
)) +
  geom_point(alpha = 0.8) +   # Transparency for overlays
  scale_color_viridis_c(
    name = "-log10(p)",
    option = "D",             
    guide = guide_colorbar(reverse = FALSE)
  ) +
  scale_size_continuous(
    name = "Nº of genes",
    range = c(4, 12),         # Range of dot sizes
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(
    title = "GO Enrichment: Molecular Funtion",
    x = "GO Term",
    y = "-log10(p-value)",
    caption = "Size: Number of significant genes\nColor: Significance level"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  coord_flip()

ggsave("Your_Folder/Your_Grafics.png", gg1, width = 40, height = 30, units = "cm", dpi = 300)
