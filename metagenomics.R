################################################################################
# Install packages (if needed)
################################################################################
install.packages("BiocManager",
                 "devtools",
                 "RColorBrewer",
                 "ape", 
                 "tidyverse",
                 "ggplot2",
                 "phytools",
                 "stringr", 
                 "reshape2",
                 "cowplot",
                 "grid",
                 "circlize",
                 "dplyr",
                 "devtools")

BiocManager::install("pathview")

library(devtools)
install_github("YuLab-SMU/ggtree")

################################################################################
# Load packages
################################################################################


# Working with trees
library(ape)
library(phytools)
# Plotting Grammar of Graphics
library(ggplot2)
# Tree visualization
library(ggtree)
# String maniputation
library(stringr)
# Data structure manipulation
library(reshape2)
# Multi-panel figure layouts
library(cowplot)
library(grid)
# Color maniputation
library(RColorBrewer)
# Pathway analysis with KEGG Orthologs
library(pathview)
#plot circular graphs and contig maps
library(circlize)
# Suite of packages for data manipulation and visualization
library(tidyverse)


################################################################################
# Downloading data
################################################################################

#All data can be found in the directory 
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/checkm.qa.tsv",
"data/checkm.qa.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/gtdb.classification.tree",
"data/gtdb.classification.tree")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/gtdb.classification.tsv",
"data/gtdb.classification.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/gtdb.taxonomy.tsv",
"data/gtdb.taxonomy.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/metagenome.counts.tsv",
"data/metagenome.counts.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/KO.genomes.tsv",
"data/KO.genomes.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/KO.pathways.tsv",
"data/KO.pathways.tsv")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/LE.CH5.EKL.cov",
"data/LE.CH5.EKL.cov")
download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/data/LE.CH5.gff",
"data/LE.CH5.gff")


################################################################################
# Downloading custom functions
################################################################################

download.file(
"https://raw.githubusercontent.com/davidlevybooth/Metagenomics_Workshop/master/scripts/metagenome.utilities.R",
"scripts/metagenome.utilities.R")

#Loading custom functions for metagenomics and metatranscriptomics
source("scripts/metagenome.utilities.R")

################################################################################
# Section 1. Basic Tree Plotting
################################################################################

################################################################################
# Loading GTDB Tree
################################################################################

tree.path <- "data/gtdb.classification.tree"

#Use read.tree(<tree.path>) to load tree file from path.
#Important: Make sure you use the "ape" version of the read.tree() function
tree <- ape::read.tree(tree.path)

#Inspect tree object
tree


################################################################################
# Loading GTDB Taxonomy
################################################################################

#Set path
GTDB.taxonomy.path <- "data/gtdb.taxonomy.tsv"

#Load data
GTDB.taxonomy <- read.csv(GTDB.taxonomy.path, sep="\t", stringsAsFactors = FALSE)

#Inspect data
head(GTDB.taxonomy)

#Isolate and inspect tip labels. Each tip is a genomic placement. 
tips <- tree$tip.label
head(tips)

#Randomly take a sample of 1000 tips 
#Let's set "seed" to ensure that we all end up with the same list of tips to keep
set.seed(69420)
tips.keep <- sample_gtdb_tips(tips, 1000)

GTDB.taxonomy <- GTDB.taxonomy[which(GTDB.taxonomy$GenomeID %in% tips.keep), ]
#Use keep.tips() to subset tree.

tree.subset <- ape::keep.tip(tree, tips.keep)

#if you are having difficulty with you subset tree, please uncomment to load from the following file:
#tree.subset <- read.tree("results/tree.subset")


################################################################################
# Plot the reduced tree with ggtree
################################################################################

#Plot the tree using rectangular layout
tree.plot <- ggtree(tree.subset) + xlim(0,4)
tree.plot

#Plot the tree using ciculat layout
ggtree(tree.subset, layout="circular")


#add tip labels with geom_tiplab()
tree.plot + geom_tiplab(size = 2)

################################################################################
# Add phyla labels with geom_cladelabel()
################################################################################

#Assign our phylogenetic clade of interest to a variable called "Clade"
Clade <- "p__Firmicutes_A"

#Here we combine a few base R functions to subset the GTDB taxonomy table to only: 
GTDB.taxonomy.clade <- GTDB.taxonomy[which(GTDB.taxonomy[,"Phylum"] == Clade),]
Firmicutes_A_tips <- GTDB.taxonomy.clade$GenomeID

#Inspect tips
head(Firmicutes_A_tips)


#find the common ancestor nodes for a list of tree tips. 
Firmicutes_A_node <- findMRCA(tree, Firmicutes_A_tips, type="node")
Firmicutes_A_node

#Use custom function: collapse_nodes() to apply findMRCA() over the whole tree to collect nodes for all phyla.
Phyla.nodes <- collapse_nodes("Phylum", tree.subset, GTDB.taxonomy)


#Almost done, next we can highlight clades of interest. 
#First, select Firmicutes node again from our handy dataframe.
Firmicutes_A_node <- Phyla.nodes[which(Phyla.nodes$Group == "p__Firmicutes_A"), ]$Node
Firmicutes_A_node

################################################################################
# Plot the tree with ggtree
################################################################################

tree.plot <- ggtree(tree.subset) + xlim(0,4)
tree.plot + geom_balance(node=Firmicutes_A_node, 
                         fill='darkgreen', 
                         color='white', 
                         alpha=0.6, 
                         extend=1)  +
  geom_cladelabel(node = Firmicutes_A_node, 
                  "Firmicutes A", 
                  offset = 1, 
                  fontsize = 2)

#We can run geom_balance() to produce colored ranges for all phyla in our dataframe in a loop. 
#Select colors
colors <- some_color(Phyla.nodes)

Phyla.nodes <- collapse_nodes("Phylum", tree.subset, GTDB.taxonomy)

for(i in 1:nrow(Phyla.nodes)) {
  if(Phyla.nodes$TipCount[i] > 10) {
    tree.plot <- tree.plot + 
      geom_balance(node=Phyla.nodes$Node[i], fill=colors[i], color='white', alpha=0.6, extend=0.6)  +
      geom_cladelabel(node = Phyla.nodes$Node[i], Phyla.nodes$Group[i], offset = 0.6, fontsize = 2)
  }
}
tree.plot


#Can we do this for a circular plot?
tree.plot <- ggtree(tree.subset, layout="circular")
tree.plot <- clade_colors(tree.plot, Phyla.nodes, colors)
tree.plot

#(Looks pretty badass to me!)

#One more thing: adding bootstrap values to trees
tree.subset$node.label <- parse_bootstraps(tree.subset, method = "parse")

#Find and parse tree boostraps
bs_count <- parse_bootstraps(tree.subset, method = "count")
tail(bs_count)

#Plot tree with bootstraps as text
ggtree(tree.subset) + xlim(0,3) + geom_text(aes(label = bs_count), size = 2)

#Plot tree with bootstraps as nodepoint
tree.plot <- ggtree(tree.subset) + 
  xlim(0,4) + 
  geom_nodepoint(aes(subset= bs_count >= 90), 
                 fill = "cadetblue",
                 size=1.5, 
                 alpha = 0.5, 
                 shape = 21)

tree.plot


# And add phyla labels:
tree.plot <- clade_labels(tree.plot, Phyla.nodes, tiplimit = 5, fontsize = 2)
tree.plot


#add a scale representing the substitution distance

tree.plot + geom_treescale(x=0, y=length(tree.subset$tip.label)-50, width=0.2, offset = 10)


################################################################################
# Section 2. Plotting MAG phylogenies
################################################################################

#Set path for MAG taxonomy table and import. 
mag.taxonomy.path <- "data/gtdb.classification.tsv"
mag.taxonomy <- read.csv(mag.taxonomy.path, sep="\t")

#Subset tree to show only MAG placements using keep.tips().
tree.mags <- keep.tip(tree, as.character(mag.taxonomy$GenomeID))
tree.mags

#Plot MAG tree with ggtree()
tree.plot <- ggtree(tree.mags, ladderize=F) + geom_tiplab(size = 2) + xlim(0,2.6)
tree.plot

#Add arbitrary color ranges with the predefined some_color() function. 
Phyla.nodes <- collapse_nodes("Phylum", tree.mags, mag.taxonomy)
colors <- some_color(Phyla.nodes)

#Create tree plot with custom clade colors: 
tree.plot <- clade_colors(tree.plot, Phyla.nodes, colors, extend = 0.6)

#Add bootstraps with the parse_bootstraps() function: 
tree.mags$node.label <- parse_bootstraps(tree.mags, method = "parse")
bs_count <- parse_bootstraps(tree.mags,method = "count")

#Finally: Add it all together for a great looking phylogeny of our MAGs:
tree.plot <- tree.plot + 
geom_nodepoint(aes(subset= bs_count >= 75), fill = "cadetblue", size=2, alpha = 0.5, shape = 21) + 
geom_treescale(x=0, y=1, width=0.1, offset = 0.5)

tree.plot


################################################################################
# Section 3. MAG statstic plots
################################################################################

#Load MAG statistics data from CheckM
MAG.checkM.path <- "data/checkm.qa.tsv"
MAG.checkM <- read.csv(MAG.checkM.path, sep="\t")

#Inspect data
head(MAG.checkM)

#To plot data we must organize data in the 'long' format using reshape2::melt()
MAG.checkM.long <- melt(MAG.checkM, id.vars = "GenomeID")
head(MAG.checkM.long)


#Set GenomeID as an ordered factor to match tree tip order. 
MAG.checkM.long$GenomeID <- factor(MAG.checkM.long$GenomeID, levels = tree.mags$tip.label)
MAG.checkM.long$variable <- factor(MAG.checkM.long$variable, levels = rev(unique(MAG.checkM.long$variable)))

#Create separate 'completeness' data.frame (still in "long" format). 
MAG.checkM.long.completeness <- MAG.checkM.long[grepl("Completeness", MAG.checkM.long$variable), ]
MAG.checkM.long.completeness$variable <- factorize(MAG.checkM.long.completeness$variable, rev = TRUE)

#Plot basic completeness plot
ggplot(MAG.checkM.long.completeness,
  aes(x = GenomeID, y = value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + #The "black" color provides the border. 
  scale_fill_manual(values = c("cadetblue", "gold")) + #Manually assign fill for residual and completeness
  coord_flip() + #Switches the x and y axis to produce vertical stacked bar plot. 
  theme(legend.position = "none") + #We don't need a legend for these data
  ylab("Completeness") #Add label

#Final completeness plot
MAG.completeness <- tree_bar_plot(MAG.checkM.long.completeness, "Completeness")
MAG.completeness

#Parse data for, and plot, additional statistics plots
MAG.checkM.long.contamination <- MAG.checkM.long[grepl("Contamination", MAG.checkM.long$variable), ]
MAG.checkM.long.contamination$variable <- factorize(MAG.checkM.long.contamination$variable, rev = TRUE)

MAG.checkM.long.coverage <- MAG.checkM.long[grepl("Coverage", MAG.checkM.long$variable), ]
MAG.checkM.long.coverage$variable <- factorize(MAG.checkM.long.coverage$variable, rev = TRUE)


MAG.contamination <- tree_bar_plot(MAG.checkM.long.contamination, "Contamination", 
                                   colors = c("cadetblue", "firebrick"))

MAG.coverage <- tree_bar_plot(MAG.checkM.long.coverage, "Coverage", 
                              colors = c("cadetblue", "darkblue"))
                              
#Put it all together to plot genome statistics with the tree we created earlier: 
tree.plot <- tree.plot + scale_y_tree()

MAG.title <- plot_title("Liard River Phylogeny and Statistics")

MAG.plots <- plot_grid(tree.plot, 
                       MAG.completeness, 
                       MAG.contamination, 
                       MAG.coverage, 
                       align = "h", 
                       axis = "b", 
                       ncol = 4, 
                       rel_widths = c(1.25, 0.25, 0.25, 0.25))

#Nesting the plot_grid() function allows us to easily add titles to our groups of plots!
plot_grid(MAG.title, 
          MAG.plots, 
          ncol = 1, 
          rel_heights = c(0.1, 1))

#Looks great!
  

################################################################################
# Section 4. Aligning abundance plots with phylogenetic trees.
################################################################################

#use a heatmap to explore MAG environmental abundance and temperature range.

#Load count data.
counts.path <- "data/metagenome.counts.tsv"
counts <- read.table(counts.path, sep="\t", header = TRUE)

#parse the lowest level of taxonomic classification 
counts$Classification <- tax_class(counts$Taxonomy)

#Optionally, we can calculate relative abundance for each sample:
counts.ra <- relative_abundance(counts)

#Or normalize across taxa:
counts.norm <- normalize(counts)

#We once again match the counts.norm and mags.taxonomy.norm data.frames with the order of tips in the tree.
#It can be confusing that the first vector in the match function is actually
#what you're trying to match to, so be careful with the order. 
counts.norm <- counts.norm[match(tree.mags$tip.label, counts.norm$GenomeID), ]
mag.taxonomy.norm <- mag.taxonomy[match(tree.mags$tip.label, mag.taxonomy$GenomeID), ]

#Now that the tree tip, counts.norm and mag.taxonomy.norm all match we can consolidat the variables we want to retain in the counts.norm data.frame. 
counts.norm$GenomeID <- tree.mags$tip.label
counts.norm$Taxonomy <- mag.taxonomy.norm$Taxonomy

#Finally we use the tax_class() function to pull out the lowest level of GTDB classification for each genome in the counts.norm data.frame. 
counts.norm$Classification <- tax_class(counts.norm$Taxonomy)

#Melt data.frame for plotting and do some more formatting:
#Transform the count data to the "long" format for plotting. 
#reshape2::melt() creates three types of columns: IDs, variables and values. 

#We want to retain "GenomeID", "Taxonomy", and "Classification" as IDs. 
#And we keep the actual counts as the "value". Making each sample a "variable"
counts.long <- melt(counts.norm, 
                    id.vars = c("GenomeID", "Taxonomy", "Classification"), 
                    value.name = "norm.count")

#Use factorize() again to lock in the current variable order when plotting
counts.long$variable <- factorize(counts.long$variable)

#We use the paste0 function to create a new label by combining the genome name and 
#the classification. 
counts.long$Label <- paste0(counts.long$Genome, " (", counts.long$Classification, ")")
counts.long$Label <- factorize(counts.long$Label)

#Finally, we use cut_last(), which is a wrapper for the substr() function. 
#This allows us to trim characters from the end of our sample variables 
#removing the replicates and keeping the pool IDs. 
counts.long$Pool <- cut_last(counts.long$variable, 3)
counts.long$Pool <- factorize(counts.long$Pool)

#Let's take a look:
head(counts.long)

#Plot Gene Count Heatmap with ggplot2 using geom_tile()

counts.heatmap <- ggplot(data = counts.long, mapping = aes(x = variable, y = Label, fill = norm.count)) +
  geom_tile(color = "grey50") +
  scale_fill_gradient(low = "#FFFFFF", high = "#012345", na.value = "#FFFFFF") +
  scale_y_discrete(drop=FALSE) + #keep even "empty" genomes to preserve all elements in the tree
  theme(axis.text = element_text(size = 8, vjust = 0.4),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank()) +
  facet_grid(.~Pool, scales = "free_x")

#Let's take a look at the plot!
counts.heatmap

#Finally, combine tree and heatmap plot with cowplot::plot_grid()
counts.long$Genome <- factor(counts.long$GenomeID, levels = tree.mags$tip.label)

#Remove y-axis labels
new.count.heatmaps <- heatmap_plot(counts.long) + theme(axis.text.y = element_blank())

#Put them together with cowplot::plot_grid()
MAG.title <- plot_title("Liard River MAG Phylogeny and Abundance")

MAG.plots <- plot_grid(tree.plot, 
                       new.count.heatmaps, 
                       align = "h", 
                       axis = "bt", 
                       ncol = 2, 
                       rel_widths = c(1, 1))

plot_grid(MAG.title, 
          MAG.plots, 
          ncol = 1, 
          rel_heights = c(0.1, 1))


################################################################################
# Section 5. Functional Annotation 
################################################################################

# Load KEGG pathways of interest
KO.pathways <- read.table("data/KO.pathways.tsv", sep = "\t", header = TRUE)
KO.pathways$Function <- factorize(KO.pathways$Function) 
#This is important to keep the pathway data organized for later!

#Load MAG kofamscan annotations
KO.mag <- read.delim("data/KO.genomes.tsv", sep = "\t", header = TRUE)
head(KO.mag)

#Select only genomes KOs in pathways of interest:
# N-cycling
# S-cycling
# CH4-cycling
# Aromatic degradation
# Photosynthesis
# Etc. 

#use match again to ensure that the MAG and Pathway data are in the same order. 
KO.mag <- KO.mag[which(KO.mag$KO %in% KO.pathways$KO), ]
KO.mag$Function <- KO.pathways[match(KO.mag$KO, KO.pathways$KO), ]$Function
KO.mag$Pathway <- KO.pathways[match(KO.mag$KO, KO.pathways$KO), ]$Pathway

#Select pathways to plot: 
keep.pathways <- c("Nitrogen Metabolism", "Sulphur Metabolism", "Photosynthesis", "Hydrogen Metabolism", "Carbon Metabolism")

#Isolate selected pathways in data.frame
KO.mag <- KO.mag[which(KO.mag$Pathway %in% keep.pathways), ]

#Order the genomes for plotting based on order of tips in tree using a factor:
KO.mag$Genome <- factor(KO.mag$GenomeID, levels = as.character(tree.mags$tip.label))

#We also apply facetting to divide the individual functions. 
pathway.colors <- set_path_colors()


#Plot functional annotations with ggplot2
function.plot <- ggplot(KO.mag, aes(x = KO, y = Genome, fill = Function)) +
  geom_tile(color = "black") + #geom_tile creates each square for an annotated gene
  scale_fill_manual(values = pathway.colors) + #I use a named vector to apply fills to specific conditions
  theme_bw() + # this is just making it pretty
  theme(axis.title = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, size = 6),
        legend.position = "none", 
        strip.text.x = element_text(size = 8, angle = 90, hjust = 0, vjust = 0.6), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  scale_x_discrete(position = "top") +
  facet_grid(.~Function, scales="free_x", space="free_x") #facetting by Function 

function.plot


################################################################################
# Section 6. Pathway Analysis 
################################################################################

#Optionally, we can use the _pathview_ package to create Pathway Analysis maps. 
#We will look at the benzoate degradation map for the Chloroflexi MAG L.E.CH.5. 

#Create labelled benzoate degradation pathway map for genome: CH.5
CH5 <- KO.mag[which(KO.mag$GenomeID == "L.E.CH.5"), ]

pathview(gene.data = as.character(CH5$KO),
  pathway.id = "00362", 
  species = "ko",
  gene.idtype = "KEGG",
  kegg.native = TRUE,
  out.suffix = "CH5")


################################################################################
# Section 7. Mapping Metatranscriptomic Data
################################################################################


#The following data comes from incubating Epsilon sediment at 45^o^C with the following substrates: 
# 0.1% Kraft Lignin (EKL)
# 0.1% Milled Wood Lignin (EMWL)
# 0.1% DHP Lignin (DHP)
# 0.1% Vanillate (VAN)
# No exogenous carbon controls (NC)


#Load data file already in the "long" format
metatrans.long <- read.table("data/metatranscriptome.counts.long.tsv", sep = "\t")

#factorize the labels to maintain order on the y-axis
metatrans.long$Label <- factorize(metatrans.long$Label)

#Inspect
head(metatrans.long)

#Plot as a heatmap
trans.heatmap <- ggplot(data = metatrans.long, mapping = aes(x = variable, y = Label, fill = value)) +
  geom_tile(color = "grey50") +
  scale_fill_gradient(low = "#FFFFFF", high = "#012345", na.value = "#FFFFFF") +
  scale_y_discrete(drop=FALSE) +
  theme(axis.text = element_text(size = 8, vjust = 0.4),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank()) +
  facet_grid(.~Substrate, scales = "free_x")
trans.heatmap


################################################################################
# Section 8. Circular genome plots with gene expression data 
################################################################################

# Plot the metatranscriptomics data for a single genome bin. 
#E.g., the Chloroflexi MAG L.E.CH.5 grown on Kraft Lignin. 

#Prepare the data: 
#Load the coverage table generated with BBMap
CH5.EKL.COV <- read.table("data/LE.CH5.EKL.cov", sep = "\t", header = TRUE)
#calculate maxium expression value
cmax <- max(CH5.EKL.COV$cov)

#Load the genome features file (gff) for L.E.CH.5. 
# We have already converted the gff to a tab-delimited table for easy loading. 
CH5.gff <- read.table("data/LE.CH5.gff", sep = "\t", header = TRUE)

#Turn the genome name into a factor for circlize plotting to work 
CH5.gff$seqnames <- factor(CH5.gff$seqnames)

#Use circlize to plot genome map
#We begin bit initializing the plotting parameters
circos.par(start.degree = 0, track.margin = c(0,0), cell.padding = c(0,0), gap.degree = 5)
circos.initializeWithIdeogram(CH5.gff, plotType = c("axis", "labels"))

#Plotting Gene Expression
circos.genomicTrackPlotRegion(CH5.EKL.COV, stack = FALSE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, ytop.column = 2, ybottom = 0, col = "#F24405", border = "#F24405")
}, track.height = 0.2, bg.border = "grey80", ylim = c(0,cmax))

#Coding Regions
circos.genomicTrackPlotRegion(CH5.gff, stack = FALSE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value,
                     ytop = value[,2], ybottom = value[,3], #Plots ORF positions by strand 
                     col = "grey40", border = NA)
}, track.height = 0.1, bg.border = "grey50", ylim = c(-1,1))

circos.clear()

# That's it! 