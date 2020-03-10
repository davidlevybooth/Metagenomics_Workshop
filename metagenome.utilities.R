#Metagenome Utilities v0.1
#David Levy-Booth, Feb 2020
#These functions are designed to help with the "Metagenomic Data Visualization" Workshop
########################################################################################


#sample_gtdb_tips() v0.1
#---------------------------------------------------------------------------------------
#Take a list of GTDB-Tk tips and sample number 
#Seperate GTDB accessions from placed genomes
#Account for # of placed genomes when sampling GTDB accessions

sample_gtdb_tips <- function(tree.tips, sample_n) {
  if(class(tree.tips) != "character") { 
    stop(paste("Error: Input must be character vector"))
  }
  else {
    #Use grep with or statement "|" to isolate placed geneomes
    #from GTDB accesssions. Recall that "-" reverses grep boolean
    placements <- tree.tips[-grep("GB_|RS_", tree.tips)]
    sample_n <- sample_n-length(placements)
    gtdb_genomes <- tree.tips[grep("GB_|RS_", tree.tips)]
    #randomly sample GTDB genomes
    gtdb_genomes <- sample(gtdb_genomes, sample_n) 
    sampled_tips <- c(placements, gtdb_genomes) #recombine
    return(sampled_tips)
  }
}

#collapse_nodes() v0.1
#---------------------------------------------------------------------------------------
#Taking a phylogenetic tree, a data.frame of taxonomy strings and a taxonomic level
#Select tips in each taxonomic group, then calculate most recent common ancestor 
#(MRCA) node and # of tips in each node. Return data.frame.

collapse_nodes <- function(level, tree, taxonomy) {
  require(ape)  #requires tree structure from ape
  require(phytools)  #requires the findMRCA function from phytools
  nodelist<-list() #create empty vectors
  nodedist<-list()
  grouplist<-list()
  groupcount<-list()
  i <- 1
  #Select unique taxonomic levels
  #Ensure data.frame entries only include tip names
  taxonomy <- taxonomy[which(taxonomy[,1] %in% tree$tip.label), ]
  groups <- as.character(unique(taxonomy[,level]))

  for(group in groups) {
    print(group)
    #For each group, find all tips
    group_tips <- as.character(taxonomy[which(taxonomy[,level] == group),]$GenomeID)
    #Stuff will break if you call MRCA on groups with only 1 tip.  
    if(length(group_tips) > 1) {
      node <- findMRCA(tree, group_tips, type="node")
      provtree <- ape::keep.tip(tree, group_tips)
      justthetip <- provtree$edge.length
      groupcount[i] <- length(group_tips)
      grouplist[i] <- group
      nodelist[i] <- node
      nodedist[i] <- max(justthetip)
      i <- i + 1
    }
  }
  #Unlisting a list of list is a hacky way to mimic the python "append" function in R. 
  return(data.frame(Group = unlist(grouplist), TipCount = unlist(groupcount), Node = unlist(nodelist), Distance = unlist(nodedist)))
}


#clade_labels() v0.1
#---------------------------------------------------------------------------------------
#Takes a ggtree tree object and a node taxonomy data.frame. For each taxonomic group
#add label to tree. Additional variables control label placement. Returns labelled tree. 

clade_labels <- function(tree_plot, node_table, offset = 0.3, fontsize = 4, tiplimit = 10) {
  if("Node" %in% colnames(node_table) & 
     "Group" %in% colnames(node_table) & 
     "TipCount" %in% colnames(node_table)) {
    for(i in 1:nrow(node_table)) {
      if(Phyla.nodes$TipCount[i] > tiplimit)
        tree_plot <- tree_plot + geom_cladelabel(node = node_table$Node[i],
                                                 node_table$Group[i],
                                                 offset = offset,
                                                 fontsize = fontsize)
    }
    return(tree_plot)
  }
  else {
    stop(paste("ERROR: Ensure Node Table contains columns: Node, Group & TipCount"))
  }
}


#clade_colors() v0.1 
#Take a ggtree object, taxonomy data.frame and color vector (each entry == one 
#taxonomic group). Iterate through taxonomic groups and add color range to plot. 
#Returns ggtree plot with color ranges. 
#---------------------------------------------------------------------------------------

clade_colors <- function(tree_plot, node_table, colors, alpha = 0.5, extend = 0.5, tiplimit = 1) {
  if("Node" %in% colnames(node_table) & 
     "Group" %in% colnames(node_table) & 
     "TipCount" %in% colnames(node_table)) {
        for(i in 1:nrow(node_table)) {
          if(node_table$TipCount[i] > tiplimit) {
            tree_plot <- tree_plot + 
              geom_balance(node=node_table$Node[i], fill=colors[i], color='white', alpha=alpha, extend=extend)
          }
        }
    return(tree_plot)
  }
  else {
    stop(paste("ERROR: Ensure Node Table contains columns: Node, Group & TipCount"))
  }
}


#parse_bootstraps() v0.1
#---------------------------------------------------------------------------------------
#GTDB phylogenetic tree node labels can contain both MRCA labels as well as bootstrap
#decimal value. First we determine what information the node label contains, then if
#necissary, we extract the bootstrap value. Values are returned as integers from 1-100.  
#When using the "count" method this function returns a vector of boostraps with NA values
#for tip labels. Because ggtree is weird like that. 

parse_bootstraps <- function(tree, method = "parse") {
  require(stringr)
  if(method == "parse") {
    #ifelse is one of the most useful functions in R, here we combine it with the boolean 
    #vector version of grep (grepl) to act only on node labels with an appended MRCA string
    #Recall that the normal grep returns index vectors, which can't be used within ifelse()
    tree.bootstraps <- ifelse(grepl(":", tree$node.label), str_split_fixed(tree$node.label, ":", 2)[,1], tree$node.label)
    tree.bootstraps <- round(as.numeric(tree.bootstraps)*100,0)
    return(tree.bootstraps) 
  }
  else if(method == "count")
  {
    tl <- length(tree$tip.label)
    bs <- as.numeric(c(rep(NA, tl), tree$node.label))
    #bs[is.na(bs)] <- 0
    return(bs)
  }
  else { stop(paste("ERROR: Unsupported method."))}
}


some_color <- function(nodes, colors = c("darkgreen", "steelblue")) {
  return(rep(colors, round(nrow(nodes)/2,0),length.out = nrow(nodes)))
}


#tree_bar_plot() v0.1
#---------------------------------------------------------------------------------------
#Creates a nicely formatted two-color stacked barplot from a melted data.frame. 
#Can change x label and colors. Recall that b/c of coord_flip(), x and y axis are
#switched. Default "variable" and "value" column names are used. 

tree_bar_plot <- function(long.df, label, colors = c("cadetblue", "gold")) {
  p <- ggplot(long.df, 
         aes(x = GenomeID, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.5) + 
    theme(legend.position = "none", #Don't need a legend b/c of y-label
          axis.title.y = element_blank(), #remove y-axis title
          axis.text.y = element_blank(), #remove y-axis title (GenomeIDs)
          panel.background = element_blank(), #remove ugly default ggplot2 background
          axis.line.x = element_line(size = 0.5, color = "black")) + #add x-axis line
    coord_flip() + # Switch axis of horizontal barplot to create vertical plot
    scale_fill_manual(values = colors) + #apply color vector. 
    #This squishes the default x-axis spacing that ggplot likes to use. 
    scale_y_continuous(limits=c(min(long.df$value)-1, max(long.df$value)+1), expand = c(0, 0)) +
    ylab(label)
  return(p)
}


#scale_y_tree v0.2
#from: https://thackl.github.io/ggtree-composite-plots
#---------------------------------------------------------------------------------------
#Change the second element of tupple to match vertical position when combining
#with other plots using cowplot::plot_grid()


scale_y_tree <- function(expand=expand_scale(0.01, 0.1), ...){
  scale_y_continuous(expand=expand, ...)
} 

scale_y_tree2 <- function(expand=expand_scale(0.01, -0.01), ...){
  scale_y_continuous(expand=expand, ...)
}

#plot_title v0.2
#---------------------------------------------------------------------------------------
#Create plot title element for use in cowplot::plot_grid()

plot_title <- function(title_text) {
  title <- ggdraw() + 
    draw_label(
      title_text,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  return(title)
}


relative_abundance <- function(count_table) {
  nums <- unlist(lapply(count_table, is.numeric)) 
  if(all(nums) == TRUE) {
    count_table.t <- t(count_table)
    count_table.t.ra <- sweep(count_table.t, 1, rowSums(count_table.t), '/')
    count_table.ra <- data.frame(t(count_table.t.ra))
  }
  else {
    count_table.num <- count_table[ , nums]
    count_table.nonnum <- count_table[ , !nums]
    count_table.t <- t(count_table.num)
    count_table.t.ra <- sweep(count_table.t, 1, rowSums(count_table.t), '/')
    count_table.ra <- t(count_table.t.ra)
    count_table.ra <- cbind(count_table.nonnum, count_table.ra)
  }
  return(count_table.ra)
}

normalize <- function(count_table) {
  nums <- unlist(lapply(count_table, is.numeric))
  if(all(nums) == TRUE) {
    count_table.norm <- data.frame(t(apply(count_table, 1, function(x)(x-min(x))/(max(x)-min(x)))))
  }
  else {
    count_table.num <- count_table[ , nums]
    count_table.nonnum <- count_table[ , !nums]
    count_table.norm <- t(apply(count_table.num, 1, function(x)(x-min(x))/(max(x)-min(x))))
    count_table.norm <- cbind(count_table.nonnum, count_table.norm)
  }
  return(count_table.norm)
}



tax_class <- function(taxonomy) {
  empties <- c("p__", "c__", "o__", "f__", "g__", "s__")
  taxonomy_list <- rep(NA, length(taxonomy))
  i = 1
  for(tax in taxonomy) {
    spl <- strsplit(as.character(tax), ";")
    spl <- setdiff(spl[[1]], empties)
    cl <- tail(spl, n=1)
    taxonomy_list[i] <- cl
    i = i + 1
  }
  return(taxonomy_list)
}



cut_last <- function(char_list, n) {
  char_list <- substr(as.character(char_list),1,nchar(as.character(char_list))-n)
  return(char_list)
}


factorize <- function(char_list, rev = FALSE) {
  if(rev) {
    return(factor(char_list, levels = rev(unique(char_list))))
  }
  else { return(factor(char_list, levels = unique(char_list))) }
  
}


heatmap_plot <- function(heatmap_data.long) {
  counts.heatmap <- ggplot(data = counts.long, mapping = aes(x = variable, y = Genome, fill = norm.count)) +
    geom_tile(color = "grey50") +
    scale_fill_gradient(low = "#FFFFFF", high = "#012345", na.value = "#FFFFFF") +
    theme(axis.text = element_text(size = 8, vjust = 0.4),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank()) +
    facet_grid(.~Pool, scales = "free_x")
  return(counts.heatmap)
}

'%!in%' <- Negate('%in%')

set_path_colors <- function() {
  colvec <- c("#005824","#005824","#005824","#005824","#005824", "#005824", "#005824", #N
           "#7a0177", "#7a0177", "#7a0177", #S
           "#0c2c84", "#0c2c84", #CH4
           "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", "#cc4c02", #Aromatic
           "#FA0000", #halogenated
           "#dd1c77", "#dd1c77", "#dd1c77", "#dd1c77", #Photosynthesis
           "#a6761d", "#a6761d", #H
           "#1b9e77", "#1b9e77","#1b9e77","#1b9e77","#1b9e77","#1b9e77","#1b9e77", #CO2/CO
           "#e6ab02") #Fe
  funvec <- c("N Fixation", "Urea", "NH3 Oxidation", "NO2 Oxidation", "NO3 Reduction", "Denitrification", "ANAMMOX",
  "SO4 Reduction", "H2S Oxidation", "S2O3 Oxidation", 
  "CH4", "CH4 Oxidation",          
  "Phenol", "Methylphenol", "Catechol Metacleavage", "Catechol Orthocleavage", "Vanillate", "p-HP Degradation", "PCA Orthocleavage", "bKA Pathway", "PCA Metacleavage",                   
  "Halogenated",
  "Photosystem I", "Photosystem II", "Anox. Photosystem I",  "Anox. Photosystem II",
  "H2 (FeFe)", "H2 (NiFe)", 
  "CO2", "CO", "CBB Cycle", "Reverse TCA", "W-L pathway", "3HP Bi-cycle", "3HP-4HB Cycle",                       
  "Fe(II)")
  names(colvec) <- funvec
  return(colvec)
}
