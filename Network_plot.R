
###############################################################################################
### This script creates a co-author network plot for Tina Lonsdorf: STATIC AND INTERACTIVE
# University Medical Center Hamburg-Eppendorf
# March, 2023
# written by Maren Klingelhöfer-Jens
# adjusted Code from https://datalab.ucdavis.edu/2019/08/27/creating-co-author-networks-in-r/
###############################################################################################

# Load packages
library(here)
library(scholar)
library(networkDynamic)
library(ndtv)
library(igraph)
library(statnet)
library(intergraph)
library(visNetwork)
library(dplyr)
library(stringr)

# Set a seed
set.seed(42)

# Get publications from Google Scholar
lonsdorf <- get_publications("UWLcKcIAAAAJ")

# Convert factor variables to strings
lonsdorf[,sapply(lonsdorf, class)=="factor"] = 
  as.character(lonsdorf[,sapply(lonsdorf, class)=="factor"])
saveRDS(object = lonsdorf, "lonsdorf_scholar_pull.RDS")
write.csv(lonsdorf, here("data/lonsdorf_publications.csv"), row.names = FALSE)

# The code above restrict the number of authors to 5, but we want all of them
ellipsis_indices <- grep(lonsdorf$author, pattern =  "\\.\\.\\.")
author_complete <- sapply(lonsdorf$pubid[ellipsis_indices], 
                         get_complete_authors, id = "UWLcKcIAAAAJ")
saveRDS(object = author_complete, file = "author_complete.RDS")
write.csv(lonsdorf, here("data/lonsdorf_publications_complete.csv"), row.names = FALSE)

### Convert names to a nice format
# Helper function to put last names in a regular format (first letter uppercase, rest lowercase)
lowerName <- function(x) {gsub(x, pattern = "\ ([A-Z]){1}([A-Z]+)", replacement = '\ \\1\\L\\2', perl = TRUE)}

# Helper function to convert a name to an initial (applied after lowername):
initialName <- function(x) {gsub(x, pattern = "([A-Za-z]){1}([a-z]+)", replacement = '\\1', perl = TRUE)}


# Format names
author_complete_reformat = lapply(author_complete, function(elem) {
  
  # Split strings of authors into individual authors based on commas:
  split_elem = strsplit(elem, ", ")[[1]]
  split_elem = sapply(split_elem, gsub, pattern = "(\\x.{2})+", replacement ="")
  
  # Put author names into INITIALS, Lastname format:
  rename_elem = sapply(split_elem, function(name) {
    
    #in case of name like "H J\xfc\xbe\x8d\x86\x94\xbcrvinen":
    name2 = iconv(name, "latin1", "ASCII", sub = "") 
    name2 = lowerName(name2)
    name2 = strsplit(name2, " ")[[1]]
    
    lastname = last(name2)
    
    if (length(name2) > 1) {
      name2 =  sapply(1:(length(name2)-1), function(i) {initialName(lowerName(name2[i]))})
      name2 = paste(paste(name2, collapse = ""), lastname, sep = " ")
    } else {
      name2 = lastname
    }
    return(name2)
  })
  
  # Put separated names back in single strings:
  rename_elem = paste(rename_elem, collapse = ", ")
  
  return(rename_elem)
})


### Add these names to the lonsdorf data frame
# Save original author column as "author_orig" and update the "author column"
lonsdorf$author_orig = lonsdorf$author
lonsdorf$author = as.character(lonsdorf$author)
lonsdorf$author[ellipsis_indices] = author_complete_reformat

# Convert all T Lonsdorf to TB Lonsdorf so that these entries refer to the same person
lonsdorf$author = sapply(lonsdorf$author, str_replace, pattern = "T Lonsdorf", replacement = "TB Lonsdorf")

# Save this data for further cleaning
write.csv(lonsdorf, here("data/lonsdorf_publications_2clean.csv"), row.names = FALSE)

###########################################
### Add cleaning steps in Excel!
###########################################

# Load data again
lonsdorf <- read.csv2("data/lonsdorf_publications_cleaned.csv", stringsAsFactors = FALSE, sep = ",")

# Exclude stread.csv2()# Exclude studies with sooo many authors
lonsdorf_static <- lonsdorf[-grep("A community-sourced", lonsdorf$title) , ]

# Creates a list equal in length to the number of publications,
# where each entry is a vector of length equal to the number of authors on a given publication
lonsdorf.coauthors <- sapply(as.character(lonsdorf_static$author), strsplit, ", ")

# Clean whitespace from all individual names
lonsdorf.coauthors <- lapply(lonsdorf.coauthors, trimws)

# create the node set (alphabetized)
lonsdorf.coauthors.unique <- unique(unlist(lonsdorf.coauthors))[order(unique(unlist(lonsdorf.coauthors)))]

# Create the edge set by tracking which authors appear in which papers
# These can be referred to as bipartite edges because entities (authors and papers) are only connected to entities of a different class.
# (Authors are connected to papers and vice versa.)

lonsdorf.bipartite.edges <- lapply(lonsdorf.coauthors, function(x) {lonsdorf.coauthors.unique %in% x})
lonsdorf.bipartite.edges <- do.call("cbind", lonsdorf.bipartite.edges) # dimension is number of authors x number of papers
rownames(lonsdorf.bipartite.edges) <- lonsdorf.coauthors.unique

# lonsdorf.bipartite.edges is a matrix of dimension equal to the number of papers by the number of unique authors in our data set
# Bipartiate edges (between authors and papers) can be projected into unimode edges (between authors) using matrix multiplication (%*%).

lonsdorf.mat <- lonsdorf.bipartite.edges %*% t(lonsdorf.bipartite.edges) #bipartite to unimode
mat <- lonsdorf.mat[order(rownames(lonsdorf.mat)), order(rownames(lonsdorf.mat))]

# lonsdorf.mat is a square matrix with a number of rows equal to the number of unique authors in our data set, 
# where the i,jth entry is equal to the numbers of co-authorships between author $i$ and author $j$.

# Create a network from lonsdorf.mat: edge values correspond to the number of co-authorships between authors 
lonsdorf.statnet = as.network(lonsdorf.mat, directed = FALSE, names.eval = "edge.lwd", ignore.eval = FALSE)
lonsdorf.statnet # view network summary

# Customize node color
col.1 <- adjustcolor("magenta1", alpha = 0.8)
col.2 <- adjustcolor("orchid4", alpha = 0.8)
node.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
node.col <- node.pal(50)


### Static visualization
# Add the degree of co-authoring
author.codegree = lonsdorf.mat["TB Lonsdorf",]
lonsdorf.statnet%v%"size" = log(author.codegree) + .5 #

# Plot and save it
#png("lonsdorf_network.png", width = 40, height = 40, units='cm', res = 300)
plot.network(lonsdorf.statnet, edge.col = "gray", edge.lwd = .2,
label = "vertex.names", vertex.col = node.col, vertex.border="#ffffff", label.cex = .7, label.pad = 0, label.pos = 1, vertex.cex = "size")
#dev.off()

#############################################################
### Interactive plot ###
#############################################################

# Create a data frame with nodes and edges for the graph
lonsdorf.nodes <- data.frame(id = 1:length(lonsdorf.statnet%v%"vertex.names"),
                             label = lonsdorf.statnet%v%"vertex.names",
                             title = lonsdorf.statnet%v%"vertex.names",
                             size = 5*(2+lonsdorf.statnet%v%"size"))

lonsdorf.edges <- data.frame(from=data.frame(as.edgelist(lonsdorf.statnet))$X1, 
                             to=data.frame(as.edgelist(lonsdorf.statnet))$X2)

# Plot it
lonsdorf_interactive = visNetwork(lonsdorf.nodes, lonsdorf.edges, main = "THE FABULOUS NETWORK OF TINA B. LONSDORF", width = 800, height = 800) %>% 
  visIgraphLayout(layout = "layout_nicely", type = "full")

# Customize nodes and edges
lonsdorf_interactive = lonsdorf_interactive  %>%
  visNodes(color = list(background = "white", border = "black", highlight = "#996699",
                        hover = list(background = "#996699", border = "magenta")),
           shadow = list(enabled = TRUE, size = 10),
           font = '30px') %>%
  visEdges(selectionWidth = 10, color = list(highlight = "magenta"), smooth = FALSE) %>%
  visInteraction(hover = TRUE)
lonsdorf_interactive # view interactive network

# Add an author dropdown menu
lonsdorf_interactive = lonsdorf_interactive  %>%  
  visOptions(nodesIdSelection = list(enabled  = TRUE, useLabels = TRUE, main = "Select by Author"))

# Save it
#visSave(lonsdorf_interactive, "lonsdorf_coauthor_network.html", selfcontained = T)


