########################################################################################
#                       Face-to-Face Contact Patterns in a Primary School
#Robin Beura
#
########################################################################################

# clear everything out of memory
rm(list=ls())

getwd()
# Set your working directory
# Save the data files to a location on your hard drive and specify the path here (Windows systems use forward slashes)
dir_path <-"C:\\Users\\robin\\Desktop\\Part 2"
setwd(dir_path)

## Load igraph package
library(igraph)

# load data from Stehle´ J, Voirin N, Barrat A, Cattuto C, Isella L, et al. (2011): doi:10.1371/journal.pone.0023176
infile_edges<-"Edges_sp_data_school_day_2.csv"
infile_nodes<-"Nodes_sp_data_school_day_2.csv"

###
edge_frame=read.csv(infile_edges, header = TRUE, sep = ",")
node_frame=read.csv(infile_nodes, header = TRUE, sep = ",")
g_primschool_orig<-graph.data.frame(edge_frame, directed = FALSE, vertices= node_frame)
head(node_frame)
head(edge_frame)

node_frame$label[node_frame$classname == 'Teachers']


med<-median(E(g_primschool_orig)$weight)

# This is the network that you will analyze: g_primschool_final
g_primschool_final<-delete.edges(g_primschool_orig, which(E(g_primschool_orig)$weight < med))

## For visualization purposes, the edge widths are set to the standardized value of weight.
E(g_primschool_final)$width <- (E(g_primschool_final)$weight - mean(E(g_primschool_final)$weight))/sd(E(g_primschool_final)$weight)

### This  happens to be the Fruchterman-Reingold, but you may choose any layout algorithm by changing the optional setting below
plot(g_primschool_final, layout=layout.fruchterman.reingold(g_primschool_final), 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=10, 
     vertex.label.color="black", 
     edge.width=0.5)

### Using with Kamada kawai layout
l = layout.kamada.kawai(g_primschool_final)
igraph.options(vertex.size=3, edge.arrow.size=0.5,vertex.label=NULL)
plot(g_primschool_final,edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=10, 
     vertex.label.color="black", 
     edge.width=0.5, layout =l)

##########           Basic network statistics

ecount(g_primschool_final)   # Edges for each node

vcount(g_primschool_final)   # nodes count

is.simple(g_primschool_final)  # True

is.connected(g_primschool_final)  #decides whether the graph is weakly or strongly connected.

components(g_primschool_final)  # Calculate the maximal (weakly or strongly) connected components of a graph

stgclusters <- clusters(g_primschool_final, mode="strong")$membership    #if strong cluster
plot(g_primschool_final, vertex.color = stgclusters, vertex.size=3, vertex.label= NA)

weakclusters <- clusters(g_primschool_final, mode="weak")$membership   #no members
plot(g_primschool_final, vertex.color = weakclusters, vertex.size=3, vertex.label= NA)  # can't plot because no members

reciprocity(g_primschool_final)  # measures the propensity of each edge to be a mutual edge

transitivity(g_primschool_final) # also known as clustering coefficient, measures that probability that adjacent nodes of a network are connected

mean_distance(g_primschool_final, directed=FALSE)  ## average number of edges between any two nodes

diameter(g_primschool_final, directed=FALSE, weights=NA) #length of the longest path (in number of edges) between two nodes
get_diameter(g_primschool_final, directed=FALSE, weights=NA)

### Community detection
#Used Walktrap algorithm because it finds communities through a series of short random walks
school_comm_fast <- cluster_walktrap(g_primschool_final, weights=E(g_primschool_final)$weight)
c.w <- membership(school_comm_fast)
plot(school_comm_fast,g_primschool_final, vertex.label= NA, vertex.size=2)

sort(degree(g_primschool_final))  ### Considered a measure of direct influence of each nodes
sort(strength(g_primschool_final)) ### weighted measure of degree that takes into account the number of edges that go from one node to another

list.vertex.attributes(g_primschool_final) #Listing vertex attributes
class_name = get.vertex.attribute(g_primschool_final, "classname")
table(c.w, class_name, useNA = c("no"))  ###Shows there is some interaction

### Cliques 
maximal.cliques(g_primschool_final)
clique.number(g_primschool_final)  ## Size of the largest clique = 8
table(sapply(maximal.cliques(g_primschool_final), length))
###


# plot to see represent largest clique 
clique <- maximal.cliques(g_primschool_final)

cliques_large <- largest.cliques(g_primschool_final)

cliques_18 <- c(cliques_large[[1]])#,cliques_large[[2]],cliques_large[[3]],cliques_large[[4]],cliques_large[[5]], cliques_large[[6]], cliques_large[[7]], cliques_large[[8]])
g2 <- induced.subgraph(graph=g_primschool_final,vids=(cliques_18))
plot(g2,main="cliques of size 18",vertex.size=14, edge.arrow = 0.2)

######## Betweenness
# Vertex betweenness centrality
head(sort(betweenness(g_primschool_final), decreasing = T), n = 10) ###Found betweenness and sorted based on the maximum paths they appear in
plot(g_primschool_final, layout=l, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=betweenness(g_primschool_final)/100, 
     vertex.label.color="black", 
     edge.width=0.5)

sort(edge_betweenness(g_primschool_final, e = E(g_primschool_final), directed = F, weights = g_primschool_final$weight), decreasing = T)
lay = layout_with_lgl(g_primschool_final)
plot(as.undirected(g_primschool_final), layout=lay, margin=c(-0.25,-0.25), vertex.size=betweenness(g_primschool_final)/100)

#edge.betweenness and weight
edgebetweens<-edge.betweenness(g_primschool_final, e=E(g_primschool_final), directed = F)
num_weight <- E(g_primschool_final)$weight
inv_weight <- 1/log(E(g_primschool_final)$weight  + 1)
edge_frame<-data.frame(edgebetweens, num_weight, inv_weight)
a_edge<-aggregate(edgebetweens ~ inv_weight, data=edge_frame, mean)  ## y~x : y variables are numeric data to be split into groups according to the grouping x variables (usually factors)
plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")

#edge transivity and weight
neigh <- neighbors(g_primschool_final, v=V(g_primschool_final)[V(g_primschool_final) == 1594])
overlap <- list()
valu <- list()

for (val in neigh$name){
        nh1 = neighborhood(graph = g_primschool_final, order = 1, nodes = '1594')[[1]]
        nh2 = neighborhood(graph = g_primschool_final, order = 1, nodes = val)[[1]]
        common = intersect(nh1, nh2)
        commonl = length(common)
        union = degree(g_primschool_final, v=c('1594')) + degree(g_primschool_final, v=val) - commonl - 2
        overlap <- append(overlap, list(commonl/union))
        valu <- append(valu, list(val))
}

#degree and transivity
clus_coe <- transitivity(g_primschool_final, type = "local")
ver_deg <- degree(g_primschool_final)
ver_df <- data.frame(clus_coe, ver_deg)
v_deg <- aggregate(clus_coe ~ ver_deg, data = ver_df ,mean)
plot(v_deg, col="blue", log="xy", xlab="vertices degree", ylab="Average transitivity")
plot(ver_deg, clus_coe, col = 4, xlab="degree", ylab="Cluster coefficient   /   transitivity")
