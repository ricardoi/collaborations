library(bipartite)
library(igraph)

file <- read.csv(file.choose(), row.names = 1)
dat.mat = as.matrix(file[1:6,])
plotweb(sortweb(dat.mat, sort.order="inc"), method="normal")   
visweb(sortweb(dat.mat,sort.order="inc"), type="diagonal", labsize=1,
       square="interaction", text="none", textsize = 1,circles=FALSE, frame=FALSE)# 
# bipartite 
dimdatb.net.table <- networklevel(dat.mat)
datb.sp.table  <- specieslevel(dat.mat)
datb.nest.table <- nestedcontribution(dat.mat)

bacteria <-  graph.incidence(dat.mat)
plot(bacteria)
#@ Checking type of network and names
V(bacteria)$type
V(bacteria)$name
#@ Generating attributes: colors and shapes
cols= c("#DFBCAE", "#B4DDE5", "#E485C6", "#FCC63C",  "#F5FBC5", "#475BDE", "#B8B8B6", "#50E8CC", "#9BD198", "#9F5FD3", "#579BC8", "#E77375")
shapes = c(rep("circle", length(V(bacteria)$type[V(bacteria)$type == "FALSE"])), rep("square", length(V(bacteria)$type[V(bacteria)$type == "TRUE"])))
#@ Adding attributes
V(bacteria)$color <-  c(rep("#F4A460", length(V(bacteria)$type[V(bacteria)$type == "FALSE"])), cols)
V(bacteria)$size <- c(datb.sp.table$`lower level`$degree, (datb.sp.table$`higher level`$weighted.betweenness+1)^5, 0.5)+5

#@ plot bipartite network
plot(bacteria, edge.arrow.size=1, vertex.shape=shapes, vertex.label.cex=1, vertex.label.color='black', vertex.frame.color="gray", 
     vertex.frame.color="gold",   edge.curved=F,  layout=layout_with_mds(bacteria, dim = 2))

#@ Ploting one-mode networks
onemode <- bipartite.projection(bacteria)
#@ Chosing projection
Ind_nodes <- as.list(V(onemode$proj1))
#@ Adding attributes
V(onemode$proj2)$size <- c((datb.sp.table$`higher level`$weighted.betweenness+1)^5, 0.5)+5

#@ Ploting haplotypes
plot(onemode$proj1,edge.width=log(E(onemode$proj1)$weight), vertex.label=V(onemode$proj1)$name, layout=layout_with_kk(onemode$proj1) )
#@ Plotng counties
plot(onemode$proj2, edge.width=log(E(onemode$proj2)$weight), vertex.label=V(onemode$proj2)$name, layout=layout_with_kk(onemode$proj2))




