library(bipartite)
library(igraph)

file <- read.csv(file.choose(), row.names = 1, as.is=T)
# dat.mat = as.matrix(file[1:7,])
dat.mat=file
rownames(dat.mat) <- file$cc
dat.mat <- dat.mat[,2:ncol(dat.mat)]

plotweb(sortweb(t(dat.mat), sort.order="inc"), method="normal")   
visweb(sortweb(t(dat.mat),sort.order="inc"), type="diagonal", labsize=1,
       square="interaction", text="none", textsize = 1,circles=FALSE, frame=FALSE)# 
# bipartite 
dimdatb.net.table <- networklevel(t(dat.mat))
datb.sp.table  <- specieslevel(t(dat.mat))
datb.nest.table <- nestedcontribution(t(dat.mat))

bacteria <-  graph.incidence(t(dat.mat))
plot(bacteria)
#@ Checking type of network and names
V(bacteria)$type
V(bacteria)$name
#@ Generating attributes: colors and shapes
# cols= c("#DFBCAE", "#B4DDE5", "#E485C6", "#FCC63C",  "#F5FBC5", "#475BDE", "#B8B8B6", "#50E8CC", "#9BD198", "#9F5FD3", "#579BC8", "#E77375")
cols= c("#E77375", "#9F5FD3", "#475BDE", "#E485C6", "#F5FBC5", "#B8B8B6", "#FCC63C", "#50E8CC", "#B4DDE5", "#9BD198", "#DFBCAE", "#579BC8")
shapes = c(rep("circle", length(V(bacteria)$type[V(bacteria)$type == "FALSE"])), rep("square", length(V(bacteria)$type[V(bacteria)$type == "TRUE"])))
#@ Adding attributes
V(bacteria)$color <-  c(cols, rep("#F4A460", length(V(bacteria)$type[V(bacteria)$type == "TRUE"])))
V(bacteria)$size <- c((datb.sp.table$`lower level`$proportional.similarity)*20, (datb.sp.table$`higher level`$species.strength*5)) 
# V(bacteria)$size[V(bacteria)$size == Inf ] <- 0
#@ plot bipartite network
plot(bacteria, edge.arrow.size=1, vertex.shape=shapes, vertex.label.cex=1.2, vertex.label.color='black', vertex.frame.color="gray", 
     vertex.frame.color="gold",   edge.curved=F,  layout=layout_with_kk(bacteria, dim = 2, maxiter = 1000 , kkconst = 666 * vcount(bacteria) ))

dat <- datb.sp.table$`higher level`[1]
dat <- sort(dat, decreasing = T)
barplot(height = dat$degree, names = rownames(dat), ylim = c(0,12)) 

dat <- datb.sp.table$`lower level`[1]
dat <- datb.sp.table$`lower level`[18]
dat <- sort(dat, decreasing = T)
barplot(height = dat[,1], names = rownames(dat), ylim = c(0,10), las =2, cex.names = 2, col.lab = 2) 

#@ Ploting one-mode networks 
onemode <- bipartite.projection(bacteria)
#@ Chosing projection
Ind_nodes <- as.list(V(onemode$proj1))
#@ Adding attributes
V(onemode$proj2)$size <- c((datb.sp.table$`higher level`$weighted.betweenness+1)^5)+5

#@ Ploting haplotypes
plot(onemode$proj1,edge.width=E(onemode$proj1)$weight/2, vertex.label=V(onemode$proj1)$name, layout=layout_with_kk(onemode$proj1) )
#@ Plotng counties
plot(onemode$proj2, edge.width=E(onemode$proj2)$weight/2, vertex.label=V(onemode$proj2)$name, layout=layout_with_kk(onemode$proj2))



#geoBURST

adj.mat <- read.csv(file.choose(), row.names = 1, as.is=T)
adj.mat[is.na(adj.mat)] <- 0 

colnames(adj.mat) <- c(1:ncol(adj.mat))

geoBURST <- graph_from_adjacency_matrix(as.matrix(adj.mat),diag = F)
plot(geoBURST, edge.arrow.size=0, vertex.frame.color="gray", vertex.color="#C8C8FF", edge.width = 1.5, edge.color= "black")
