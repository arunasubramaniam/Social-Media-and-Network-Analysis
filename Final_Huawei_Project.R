#Background about data

#Huawei R & D nowadays are using the social network analysis tools and techniques 
#to enhance their business positions. One of the Major success behind Huawei Company 
#is that they promote their products through Social Media.
#Huawei Facebook Communication Network is Directed and Labeled having 1000 nodes and 100,306 edges
#Huawei Twitter Communication Network is Directed and Labeled having 1000 Nodes and 500,630 edges
#Huawei Instagram Communication Network is Directed and Labeled having 1000 Nodes and edges

getwd()
setwd("E:/UIC Fall 18/IDS 564 Social Media and Network Analysis/Project/Huawei Social Network Data")
library(igraph)
library(Matrix)


##############################################################################################
###################################### Facebook ###############################################
##############################################################################################


#Load the facebook data
data = read.csv("Facebook_Data.csv",header=TRUE,row.names=1,check.names=FALSE) # choose an adjacency matrix from a .csv file
m = as.matrix(data) # coerces the data set as a matrix
g = graph.adjacency(m,mode="directed",weighted=NULL) # this will create an 'igraph object'

ecount(g)   #[1] 100306
vcount(g)   #[1]1000
E(g)$color = "gray"
E(g)$width = .5
E(g)$arrow.width = .25
V(g)$label.color = "black"
V(g)$color = "dodgerblue"
V(g)$size = 4

set.seed(40)
l <- layout.fruchterman.reingold(g)

pdf("facebook_network.pdf")
plot(g,layout=l,rescale=TRUE,axes=FALSE,ylim=c(-1,1),asp=0,vertex.label=NA)
dev.off()

tab <- data.frame(sort(colSums(g[]), decreasing = T))
colnames(tab)<- "number of connections"
tab$names <- row.names(tab)
row.names(tab) <- NULL

min(colSums(g[])) 
# 64
max(colSums(g[]))
# 127
sort(unique(colSums(g[])),decreasing = T)

V(g)$shape <- "circle"
V(g)[c("Engkos Kosasih", "Ernie")]$shape <- "square"
V(g)[c("Aleisha")]$shape <- "circle"
g_vertexShape <- V(g)$shape

V(g)$color <- "grey"
V(g)[c("Engkos Kosasih", "Ernie")]$color <- "yellow"
V(g)[c("Aleisha")]$color <- "blue"
g_vertexColor<- V(g)$color

closeness(g ,mode = "in", v=V(g))
closeness(g,mode = "out", v=V(g))


layout_fr <- layout.fruchterman.reingold(g)
plot(g,
     vertex.color = g_vertexColor, # change color of nodes
     vertex.shape = g_vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .4, # change size of labels to 75% of original size
     vertex.size = closeness(g)*10000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr)

summary(g)
is.connected(g, mode="strong")
#[1] TRUE
is.simple(g)
#[1] TRUE

reciprocity(g)  #[1] 1

diameter(g)    #[1] 3

transitivity(g)  #[1] 0.1002897 )# global clustering: the ratio of the triangles
                                 # and the connected triples in the graph.
transitivity(g, type = "average")  #[1] 0.1003296     #average clustering
mean(transitivity(g, type = "local"), na.rm = TRUE)  #[1] 0.02495558  # local clustering

average.path.length(g, directed=TRUE, unconnected=TRUE)  #[1] 1.899632
edge_density(g)  #[1] 0.1004064

############## Cliques and Clusters ####################

table(sapply(maximal.cliques(g), length))
#  2     3     4     5     6 
#  5   61153 38225   855     1 

C <- get.adjacency(g, sparse=FALSE)

# Local clustering coefficients
clustering_contact <- transitivity(g, type="local", vids=V(g))

c1 = transitivity(g, "local")
d1 = degree(g)

plot(d1, c1, xlab = "Degree", ylab = "Clustering Coefficient",
     main = "Clustering Coefficient Vs Degree", col = "orange", pch = 16,
     cex = 0.7)

assortativity_degree(g, directed = TRUE)  #[1] -0.007703903
max(degree(g))  #[1] 254

#################### centrality measures#######################
# Degree centrality
degree_contact <- degree(g)
# Node betweenness
betweens_contact <- round(betweenness(g, v=V(g), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_contact<-edge.betweenness(g, e=E(g), directed = TRUE)
# Local clustering coefficients
clustering_contact <- transitivity(g, type="local", vids=V(g))

###########################Community structures#########################################


# Fast Greedy is only for undirected networks
# comm_fs <- fastgreedy.community(g)
# c.m.fs <- membership(comm_fs)
# plot(comm_fs,g, vertex.label= NA, vertex.size=2, main = "Fast Greedy Algorithm")

# comm_eb <- cluster_edge_betweenness(g, weights = NULL, directed = TRUE,
#                          edge.betweenness = TRUE, merges = F, bridges = F,
#                          modularity = TRUE, membership = F)
# c.m.eb <- membership(comm_eb)
# modularity(c.m.lp)

# plot(comm_eb,g, vertex.label= NA, vertex.size=2, main = "Girvan Newman Algorithm")
# girvNew <- cluster_edge_betweenness(g, modularity = F)
# girvNew_sizesComm <- sizes(girvNew) 
# girvNew_numComm <- length(girvNew_sizesComm)
# girvNew_modularity <- modularity(girvNew)

# # Label Propagation is only for undirected networks
# comm_lp <- label.propagation.community(g)
# c.m.lp <- membership(comm_lp)
# plot(comm_lp,g, vertex.label= NA, vertex.size=2, main = "Label Propagation Algorithm")

# Walktrap Algorithm
comm_wt <- walktrap.community(g)
c.m.wt <- membership(comm_wt)
table(c.m.wt)
# 1   2   3   4   5   6 
# 225 281 153 147 149  45
plot(comm_wt,g, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.05, main = "Walktrap Algorithm")
modularity(comm_wt)

# # Spinglass 
# comm_sg <- spinglass.community(g)
# c.m.sg <- membership(comm_sg)
# plot(comm_sg,g, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.05, main = "Spin Glass Algorithm")

#InfoMap 
comm_infomap_fb <- cluster_infomap(g, e.weights = NULL, v.weights = NULL, nb.trials = 10,
                modularity = T)
plot(comm_infomap_fb,g, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.005, main = "InfoMap Algorithm")
modularity(comm_infomap_fb)

############################ Sub-graph for communities ##################################
fb_comm1 <- V(g)[c.m.wt==1]
fb_comm1 <- V(g)[c.m.wt==2]

sub_net_cm1 <- induced.subgraph(g, v=fb_comm1)

layout_fr <- layout.fruchterman.reingold(sub_net_cm1)

V(sub_net_cm1)$shape <- "circle" 
vertexShape <- V(sub_net_cm1)$shape

V(sub_net_cm1)$color <- "grey"
vertexColor<- V(sub_net_cm1)$color


plot(sub_net_cm1,
     vertex.color = vertexColor, # change color of nodes
     vertex.shape = vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .5, # change size of labels to 75% of original size
     vertex.size = closeness(sub_net_cm1)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr)

############################ Sub-graph: isloation top 10 nodes based on the number of their connections #################################3

#Isolate 10% of the nodes and induce a subgraph

top_10 <- tab[tab$`number of connections`>=112,]
sub_net <- induced.subgraph(g, v=top_10$names)
vcount(sub_net)
# 107
ecount(sub_net)
# 1544

diameter(sub_net)
#  3
#closeness
closeness(g ,mode = "in", v=top_10$names)
closeness(g,mode = "out", v=top_10$names)

betweenness(sub_net)
eigen_centrality(sub_net)
E(sub_net)
V(sub_net)

components(sub_net, mode = c("weak", "strong"))
# 1 component with 107 nodes

least_conn_fb <- tab[tab$`number of connections`==112,]

V(sub_net)$shape <- "circle" 
V(sub_net)[c("Engkos Kosasih", "Ernie")]$shape <- "square"
V(sub_net)[least_conn_fb$names]$shape <- "circle"
vertexShape <- V(sub_net)$shape

V(sub_net)$color <- "grey"
V(sub_net)[c("Engkos Kosasih", "Ernie")]$color <- "yellow"
V(sub_net)[least_conn_fb$names]$color <- "blue"
vertexColor<- V(sub_net)$color

layout_kk <- layout.fruchterman.reingold(sub_net)
plot(sub_net,
     vertex.color = vertexColor, # change color of nodes
     vertex.shape = vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .90, # change size of labels to 75% of original size
     vertex.size = closeness(sub_net)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_kk)

#Induced subgraph based on communities
is_fb_wt <- induced.subgraph(g, which(c.m.wt == 1))
layout_kk <- layout.fruchterman.reingold(is_fb_wt)
plot(is_fb_wt,
     # vertex.color = vertexColor, # change color of nodes
     # vertex.shape = vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .90, # change size of labels to 75% of original size
     vertex.size = closeness(is_fb_wt)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_kk)

############################## Epidemic - Info spread ######################################

simulate_sir = function(g, simlength=15, p.t=0.2,
                        display_net=TRUE, removeafter=2,
                        susceptibleafter=1000, startNode = 0)
{
  links = get.edgelist(g)
  N = vcount(g)
  time_stats = list()
  
  # Number of nodes in S, I, or R status in each round of time
  time_stats$infected_t = rep(1, simlength)
  time_stats$removed_t = rep(0, simlength)
  time_stats$susceptible_t = rep(N-1, simlength)
  time_stats$unInfected = rep(TRUE, N)
  
  infected = logical(N)
  susceptible = rep(TRUE, N)
  removed = logical(N)
  if (startNode == 0)
  {
    patientzero = sample(N,1)
  }
  else
  {
    patientzero = startNode
  }
  
  
  # Initialize a vector that keeps track of the time
  infected_time = rep(0, N)
  removed_time = rep(0, N)
  
  # Patient zero  
  infected[patientzero] = TRUE
  susceptible[patientzero] = FALSE
  time_stats$unInfected[patientzero] = FALSE
  infected_time[patientzero] = 1
  
  if (N > 50) {
    if (N < 100) {
      V(g)$size = 8
    }
    else {
      V(g)$size = 4
      V(g)$label = ""
    }
  }
  
  for (i in 1:simlength) {
    
    # Find the indeces of links that connect an
    # infected individual to an uninfected
    discordant.links = which(xor(infected[as.integer(links[,1])],
                                 infected[as.integer(links[,2])]))
    # Determine randomly which of the discordant links
    # transmit the disease
    transmit = rbinom(length(discordant.links), 1, p.t)
    transmitter.links = discordant.links[transmit==1]
    
    nodes.of.transmitter.links =
      unique(as.vector(as.integer(links[transmitter.links,1:2])))
    
    # Remove any immune nodes from the nodes to transmit
    any.immune.nodes = which(removed[nodes.of.transmitter.links]>0)
    if (length(any.immune.nodes) > 0){
      nodes.of.transmitter.links =
        nodes.of.transmitter.links[-any.immune.nodes]
    }
    
    if (display_net)
    {
      fixlayout = layout.kamada.kawai(g)
      node.colour = rep("SkyBlue2", N)
      node.colour[infected] = "red"
      node.colour[removed] = "yellow"
      node.colour[susceptible] = "SkyBlue2"
      
      bmp(paste0(i, "_Epidemics_Graph.jpg"), 1800, 1800, pointsize = 30)
      plot(g, layout = fixlayout,
           main=paste("Time =", i, " out of ", simlength),
           vertex.color=node.colour)
      dev.off()
    }
    
    infected[nodes.of.transmitter.links] = TRUE
    susceptible[nodes.of.transmitter.links] = FALSE
    time_stats$unInfected[nodes.of.transmitter.links] = FALSE
    
    old_infected = which(infected > 0)
    new_infected = setdiff(nodes.of.transmitter.links, old_infected)
    infected_time[old_infected] = infected_time[old_infected] + 1
    infected_time[new_infected] = 1
    
    remove.them = infected_time > removeafter
    removed[remove.them] = TRUE
    infected[remove.them] = FALSE
    infected_time[remove.them] = 0
    
    index_r = removed_time > 0
    removed_time[index_r] = removed_time[index_r] + 1
    removed_time[remove.them] = 1
    
    susceptible.again = removed_time > susceptibleafter
    susceptible[susceptible.again] = TRUE
    removed[susceptible.again] = FALSE
    removed_time[susceptible.again] = 0
    
    time_stats$infected_t[i] = sum(infected)
    time_stats$removed_t[i] = sum(removed)
    time_stats$susceptible_t[i] = sum(susceptible)
  }
  return(time_stats)
}

plotGraphs = function(infected_time, N = 0){
  par(mfrow=c(4,1))
  plot(infected_time$infected_t, type="l", col="red",
       ylab="Infected Over Time", xlab = "Time Index")
  plot(infected_time$removed_t, type="l", col="blue",
       ylab="Removed Over Time", xlab = "Time Index")
  plot(infected_time$susceptible_t, type="l", col="yellow",
       ylab="Susceptible Over Time", xlab = "Time Index")
  unInfected_colors = rep("white", N)
  unInfected_colors[infected_time$unInfected] = "Darkgreen"
  plot(as.integer(infected_time$unInfected), type = "p",
       col=unInfected_colors, pch=19,cex = 1,
       ylim = range(0:2), xlab = "Nodes", yaxt = 'n',
       ylab = "", main = "Uninfected Nodes")
}


V(g)$indeg <- degree(g, mode = "in")
V(g)[V(g)$indeg == max(V(g)$indeg)]  
V(g)[V(g)$indeg == min(V(g)$indeg)]  

# Consider the most & least influential nodes
most = which(V(g)$name == "Engkos Kosasih")
least = which(V(g)$name == "Aleisha")

V(g)$name = c(1:length(V(g)))
V(g)$label = c(1:length(V(g)))


FacebookModel = simulate_sir(g, simlength = 50, p.t=0.01,
                            removeafter = 2, susceptibleafter = 7,
                            display_net = F, startNode = most)
plotGraphs(FacebookModel, vcount(g))


FacebookModel1= simulate_sir(g, simlength = 50, p.t=0.01,
                           removeafter = 2, susceptibleafter = 7,
                           display_net = F, startNode = least)
plotGraphs(FacebookModel1, vcount(g))

#changed r = 5 and s=10
FacebookModel2 = simulate_sir(g, simlength = 50, p.t=0.01,
                             removeafter = 5, susceptibleafter = 10,
                             display_net = F, startNode = most)
plotGraphs(FacebookModel2, vcount(g))


FacebookModel3= simulate_sir(g, simlength = 50, p.t=0.01,
                             removeafter = 5, susceptibleafter = 10,
                             display_net = F, startNode = least)
plotGraphs(FacebookModel3, vcount(g))


##############################################################################################
###################################### Twitter ###############################################
##############################################################################################


#Loading Twitter data
data1=read.csv("Twitter_Data.csv",header=TRUE,row.names=1,check.names=FALSE) # choose an adjacency matrix from a .csv file
m1=as.matrix(data1) # coerces the data set as a matrix
g1=graph.adjacency(m1,mode="directed",weighted=NULL) # this will create an 'igraph object'

ecount(g1)  #[1] 500630
vcount(g1)  #[1] 1000
E(g1)$color = "gray"
E(g1)$width = .5
E(g1)$arrow.width = .25
V(g1)$label.color = "black"
V(g1)$color = "dodgerblue"
V(g1)$size = 1

set.seed(40)
l <- layout.fruchterman.reingold(g1)

pdf("twitter_network.pdf")
plot(g1,layout=l,rescale=TRUE,axes=FALSE,ylim=c(-1,1),asp=0,vertex.label=NA)
dev.off()

tab1 <- data.frame(sort(colSums(g1[]), decreasing = T))
colnames(tab1)<- "number of connections"
tab1$names <- row.names(tab1)
row.names(tab1) <- NULL

min(colSums(g1[])) 
# 445
max(colSums(g1[]))
# 545
sort(unique(colSums(g1[])),decreasing = T)

V(g1)$shape <- "circle"
V(g1)[c("Dililah")]$shape <- "square"
V(g1)[c("Saleem Hassan")]$shape <- "circle"
g1_vertexShape <- V(g1)$shape

V(g1)$color <- "grey"
V(g1)[c("Dililah")]$color <- "yellow"
V(g1)[c("Saleem Hassan")]$color <- "blue"
g1_vertexColor<- V(g1)$color

closeness(g1 ,mode = "in", v=V(g1))
closeness(g1,mode = "out", v=V(g1))

layout_fr1 <- layout.fruchterman.reingold(g1)
plot(g1,
     vertex.color = g1_vertexColor, # change color of nodes
     vertex.shape = g1_vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .4, # change size of labels to 75% of original size
     vertex.size = closeness(g1)*10000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr1)

is.connected(g1, mode="strong")
#[1] TRUE
is.connected(g1, mode="weak")
#[1] TRUE


is.simple(g1)  #[1] TRUE
reciprocity(g1)  #[1] 1
diameter(g1)     #[1] 2
transitivity(g1)  #[1]  0.5012008
transitivity(g1, type = "average")  #[1] 0.50121     #average clustering
mean(transitivity(g1, type = "local"), na.rm = TRUE)  #[1] 0.1251771  # local clustering

average.path.length(g1, directed=TRUE, unconnected=TRUE)  #[1] 1.498869
edge_density(g1)  #[1]  0.5011311
############## Cliques and Clusters ####################
table(sapply(maximal.cliques(g1), length))


B <- get.adjacency(g1, sparse=FALSE)

# Local clustering coefficients
clustering_contact1 <- transitivity(g1, type="local", vids=V(g1))

c2 = transitivity(g1, "local")
d2 = degree(g1)

plot(d2, c2, xlab = "Degree", ylab = "Clustering Coefficient",
     main = "Clustering Coefficient Vs Degree", col = "orange", pch = 16,
     cex = 0.7)


assortativity_degree(g1, directed = TRUE)  #[1] -0.004438593
max(degree(g1))  #[1] 1090

#################### centrality measures#######################
# Degree centrality
degree_contact <- degree(g1)
# Node betweenness
betweens_contact <- round(betweenness(g1, v=V(g1), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_contact<-edge.betweenness(g1, e=E(g1), directed = TRUE)
# Local clustering coefficients
clustering_contact <- transitivity(g1, type="local", vids=V(g1))

###########################Community structures#########################################


# Fast Greedy is only for undirected networks
# comm_fs <- fastgreedy.community(g1)
# c.m.fs <- membership(comm_fs)
# plot(comm_fs,g, vertex.label= NA, vertex.size=2, main = "Fast Greedy Algorithm")

# # Girvan Newman is only for undirected networks
# comm_eb1 <- edge.betweenness.community(g1)
# c.m.eb1 <- membership(comm_eb1)
# plot(comm_eb1,g1, vertex.label= NA, vertex.size=2, main = "Girvan Newman Algorithm")
# # girvNew <- cluster_edge_betweenness(g, modularity = TRUE)
# # girvNew_sizesComm <- sizes(girvNew) 
# # girvNew_numComm <- length(girvNew_sizesComm)
# # girvNew_modularity <- modularity(girvNew)
# 
# # Label Propagation is only for undirected networks
# comm_lp1 <- label.propagation.community(g1)
# c.m.lp1 <- membership(comm_lp1)
# plot(comm_lp1,g1, vertex.label= NA, vertex.size=2, main = "Label Propagation Algorithm")

# Walktrap Algorithm
comm_wt1 <- walktrap.community(g1)  
c.m.wt1 <- membership(comm_wt1)
table(c.m.wt1)
# 1  2 
# 74 33
plot(comm_wt1,g1, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.05, main = "Walktrap Algorithm")

# # Spinglass 
# comm_sg1 <- spinglass.community(g1)
# c.m.sg1<- membership(comm_sg1)
# plot(comm_sg1,g1, vertex.label= NA, vertex.size=2, main = "Spin Glass Algorithm")

#InfoMap 
comm_infomap_twi <- cluster_infomap(g1, e.weights = NULL, v.weights = NULL, nb.trials = 10,
                                    modularity = F)
plot(comm_infomap_twi,g1, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.005, main = "InfoMap Algorithm")

############################ Sub-graph based on community ##################################
tw_comm1 <- V(g1)[c.m.wt==1]
tw_comm1 <- V(g1)[c.m.wt==2]

sub_net1_cm1 <- induced.subgraph(g1, v=tw_comm1)

layout_fr1 <- layout.fruchterman.reingold(sub_net1_cm1)

V(sub_net1_cm1)$shape <- "circle" 
vertexShape <- V(sub_net1_cm1)$shape

V(sub_net1_cm1)$color <- "grey"
vertexColor<- V(sub_net1_cm1)$color


plot(sub_net1_cm1,
     vertex.color = vertexColor, # change color of nodes
     vertex.shape = vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .5, # change size of labels to 75% of original size
     vertex.size = closeness(sub_net1_cm1)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr1)

############################ Sub-graph: isloation top 10 nodes based on the number of their connections #################################3

#Isolate 10% of the nodes and induce a subgraph

top_10_twitter <- tab1[tab1$`number of connections`>=522,]
sub_net1 <- induced.subgraph(g1, v=top_10_twitter$names)

vcount(sub_net1)
# [1] 107
ecount(sub_net1)
#[1] 6336

diameter(sub_net1)
# 2
#closeness
closeness(g1 ,mode = "in", v=top_10_twitter$names)
closeness(g1,mode = "out", v=top_10_twitter$names)

betweenness(sub_net1)
eigen_centrality(sub_net1)

E(sub_net1)

components(sub_net1, mode = c("weak", "strong"))
# 1 component with 77 nodes
least_conn_tw <- tab1[tab1$`number of connections`==522,]

V(sub_net1)$shape <- "circle"
V(sub_net1)[c("Dililah")]$shape <- "square"
V(sub_net1)[least_conn_tw$names]$shape <- "circle"
tvertexShape <- V(sub_net1)$shape

V(sub_net1)$color <- "grey"
V(sub_net1)[c("Dililah")]$color <- "yellow"
V(sub_net1)[least_conn_tw$names]$color <- "blue"
tvertexColor<- V(sub_net1)$color

closeness(sub_net1 ,mode = "in", v=V(sub_net1))
closeness(sub_net1,mode = "out", v=V(sub_net1))

layout_fr1 <- layout.fruchterman.reingold(sub_net1)
plot(sub_net1,
     vertex.color = tvertexColor, # change color of nodes
     vertex.shape = tvertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = 0.85, # change size of labels to 75% of original size
     vertex.size = closeness(sub_net1)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr1)

############################## Epidemic - Info spread ######################################

simulate_sir = function(g1, simlength=15, p.t=0.2,
                        display_net=TRUE, removeafter=2,
                        susceptibleafter=1000, startNode = 0)
{
  links = get.edgelist(g1)
  N = vcount(g1)
  time_stats = list()
  
  # Number of nodes in S, I, or R status in each round of time
  time_stats$infected_t = rep(1, simlength)
  time_stats$removed_t = rep(0, simlength)
  time_stats$susceptible_t = rep(N-1, simlength)
  time_stats$unInfected = rep(TRUE, N)
  
  infected = logical(N)
  susceptible = rep(TRUE, N)
  removed = logical(N)
  if (startNode == 0)
  {
    patientzero = sample(N,1)
  }
  else
  {
    patientzero = startNode
  }
  
  
  # Initialize a vector that keeps track of the time
  infected_time = rep(0, N)
  removed_time = rep(0, N)
  
  # Patient zero  
  infected[patientzero] = TRUE
  susceptible[patientzero] = FALSE
  time_stats$unInfected[patientzero] = FALSE
  infected_time[patientzero] = 1
  
  if (N > 50) {
    if (N < 100) {
      V(g1)$size = 8
    }
    else {
      V(g1)$size = 4
      V(g1)$label = ""
    }
  }
  
  for (i in 1:simlength) {
    
    # Find the indeces of links that connect an
    # infected individual to an uninfected
    discordant.links = which(xor(infected[as.integer(links[,1])],
                                 infected[as.integer(links[,2])]))
    # Determine randomly which of the discordant links
    # transmit the disease
    transmit = rbinom(length(discordant.links), 1, p.t)
    transmitter.links = discordant.links[transmit==1]
    
    nodes.of.transmitter.links =
      unique(as.vector(as.integer(links[transmitter.links,1:2])))
    
    # Remove any immune nodes from the nodes to transmit
    any.immune.nodes = which(removed[nodes.of.transmitter.links]>0)
    if (length(any.immune.nodes) > 0){
      nodes.of.transmitter.links =
        nodes.of.transmitter.links[-any.immune.nodes]
    }
    
    if (display_net)
    {
      fixlayout = layout.kamada.kawai(g1)
      node.colour = rep("SkyBlue2", N)
      node.colour[infected] = "red"
      node.colour[removed] = "yellow"
      node.colour[susceptible] = "SkyBlue2"
      
      bmp(paste0(i, "_Epidemics_Graph.jpg"), 1800, 1800, pointsize = 30)
      plot(g1, layout = fixlayout,
           main=paste("Time =", i, " out of ", simlength),
           vertex.color=node.colour)
      dev.off()
    }
    
    infected[nodes.of.transmitter.links] = TRUE
    susceptible[nodes.of.transmitter.links] = FALSE
    time_stats$unInfected[nodes.of.transmitter.links] = FALSE
    
    old_infected = which(infected > 0)
    new_infected = setdiff(nodes.of.transmitter.links, old_infected)
    infected_time[old_infected] = infected_time[old_infected] + 1
    infected_time[new_infected] = 1
    
    remove.them = infected_time > removeafter
    removed[remove.them] = TRUE
    infected[remove.them] = FALSE
    infected_time[remove.them] = 0
    
    index_r = removed_time > 0
    removed_time[index_r] = removed_time[index_r] + 1
    removed_time[remove.them] = 1
    
    susceptible.again = removed_time > susceptibleafter
    susceptible[susceptible.again] = TRUE
    removed[susceptible.again] = FALSE
    removed_time[susceptible.again] = 0
    
    time_stats$infected_t[i] = sum(infected)
    time_stats$removed_t[i] = sum(removed)
    time_stats$susceptible_t[i] = sum(susceptible)
  }
  return(time_stats)
}

plotGraphs = function(infected_time, N = 0){
  par(mfrow=c(4,1))
  plot(infected_time$infected_t, type="l", col="red",
       ylab="Infected Over Time", xlab = "Time Index")
  plot(infected_time$removed_t, type="l", col="blue",
       ylab="Removed Over Time", xlab = "Time Index")
  plot(infected_time$susceptible_t, type="l", col="yellow",
       ylab="Susceptible Over Time", xlab = "Time Index")
  unInfected_colors = rep("white", N)
  unInfected_colors[infected_time$unInfected] = "Darkgreen"
  plot(as.integer(infected_time$unInfected), type = "p",
       col=unInfected_colors, pch=19,cex = 1,
       ylim = range(0:2), xlab = "Nodes", yaxt = 'n',
       ylab = "", main = "Uninfected Nodes")
}


V(g1)$indeg <- degree(g1, mode = "in")
V(g1)[V(g1)$indeg == max(V(g1)$indeg)]  
V(g1)[V(g1)$indeg == min(V(g1)$indeg)]  

# Consider the most & least influential nodes
most = which(V(g1)$name == "410")
least = which(V(g1)$name == "561")

V(g1)$name = c(1:length(V(g1)))
V(g1)$label = c(1:length(V(g1)))


#change the r=5 and s=10 for both most and least
TwitterModel = simulate_sir(g1, simlength = 50, p.t=0.01,
                             removeafter = 2, susceptibleafter = 7,
                             display_net = F, startNode = most)
plotGraphs(TwitterModel, vcount(g1))

TwitterModel1= simulate_sir(g1, simlength = 50, p.t=0.01,
                             removeafter = 2, susceptibleafter = 7,
                             display_net = F, startNode = least)
plotGraphs(TwitterModel1, vcount(g1))

#change the r=5 and s=10 for both most and least
TwitterModel2 = simulate_sir(g1, simlength = 50, p.t=0.01,
                            removeafter = 5, susceptibleafter = 10,
                            display_net = F, startNode = most)
plotGraphs(TwitterModel2, vcount(g1))

TwitterModel3= simulate_sir(g1, simlength = 50, p.t=0.01,
                            removeafter = 5, susceptibleafter = 10,
                            display_net = F, startNode = least)
plotGraphs(TwitterModel3, vcount(g1))



##############################################################################################
###################################### Instagram ###############################################
##############################################################################################


#Load the Instagram data
data2 = read.csv("Instagram_Data.csv",header=TRUE,row.names=1,check.names=FALSE) # choose an adjacency matrix from a .csv file
m2 = as.matrix(data2) # coerces the data set as a matrix
g2 = graph.adjacency(m2,mode="directed",weighted=NULL) # this will create an 'igraph object'

ecount(g2)   #[1] 9866
vcount(g2)   #[1]1000
E(g2)$color = "gray"
E(g2)$width = .3
E(g2)$arrow.width = .15
V(g2)$label.color = "black"
V(g2)$color = "dodgerblue"
V(g2)$size = 1

set.seed(40)
l <- layout.fruchterman.reingold(g2)

pdf("Instagram_network.pdf")
plot(g2,layout=l,rescale=TRUE,axes=FALSE,ylim=c(-1,1),asp=0,vertex.label=NA)
dev.off()

tab2 <- data.frame(sort(colSums(g2[]), decreasing = T))
colnames(tab2)<- "number of connections"
tab2$names <- row.names(tab2)
row.names(tab2) <- NULL

min(colSums(g2[])) 
# 2
max(colSums(g2[]))
# 20
sort(unique(colSums(g2[])),decreasing = T)

V(g2)$shape <- "circle"
V(g2)[c("Alexis", "Ishku Ishku", "Alveena")]$shape <- "square"
V(g2)[c("Denno", "Rohan", "Fahad Rehman")]$shape <- "circle"
g2_vertexShape <- V(g2)$shape

V(g2)$color <- "grey"
V(g2)[c("Alexis", "Ishku Ishku", "Alveena")]$color <- "yellow"
V(g2)[c("Denno", "Rohan", "Fahad Rehman")]$color <- "blue"
g2_vertexColor<- V(g2)$color

closeness(g2 ,mode = "in", v=V(g2))
closeness(g2,mode = "out", v=V(g2))

layout_fr2 <- layout.fruchterman.reingold(g2)
plot(g2,
     vertex.color = g2_vertexColor, # change color of nodes
     vertex.shape = g2_vertexShape, #vertex
     vertex.label = NA,
     # vertex.label.color = "black", # change color of labels
     # vertex.label.cex = .5, # change size of labels to 75% of original size
     vertex.size = 5,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr2)

is.connected(g2, mode="strong")
#[1] TRUE
is.simple(g2)
#[1] TRUE

reciprocity(g2)  #[1] 1

diameter(g2)    #[1] 5

transitivity(g2)  #[1]  0.008399037
transitivity(g2, type = "average")  #[1] 0.009066097     #average clustering
mean(transitivity(g2, type = "local"), na.rm = TRUE)  #[1] 0.002121212  # local clustering

average.path.length(g2, directed=TRUE, unconnected=TRUE)  #[1]  3.273137
edge_density(g2)  #[1] 0.009875876

############## Cliques and Clusters ####################

table(sapply(maximal.cliques(g2), length))
#2      3 
#4537  136 

C2 <- get.adjacency(g2, sparse=FALSE)

# Local clustering coefficients
clustering_contact <- transitivity(g2, type="local", vids=V(g2))
c3 = transitivity(g2, "local")
d3 = degree(g2)

plot(d3, c3, xlab = "Degree", ylab = "Clustering Coefficient",
     main = "Clustering Coefficient Vs Degree", col = "orange", pch = 16,
     cex = 0.7)


assortativity_degree(g2, directed = TRUE)  #[1] 0.005353542
max(degree(g2))  #[1] 40


#################### centrality measures#######################
# Degree centrality
degree_contact <- degree(g2)
# Node betweenness
betweens_contact <- round(betweenness(g2, v=V(g2), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_contact<-edge.betweenness(g2, e=E(g2), directed = TRUE)

###########################Community structures#########################################


# Fast Greedy is only for undirected networks
# comm_fs <- fastgreedy.community(g2)
# c.m.fs <- membership(comm_fs)
# plot(comm_fs,g, vertex.label= NA, vertex.size=2, main = "Fast Greedy Algorithm")

# # Girvan Newman is only for undirected networks
# comm_eb2 <- edge.betweenness.community(g2)
# c.m.eb2 <- membership(comm_eb2)
# plot(comm_eb2,g2, vertex.label= NA, vertex.size=2, main = "Girvan Newman Algorithm")
# # girvNew <- cluster_edge_betweenness(g, modularity = TRUE)
# # girvNew_sizesComm <- sizes(girvNew) 
# # girvNew_numComm <- length(girvNew_sizesComm)
# # girvNew_modularity <- modularity(girvNew)
# 
# # Label Propagation is only for undirected networks
# comm_lp2 <- label.propagation.community(g2)
# c.m.lp2 <- membership(comm_lp2)
# plot(comm_lp2,g2, vertex.label= NA, vertex.size=2, main = "Label Propagation Algorithm")

# Walktrap Algorithm
comm_wt2 <- walktrap.community(g2)  
c.m.wt2 <- membership(comm_wt2)
table(c.m.wt2)
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
# 9  5  8 12  4  9 12  5 10 11  4 12  4  6  5  7  3  2  2  2  2  1  1  1  1
plot(comm_wt2,g2, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.005,main = "Walktrap Algorithm")

# Spinglass does not work with unconnect graph. isolate the big component and re-do 
# comm_sg2 <- spinglass.community(g2)
# c.m.sg2 <- membership(comm_sg2)
# plot(comm_sg2,g2, vertex.label= NA, vertex.size=2, main = "Spin Glass Algorithm")

#InfoMap 
comm_infomap_ins <- cluster_infomap(g2, e.weights = NULL, v.weights = NULL, nb.trials = 10,
                                    modularity = F)
plot(comm_infomap_ins,g2, vertex.label= NA, vertex.size=2,edge.arrow.size = 0.005, main = "InfoMap Algorithm")



############################ Sub-graph for communities ##################################
in_comm1 <- V(g)[c.m.wt2==1]
in_comm2 <- V(g)[c.m.wt2==2]

sub_net2_cm1 <- induced.subgraph(g2, v=in_comm1)

layout_fr2 <- layout.fruchterman.reingold(sub_net2_cm1)

V(sub_net2_cm1)$shape <- "circle" 
vertexShape <- V(sub_net2_cm1)$shape

V(sub_net2_cm1)$color <- "grey"
vertexColor<- V(sub_net2_cm1)$color


plot(sub_net2_cm1,
     vertex.color = vertexColor, # change color of nodes
     vertex.shape = vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .5, # change size of labels to 75% of original size
     vertex.size = closeness(sub_net2_cm1)*1000,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_fr2)



############################ Sub-graph: isloation top 10 nodes based on the number of their connections #################################3

#Isolate 10% of the nodes and induce a subgraph

top_10_insta <- tab2[tab2$`number of connections`>=14,]
sub_net2 <- induced.subgraph(g2, v=top_10_insta$names)

vcount(sub_net2)
# [1] 138
ecount(sub_net2)
# [1] 440

diameter(sub_net2)
# 9
#closeness
closeness(sub_net2 ,mode = "in", v=top_10_insta$names)
closeness(sub_net2,mode = "out", v=top_10_insta$names)

betweenness(sub_net2)
eigen_centrality(sub_net2)
E(sub_net2)

components(sub_net2, mode = c("weak", "strong"))
# 11 unconnected components with 68,2,1,1,1,1,2,1,1,1,1 nodes

least_conn_in <- tab2[tab2$`number of connections`==14,]

V(sub_net2)$shape <- "circle"
V(sub_net2)[c("Alexis", "Ishku Ishku", "Alveena")]$shape <- "square"
V(sub_net2)[least_conn_in$name]$shape <- "circle"
sb2_vertexShape <- V(sub_net2)$shape

V(sub_net2)$color <- "grey"
V(sub_net2)[c("Alexis", "Ishku Ishku", "Alveena")]$color <- "yellow"
V(sub_net2)[least_conn_in$name]$color <- "blue"
sb2_vertexColor<- V(sub_net2)$color


layout_frsb2 <- layout.fruchterman.reingold(sub_net2)
plot(sub_net2,
     vertex.color = sb2_vertexColor, # change color of nodes
     vertex.shape = sb2_vertexShape, #vertex
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .5, # change size of labels to 75% of original size
     vertex.size = 5,#size of the nodes proportional to betweeness
     edge.color="grey", # change edge color to grey
     edge.arrow.size = 0.005,
     layout = layout_frsb2)


############################## Epidemic - Info spread ######################################

simulate_sir = function(g2, simlength=15, p.t=0.2,
                        display_net=TRUE, removeafter=2,
                        susceptibleafter=1000, startNode = 0)
{
  links = get.edgelist(g2)
  N = vcount(g2)
  time_stats = list()
  
  # Number of nodes in S, I, or R status in each round of time
  time_stats$infected_t = rep(1, simlength)
  time_stats$removed_t = rep(0, simlength)
  time_stats$susceptible_t = rep(N-1, simlength)
  time_stats$unInfected = rep(TRUE, N)
  
  infected = logical(N)
  susceptible = rep(TRUE, N)
  removed = logical(N)
  if (startNode == 0)
  {
    patientzero = sample(N,1)
  }
  else
  {
    patientzero = startNode
  }
  
  
  # Initialize a vector that keeps track of the time
  infected_time = rep(0, N)
  removed_time = rep(0, N)
  
  # Patient zero  
  infected[patientzero] = TRUE
  susceptible[patientzero] = FALSE
  time_stats$unInfected[patientzero] = FALSE
  infected_time[patientzero] = 1
  
  if (N > 50) {
    if (N < 100) {
      V(g2)$size = 8
    }
    else {
      V(g2)$size = 4
      V(g2)$label = ""
    }
  }
  
  for (i in 1:simlength) {
    
    # Find the indeces of links that connect an
    # infected individual to an uninfected
    discordant.links = which(xor(infected[as.integer(links[,1])],
                                 infected[as.integer(links[,2])]))
    # Determine randomly which of the discordant links
    # transmit the disease
    transmit = rbinom(length(discordant.links), 1, p.t)
    transmitter.links = discordant.links[transmit==1]
    
    nodes.of.transmitter.links =
      unique(as.vector(as.integer(links[transmitter.links,1:2])))
    
    # Remove any immune nodes from the nodes to transmit
    any.immune.nodes = which(removed[nodes.of.transmitter.links]>0)
    if (length(any.immune.nodes) > 0){
      nodes.of.transmitter.links =
        nodes.of.transmitter.links[-any.immune.nodes]
    }
    
    if (display_net)
    {
      fixlayout = layout.kamada.kawai(g)
      node.colour = rep("SkyBlue2", N)
      node.colour[infected] = "red"
      node.colour[removed] = "yellow"
      node.colour[susceptible] = "SkyBlue2"
      
      bmp(paste0(i, "_Epidemics_Graph.jpg"), 1800, 1800, pointsize = 30)
      plot(g, layout = fixlayout,
           main=paste("Time =", i, " out of ", simlength),
           vertex.color=node.colour)
      dev.off()
    }
    
    infected[nodes.of.transmitter.links] = TRUE
    susceptible[nodes.of.transmitter.links] = FALSE
    time_stats$unInfected[nodes.of.transmitter.links] = FALSE
    
    old_infected = which(infected > 0)
    new_infected = setdiff(nodes.of.transmitter.links, old_infected)
    infected_time[old_infected] = infected_time[old_infected] + 1
    infected_time[new_infected] = 1
    
    remove.them = infected_time > removeafter
    removed[remove.them] = TRUE
    infected[remove.them] = FALSE
    infected_time[remove.them] = 0
    
    index_r = removed_time > 0
    removed_time[index_r] = removed_time[index_r] + 1
    removed_time[remove.them] = 1
    
    susceptible.again = removed_time > susceptibleafter
    susceptible[susceptible.again] = TRUE
    removed[susceptible.again] = FALSE
    removed_time[susceptible.again] = 0
    
    time_stats$infected_t[i] = sum(infected)
    time_stats$removed_t[i] = sum(removed)
    time_stats$susceptible_t[i] = sum(susceptible)
  }
  return(time_stats)
}

plotGraphs = function(infected_time, N = 0){
  par(mfrow=c(4,1))
  plot(infected_time$infected_t, type="l", col="red",
       ylab="Infected Over Time", xlab = "Time Index")
  plot(infected_time$removed_t, type="l", col="blue",
       ylab="Removed Over Time", xlab = "Time Index")
  plot(infected_time$susceptible_t, type="l", col="yellow",
       ylab="Susceptible Over Time", xlab = "Time Index")
  unInfected_colors = rep("white", N)
  unInfected_colors[infected_time$unInfected] = "Darkgreen"
  plot(as.integer(infected_time$unInfected), type = "p",
       col=unInfected_colors, pch=19,cex = 1,
       ylim = range(0:2), xlab = "Nodes", yaxt = 'n',
       ylab = "", main = "Uninfected Nodes")
}


V(g2)$indeg <- degree(g2, mode = "in")
V(g2)[V(g2)$indeg == max(V(g2)$indeg)]  
V(g2)[V(g2)$indeg == min(V(g2)$indeg)]  

# Consider the most & least influential nodes
most = which(V(g2)$name == c("Alexis"))
least = which(V(g2)$name == c("Denno"))

V(g2)$name = c(1:length(V(g2)))
V(g2)$label = c(1:length(V(g2)))


#change the r=5 and s=10 for both most and least
InstaModel = simulate_sir(g2, simlength = 50, p.t=0.01,
                             removeafter = 2, susceptibleafter = 7,
                             display_net = F, startNode = most)
plotGraphs(InstaModel, vcount(g2))

InstaModel1= simulate_sir(g2, simlength = 50, p.t=0.01,
                             removeafter = 2, susceptibleafter = 7,
                             display_net = F, startNode = least)
plotGraphs(InstaModel1, vcount(g2))


##r=5, s=10

InstaModel2 = simulate_sir(g2, simlength = 50, p.t=0.01,
                          removeafter = 5, susceptibleafter = 10,
                          display_net = F, startNode = most)
plotGraphs(InstaModel2, vcount(g2))

InstaModel3= simulate_sir(g2, simlength = 50, p.t=0.01,
                          removeafter = 5, susceptibleafter = 10,
                          display_net = F, startNode = least)
plotGraphs(InstaModel3, vcount(g2))
