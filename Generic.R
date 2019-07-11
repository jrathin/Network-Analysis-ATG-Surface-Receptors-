#Set path to where all related files are maintained
setwd("C:/Users/jeyaramr/Desktop/Arabidopsis thaliana - Network Analysis")

library(dplyr)
library(igraph)
library(readxl)
library(ggplot2)
library(grid)
library(plotly)
library(ggnetwork)
library(visNetwork)
library(d3heatmap)

#Given a graph and a subset of vertices, returns a list indexing the subset with all the vertices of the graph.
#Used for community plots
high_nodes = function(graph, list) 
{
  d = c()
  for(i in 1:length(V(graph)$name))
  {
    if(V(graph)$name[i] %in% list)
      d = c(d,2)
    else
      d = c(d,1)
  }
  return(d)
}
  
#Indexes the ATGs with known immune-role. Cutoff = 2 is for strong confidence.
k_list = function(graph, known, cutoff=2)
{
  d = c()
  for(i in 1:length(V(graph)$name))
  {
    if(V(graph)$name[i] %in% known$ATG[known$Code>=cutoff])
      d = c(d,2)
    else
      d = c(d,1)
  }
  return(d)
}

#Plot a Network. Gives edge categorisation too
visPlot2 = function(subgraph, communities = rep(1, length(E(subgraph))), nodesize = 0, v_comm = rep(1, length(V(subgraph))))
{
  nodes <- data.frame(id = V(subgraph)$name, group = v_comm)
  nodes$font.size<-20
  nodes$value = (nodesize)*3000 + 10
  edges <- data.frame(from = get.edgelist(subgraph)[,1], to = get.edgelist(subgraph)[,2], value = communities)
  edges$width = 10
  edges$color <- edges$value
  plt = visNetwork(nodes, edges, height = "600px")%>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(selectedBy = "group",
              highlightNearest = TRUE, 
               nodesIdSelection = TRUE) %>%
    visInteraction(keyboard = TRUE,
                   dragNodes = T, 
                   dragView = T, 
                   zoomView = T)
  return(plt)
}

#Plot a network, no edge categorisation
visPlot = function(subgraph, communities = rep(1, length(V(subgraph))), nodesize = 0)
{
  nodes <- data.frame(id = V(subgraph)$name, group = communities)
  nodes$font.size<-20
  nodes$value = (nodesize)*3000 + 10
  edges <- data.frame(get.edgelist(subgraph))
  colnames(edges)<-c("from","to")
  plt = visNetwork(nodes, edges, height = "600px")%>%
         visIgraphLayout(layout = "layout_nicely") %>%
         visOptions(selectedBy = "group", 
                    highlightNearest = TRUE, 
                    nodesIdSelection = TRUE) %>%
         visInteraction(keyboard = TRUE,
                        dragNodes = T, 
                        dragView = T, 
                        zoomView = T)
  return(plt)
}

#Arbitrary function for some random checks. Probably not used anywhere
check = function(df_base, df_restrict)
{
  l1 = intersect(df_restrict$Gene, df_base$Bait)
  l2 = intersect(df_restrict$Gene, df_base$Prey)
  li = unique(append(l1,l2))
  return(length(li))
}

#Plot, but doesnt use visNetwork library
m_plot = function(graph)
{
  comm = cluster_walktrap(graph, steps = 8)
  layout = layout_nicely(graph, niter=10)
  colours = rainbow(max(membership(comm)))
  plot(graph, layout = layout, vertex.size = 3, vertex.label = NA, vertex.color = colours[membership(comm)])
}

#Difference Plot, doesnt use visNetwork library
m_diff_plot = function(graph)
{
  layout = layout_nicely
  plot(graph, layout = layout, vertex.size = 3, vertex.color = "slategrey", vertex.label = NA, edge.color = c("black","red")[(E(graph)$type1 == 1)+1])
}

#Normalises the Data frame read from excel. Renames column headings
data_norm = function(df, base = FALSE)
{
  
  df = mutate_if(df, is.character, list(toupper))
  
  if (base == TRUE)
  {
    names(df) = c("Bait","Prey","Score")
  }
  
  else
  {
    names(df) = c("Gene", "Value")
  }
  
  return(df)
}

#Creates induced subgraph of a graph with specified subset. Also includes the vertices NOT present in the induced set
rest_subgraph_with_isolated = function(graph, df_base, df_restrict, cutoff = median(df_restrict$Value))
{
  l1 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Bait)
  l2 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Prey)
  li = unique(append(l1,l2))
  ld = setdiff(unique(append(df_base$Bait, df_base$Prey)), li)
  sg = induced_subgraph(graph,li) + vertices(ld)
  return(sg)
}

#Creates induced subgraph of a graph with specified subset
rest_subgraph = function(graph, df_base, df_restrict, cutoff = median(df_restrict$Value), simpl = F)
{
  l1 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Bait)
  l2 = intersect(df_restrict$Gene[df_restrict$Value>=cutoff], df_base$Prey)
  li = unique(append(l1,l2))
  sg = induced_subgraph(graph,li)
  if (simpl == T)
  {
    return(simplify(sg))
  }
  return(sg)
}

#Maps Gene name to ID
mapper = function(df, map)
{
  for (i in 1:nrow(map))
  {
    df$Bait[df$Bait == map$GeneName[i]] = map$GeneID[i]
    df$Prey[df$Prey == map$GeneName[i]] = map$GeneID[i]
  }
  return(df)
}

#Splits a long Dataframe into groups of 2, each with gene name and a data column
split_df = function(df)
{
  n = ncol(df)
  x = list()
  for (i in 2:n)
  {
    x[[i-1]] = data_norm(df[,c(1,i)])
  }
  return(x)
}

#Creates difference plots between 2 graphs. Edges in first but not in second are coloured black. The alternate is red
diff_plot = function(g1,g2)
{
  dg1 = difference(g1,g2)
  E(dg1)$type1 = 'black'
  dg2 = difference(g2,g1)
  E(dg2)$type2 = 'red'
  dg = union(dg1,dg2)
  temp = c()
  E(dg)$type1[is.na(E(dg)$type1)] = 'red'
  E(dg)$type2[is.na(E(dg)$type2)] = 'black'
  return(dg)
}

