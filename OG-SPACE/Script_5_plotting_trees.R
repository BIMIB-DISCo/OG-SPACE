



do_geneaology_tree <-
  as.numeric(inputs[inputs[, 1] == "do_geneaology_tree", ][2])

do_phylo_tree <- as.numeric(inputs[inputs[, 1] == "do_phylo_tree", ][2])

size_of_points_trees  <- as.numeric(inputs[inputs[, 1] == "size_of_points_trees", ][2])

if (do_geneaology_tree == 1) {
  matrix_of_cell_geneaology <-
    read.table(paste(directory, "output/matrix_of_cell_geneaology.txt", sep =
                       "/"),
               sep = "\t", colClasses = "character")
  
  colnames(matrix_of_cell_geneaology) <- matrix_of_cell_geneaology[1, ]
  matrix_of_cell_geneaology <- matrix_of_cell_geneaology[-1, ]
  
  
  
  ### Igraph object
  for_graph_genealogy <- matrix_of_cell_geneaology[, c("father", "son")]
  
  
  cell_genealogy_graph <-
    graph_from_edgelist(as.matrix(for_graph_genealogy),
                        directed = T)
  
  ### COloring vertex
  
  
  color_vertex <-
    matrix(0, ncol = 3, nrow = length(names(V(
      cell_genealogy_graph
    ))))
  
  row.names(color_vertex)<-names(V(
    cell_genealogy_graph
  ))
  
  colnames(color_vertex) <- c("node","sub_pop","color")
  color_vertex <- as.data.frame(color_vertex)
  color_vertex$node <- row.names(color_vertex)
  
  
  for (i in 1:length(color_vertex$sub_pop)) {
    if(color_vertex$node[i]!="root"){
    color_vertex$sub_pop[i] <-
     unlist( matrix_of_cell_geneaology[matrix_of_cell_geneaology[, 2] ==
                                  color_vertex$node[i], ][4])

    }else{
      color_vertex$sub_pop[i] <- 0
    }
    }
  
  color_vertex$color <-as.numeric(unlist(color_vertex$sub_pop))
  color_vertex$color[color_vertex$color!=0]<- color_vertex$color[color_vertex$color!=0]+20
  
  color_vertex$color[color_vertex$color==0]<- color_vertex$color[color_vertex$color==0]+1
  color_vertex$color <- colors()[color_vertex$color]
  V(cell_genealogy_graph)$color <- color_vertex$color 
  ### SETTING Y COORDINATES AS ABS TIME
  
  coords <- layout_as_tree(cell_genealogy_graph)
  
  row.names(coords) <- names(V(cell_genealogy_graph))
  ### note the time is backward for the plot!
  time_of_root <- max(as.numeric(matrix_of_cell_geneaology[, 3])) + 1
  
  coords[rownames(coords) == "root", 2] <- time_of_root
  
  for (k in 1:nrow(coords)) {
    if (row.names(coords)[k] != "root") {
      coords[k, 2] <-
        time_of_root - as.numeric(matrix_of_cell_geneaology[matrix_of_cell_geneaology[, 2] == row.names(coords)[k] , 3])
    }
    
  }
  
  
  
  #### Plotting
  
  
  pdf(
    file = paste(
      directory,
      "output/plots/cell_genealogy_graph_with_branch_length.pdf",
      sep = "/"
    ),
    width = 2 * 19,
    height = 20
  )
  
  plot(
    cell_genealogy_graph,
    edge.arrow.size = 0.05,
    vertex.size =  size_of_points_trees,
    layout = coords,
    # layout=layout_as_tree,
    vertex.label = NA,
    vertex.color = V(cell_genealogy_graph)$color
  )
  
  dev.off()
  
  pdf(
    file = paste(directory, "output/plots/cell_genealogy_graph.pdf", sep = "/"),
    width = 2 * 19,
    height = 20
  )
  
  plot(
    cell_genealogy_graph,
    edge.arrow.size = 0.05,
    vertex.size =  size_of_points_trees,
    #layout=coords,
    layout = layout_as_tree,
    vertex.label = NA,
    vertex.color = V(cell_genealogy_graph)$color
  )
  
  dev.off()
  
}
#################################
if (do_phylo_tree == 1) {
  
  matrix_of_phylogenetic_relations <- read.table(paste(
    directory,
    "output/matrix_of_phylogenetic_relations.txt",
    sep = "/"
  ),
  sep = "\t", colClasses = "character")
  
  
  
  colnames(matrix_of_phylogenetic_relations) <-
    matrix_of_phylogenetic_relations[1, ]
  matrix_of_phylogenetic_relations <-
    matrix_of_phylogenetic_relations[-1, ]
  
  
  
  ### Igraph object
  
  for_graph_phylo <-
    as.matrix(matrix_of_phylogenetic_relations[, c("father", "son")])
  
  
  
  phylo_graph <- graph_from_edgelist(for_graph_phylo, directed = T)
  
  
 
  
  ### SETTING Y COORDINATES AS ABS TIME
  
  coords <- layout_as_tree(phylo_graph)
  
  row.names(coords) <- names(V(phylo_graph))
  ### note the time is backward for the plot!
  time_of_root <-
    max(as.numeric(matrix_of_phylogenetic_relations[, 3])) + 1
  
  root <-
    unique(matrix_of_phylogenetic_relations[!(matrix_of_phylogenetic_relations[, 1] %in% matrix_of_phylogenetic_relations[, 2]), 1])
  
  
  
  
  coords[rownames(coords) == root , 2] <- time_of_root
  
  for (k in 1:nrow(coords)) {
    if (row.names(coords)[k] != root) {
      coords[k, 2] <-
        time_of_root - as.numeric(matrix_of_phylogenetic_relations[matrix_of_phylogenetic_relations[, 2] == row.names(coords)[k] , 3])
      
    }
  }
  
  coords[, 1] <- coords[, 1] * 2
  
  
 
  
  
  
  ### COloring vertex
  
  
  color_vertex <-
    matrix(0, ncol = 3, nrow = length(names(V(
      phylo_graph
    ))))
  
  row.names(color_vertex)<-names(V(
    phylo_graph
  ))
  
  colnames(color_vertex) <- c("node","sub_pop","color")
  color_vertex <- as.data.frame(color_vertex)
  color_vertex$node <- row.names(color_vertex)
  
  root_real <-   unique(matrix_of_phylogenetic_relations[!(matrix_of_phylogenetic_relations[, 1] %in% matrix_of_phylogenetic_relations[, 2]), 1])

  
  
  for (i in 1:length(color_vertex$sub_pop)) {
    if(color_vertex$node[i]!="root" && color_vertex$node[i]!=root_real){
      color_vertex$sub_pop[i] <-
        unlist( matrix_of_phylogenetic_relations[matrix_of_phylogenetic_relations[, 2] ==
                                            color_vertex$node[i], ][4])
      
    }else{
      if(color_vertex$node[i]=="root" ){ color_vertex$sub_pop[i] <- 0 }else{ color_vertex$sub_pop[i] <- 1}
    }
  }
  
  color_vertex$color <-as.numeric(unlist(color_vertex$sub_pop))
  color_vertex$color[color_vertex$color!=0]<- color_vertex$color[color_vertex$color!=0]+20
  
  color_vertex$color[color_vertex$color==0]<- color_vertex$color[color_vertex$color==0]+1
  color_vertex$color <- colors()[color_vertex$color]
  V(phylo_graph)$color <- color_vertex$color
  #### Plotting
  
  pdf(
    file = paste(
      directory,
      "output/plots/phylogeny_of_samples_with_branch_length.pdf",
      sep = "/"
    ),
    width = 40,
    height = 20
  )
  
  plot(
    phylo_graph,
    edge.arrow.size = 0.05,
    vertex.size = size_of_points_trees,
    layout = coords,
    #   layout= layout_as_tree,
    vertex.label = NA,
    vertex.color = V(phylo_graph)$color
  )
  
  dev.off()
  
  
  pdf(
    file = paste(directory, "output/plots/phylogeny_of_samples.pdf", sep = "/"),
    width = 40,
    height = 20
  )
  
  plot(
    phylo_graph,
    edge.arrow.size = 0.05,
    vertex.size = size_of_points_trees,
    #  layout=coords,
    layout = layout_as_tree,
    vertex.label = NA,
    vertex.color = V(phylo_graph)$color
  )
  
  dev.off()
}
