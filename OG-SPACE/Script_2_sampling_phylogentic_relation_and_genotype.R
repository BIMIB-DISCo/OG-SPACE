


######## #### #### #### #### #### #### #### #### #### ####
#### PART 2
#### #### #### #### #### #### #### #### #### #### #### ####

#loading rdata from previous part
load(paste(directory, "output/events.Rdata", sep = "/"))
load(paste(directory, "output/last_state.Rdata", sep = "/"))



#### temporary random sampling of nodes

## Fixed parameters  reading them from external file

## lenght of the sequence (haploid to be modified)

genomic_seq_length <-
  as.numeric(inputs[inputs[, 1] == "genomic_seq_length", ][2])


# rate of neutral mut per site
neutral_mut_rate <-
  as.numeric(inputs[inputs[, 1] == "neutral_mut_rate", ][2])
## per generation per site


detected_vaf_thr <-
  as.numeric(inputs[inputs[, 1] == "detected_vaf_thr", ][2])


## number of samples


n_sample <- as.numeric(inputs[inputs[, 1] == "n_sample", ][2])



do_random_sampling <-
  as.numeric(inputs[inputs[, 1] == "do_random_sampling", ][2])

dist_sampling <- as.numeric(inputs[inputs[, 1] == "dist_sampling", ][2])


N_elements <- dim(new_state)[[2]]

### number of events
number_of_events <- length(events)


# selecting random the sampled nodes
if (do_random_sampling == 1) {
  sampled_nodes <-
    as.character(sample(colnames(new_state)[which(new_state != 0)]
                        , min(n_sample, length(
                          which(new_state !=
                                  0)
                        )), replace
                        = FALSE))
  #storing the leaf for part 3
  stored_leaf <- sampled_nodes
} else{
  load(paste(directory,"input/graph.Rdata", sep = "/"))
  
  ## renaming nodes
  V(g)$label <-  colnames(new_state)
  
  ##
  random_vertex <-
    as.character(sample(colnames(new_state)[which(new_state != 0)]
                        , 1 , replace
                        = FALSE))
  
  random_vertex <- V(g)[V(g)$label == random_vertex]
  
  
  vector_distance <- distances(
    graph = g,
    v = V(g),
    to = as.character(random_vertex),
    mode = c("all"),
    weights = NULL,
    algorithm = c("automatic")
  )
  
  ## sampling at certain distance
  
  sampled_nodes <-
    V(g)$label[(which(vector_distance <= dist_sampling))]
  # discardinf empty nodes
  sampled_nodes <-
    sampled_nodes[sampled_nodes %in% colnames(new_state)[which(new_state != 0)]]
  stored_leaf <- sampled_nodes
  
}

relations <- list()


count_of_relation <- 0

# for backward in time RECONSTRUCT THE GENEAOLOGY OF CELLS
for (k in (number_of_events):1) {
  if (as.character(events[[k]][4]) %in% sampled_nodes &&
      events[[k]][1]  != "D") {
    temp <- as.character(events[[k]][4])
    
    count_of_relation <- count_of_relation + 1
    
    
    
    
    
    relations[[count_of_relation]] <- c(
      as.numeric(events[[k]][6]),
      as.character(events[[k]][3]),
      as.character(events[[k]][4]),
      events[[k]][1],
      as.character(events[[k]][5])
    )
    
    
    sampled_nodes <- sampled_nodes[which(sampled_nodes != temp)]
    
    sampled_nodes <- c(sampled_nodes, as.character(events[[k]][3]))
    
    sampled_nodes <- unique(sampled_nodes)
    
    
  }
}






## generating edges list matrix format for igraph, it should be changed
mat_relation = matrix(0 , ncol = 5 , nrow = length(relations))
for (k in 1:length(relations)) {
  mat_relation[k, 1] = as.character(relations[[k]][[2]])
  mat_relation[k, 2] = as.character(relations[[k]][[3]])
  mat_relation[k, 3] = as.numeric(relations[[k]][[1]])
  mat_relation[k, 4] = as.character(relations[[k]][[5]])
  
}






colnames(mat_relation) <- c("father", "son", "abs_time", "sub_pop_son",
                            "branch_length")
#fake_root
root <-
  unique(mat_relation[!(mat_relation[, 1] %in% mat_relation[, 2]), 1])

for (i in 1:length(root)) {
  mat_relation <- rbind(mat_relation, c("root", root[i], 0, 1, 0))
}



cell_genealogy_graph <-
  graph_from_edgelist(mat_relation[, c(1, 2)], directed = T)




##########################
## transform timing to length



for (k in length(relations):1) {
  son <- mat_relation[k, 2]
  father <- mat_relation[k, 1]
  
  
  if (father %in% mat_relation[, 2]) {
    mat_relation[k, 5] <-
      as.numeric(mat_relation[k, 3]) - min(
        as.numeric(mat_relation[which(mat_relation[, 2] == father), 3]))
    
    
  } else{
    mat_relation[k, 5] <-  mat_relation[k, 3]
  }
  
  
}







#### mat_relation contains every information for a plot
# of the node/cell geneaology



###########
# Part 3
########
#### CONSTRUCTING PHYLOGENY OF SAMPLES




root <-
  unique(mat_relation[!(mat_relation[, 1] %in% mat_relation[, 2]), 1])



degree_nodes <-
  degree(
    cell_genealogy_graph,
    v = V(cell_genealogy_graph),
    mode = "out",
    loops = TRUE,
    normalized = FALSE
  )

# internal vertex
vertex_cell_phylo <-
  row.names(as.matrix(degree_nodes[degree_nodes == 0 |
                                     degree_nodes > 1 |
                                     row.names(as.matrix(degree_nodes)) %in% stored_leaf]))

num_degree_nodes <- as.numeric(degree_nodes[degree_nodes == 0 |
                                              degree_nodes > 1 |
                                              row.names(as.matrix(degree_nodes)) %in% stored_leaf])

leaf_cell_phylo <-
  row.names(as.matrix(degree_nodes[degree_nodes == 0]))
#######
counter_of_relation <- 0

relation_node_leaf <- list()

paths <- as.list(all_simple_paths(cell_genealogy_graph, from = root, to =
                                    leaf_cell_phylo))


for (i in 1:length(paths)) {
  nodes_of_interest <- names(paths[[i]])
  
  temp <- intersect(nodes_of_interest, vertex_cell_phylo)
  
  for (k in 2:length(temp)) {
    temp_father <- temp[k - 1]
    
    temp_son <- temp[k]
    
    counter_of_relation <- counter_of_relation + 1
    
    relation_node_leaf[[counter_of_relation]] <-
      c(temp_father, temp_son)
    
    
  }
  
}

relation_node_leaf <- unique(relation_node_leaf)

####### storing of leaf useful later
stored_stored_leaf <- stored_leaf


##################
#EDGE LIST for phylogenetic tree of the samples
#################

mat_phylo_relation = matrix(0 , ncol = 5 , nrow = length(relation_node_leaf))


for (k in 1:length(relation_node_leaf)) {
  if (length(as.character(relation_node_leaf[[k]])) > 1) {
    mat_phylo_relation[k, 1] = as.character(relation_node_leaf[[k]][[1]])
    mat_phylo_relation[k, 2] = as.character(relation_node_leaf[[k]][[2]])
  }
}

colnames(mat_phylo_relation) <-
  c("father", "son", "abs_time", "sub_pop_son",
    "branch_length")

# CLEANING FROM REPEATED EDGES AND DUMMY
mat_phylo_relation <- unique(mat_phylo_relation)

mat_phylo_relation <- mat_phylo_relation[mat_phylo_relation[1:dim(mat_phylo_relation)[1], 1] !=
                                           c("0"), ]


### graph object
cell_phylo_graph <-
  graph_from_edgelist(mat_phylo_relation[, c(1, 2)],
                      directed = T)

#### ASSOCIATE SUB_POP label TO NODES




for (i in 1:length(mat_phylo_relation[, 4])) {
  mat_phylo_relation[i, 4] <-  mat_relation[mat_relation[, 2] ==
                                              mat_phylo_relation[i, 2], 4]
  
  
}

mat_phylo_relation[mat_phylo_relation[, 1] == "root", 4] <- "0"




##### CALCULATE BRANCH LENGTH

for (i in 1:length(mat_phylo_relation[, 3])) {
  if (mat_phylo_relation[i, 1] != "root") {
    mat_phylo_relation[i, 3] <-  mat_relation[mat_relation[, 2] ==
                                                mat_phylo_relation[i, 2], 3]
  } else{
    mat_phylo_relation[i, 3] <- 0
  }
  
  
}


##########################
## transform timing to length



for (k in length(mat_phylo_relation[, 2]):1) {
  son <- mat_phylo_relation[k, 2]
  father <- mat_phylo_relation[k, 1]
  
  
  if (father %in% mat_phylo_relation[, 2]) {
    mat_phylo_relation[k, 5] <- as.numeric(mat_phylo_relation[k, 3]) -
      min(as.numeric(mat_phylo_relation[which(mat_phylo_relation[, 2] ==
                                                father), 3]))
    
    
  } else{
    mat_phylo_relation[k, 5] <-  mat_phylo_relation[k, 3]
  }
  
  
}






######
# DRIVERS MUTATIONAL TREE
######
# reconstructing the mutational tree for drivers
count_of_drivers_relation <- 0
time_drivers <- list()
mut_events <- list()

# for backward in time for times of driver events
for (k in (number_of_events):1) {
  if (events[[k]][1]  == "MUT") {
    count_of_drivers_relation <- count_of_drivers_relation + 1
    
    
    
    mut_events[[count_of_drivers_relation]] <- events[[k]]
    
    time_drivers[[count_of_drivers_relation]] <-
      as.numeric(events[[k]][6])
    
    
  }
}


##
### sKIPPING THIS PART IF IT IS PRESENT ONLY ONE DRIVER EVENT
####
if (count_of_drivers_relation > 1) {
  list_driv_relation <- list()
  
  count_observed_drivers_rel <- 0
  
  for (i in 1:length(mut_events)) {
    count_observed_drivers_rel <- count_observed_drivers_rel + 1
    
    list_driv_relation[[count_observed_drivers_rel]] <-
      c(
        as.character(mut_events[[i]][2]),
        as.character(mut_events[[i]][5]),
        as.character(mut_events[[i]][5]),
        as.numeric(mut_events[[i]][6])
      )
    
    
    
  }
  
  
  mat_driv_relation <-
    matrix(0 , ncol = 5 , nrow = length(list_driv_relation))
  
  
  for (k in 1:length(list_driv_relation)) {
    mat_driv_relation[k, 1] = as.character(list_driv_relation[[k]][[1]])
    mat_driv_relation[k, 2] = as.character(list_driv_relation[[k]][[2]])
    mat_driv_relation[k, 4] = as.numeric(list_driv_relation[[k]][[3]])
    mat_driv_relation[k, 3] = as.character(list_driv_relation[[k]][[4]])
    
  }
  
  ## adding root
  
  
  
  colnames(mat_driv_relation) <-
    c("father", "son", "abs_time", "sub_pop_son",
      "branch_length")
  
  
  
  ##########################
  ## transform timing to length
  
  
  
  for (k in length(mat_driv_relation[, 2]):1) {
    son <- mat_driv_relation[k, 2]
    father <- mat_driv_relation[k, 1]
    
    
    if (father %in% mat_driv_relation[, 2]) {
      mat_driv_relation[k, 5] <- as.numeric(mat_driv_relation[k, 3]) -
        min(as.numeric(mat_driv_relation[which(mat_driv_relation[, 2] ==
                                                 father), 3]))
      
      
    } else{
      mat_driv_relation[k, 5] <-  mat_driv_relation[k, 3]
    }
    
    
  }
  
  
  mat_temp <- mat_driv_relation[, c("father", "son")]
  
  driv_graph <- graph_from_edgelist(mat_temp, directed = T)
  
  
  coords <- layout_as_tree(driv_graph)
  
  row.names(coords) <- names(V(driv_graph))
  ### note the time is backward for the plot!
  time_of_root <- max(as.numeric(mat_driv_relation[, 3])) + 1
  
  coords[rownames(coords) == "1", 2] <- time_of_root
  
  for (k in 1:nrow(coords)) {
    if (row.names(coords)[k] != "1") {
      coords[k, 2] <-
        time_of_root - as.numeric(mat_driv_relation[mat_driv_relation[, 2] == row.names(coords)[k] , 3])
    }
    
  }
  
  
  color_vertex <-
    matrix(0, ncol = 2, nrow = length(names(V(
      driv_graph
    ))))
  
  row.names(color_vertex)<-names(V(
    driv_graph
  ))
  
  colnames(color_vertex) <- c("node","color")
  color_vertex <- as.data.frame(color_vertex)
  color_vertex$node <- row.names(color_vertex)
  color_vertex$color <- as.numeric(color_vertex$node) + 20
  
  color_vertex$color <- colors()[color_vertex$color]
  V(driv_graph)$color <- color_vertex$color 
  
  
  
  ### Plotting
  pdf(
    file = paste(directory, "output/plots/driv_graph.pdf", sep = "/"),
    width = 2 * 19,
    height = 10
  )
  
  plot(
    driv_graph,
    edge.arrow.size = 0.05,
    vertex.size = 5,
    layout = coords,
    vertex.label = NA,
    vertex.color = V(driv_graph)$color
  )
  dev.off()
  
  ####
  # GENERATING driver genotype of cells
  ####
  
  
  sub_pop_of_leaf <- list()
  
  for (k in 1:length(stored_stored_leaf)) {
    sub_pop_of_leaf[[k]] <- c(stored_stored_leaf[k],
                              mat_relation[stored_stored_leaf[[k]] ==
                                             mat_relation[, 2], 4])
    
    
    
  }
  
  
  path_from_root_drivers <-
    all_simple_paths(driv_graph, "1", to = V(driv_graph))
  
  genotype_drivers <- list()
  
  for (i in 1:length(path_from_root_drivers)) {
    genotype_drivers[[i]] <- (names(all_simple_paths(driv_graph, "1",
                                                     to = V(driv_graph))[[i]]))
    
    
  }
  
  genotype_drivers[[length(path_from_root_drivers)+1]] <- "1"
  
  
  sample_matrix_driv <-  matrix(as.numeric(0) ,
                                ncol = length(unique(unlist(genotype_drivers))) +
                                  1 ,
                                nrow = length(stored_stored_leaf))
  
  colnames(sample_matrix_driv) <-
    c("sample", unique(unlist(genotype_drivers)))
  
  sample_matrix_driv[, colnames(sample_matrix_driv) == "1"] <-
    as.numeric(1)
  
  
  ## creating the driver genotype as a matrix
  
  for (i in 1:length(stored_stored_leaf)) {
    sample_matrix_driv[i, 1] <- stored_stored_leaf[[i]]
    sub_pop <- as.character(sub_pop_of_leaf[[i]][2])
    for (k in 1:length(genotype_drivers)) {
      if (sub_pop == genotype_drivers[[k]][length(genotype_drivers[[k]])]) {
        for (l in 1:length(genotype_drivers[[k]])) {
          sample_matrix_driv[i, colnames(sample_matrix_driv)
                             == genotype_drivers[[k]][l]] <-
            as.numeric(1)
          
        }
        
        
        
      }
    }
    
    
  }
  
  
  
  
  mat_num <-
    matrix(as.numeric(sample_matrix_driv[, -1]),
           # Convert to numeric matrix
           ncol = ncol(sample_matrix_driv[, -1]))
  
  rownames(mat_num) <- sample_matrix_driv[, 1]
  colnames(mat_num) <- colnames(sample_matrix_driv[, -1])
  
  
  
  
} else{
  mat_num <-  matrix(as.numeric(1) ,
                     ncol = 1 ,
                     nrow = length(stored_stored_leaf))
  rownames(mat_num) <- stored_stored_leaf
  colnames(mat_num) <- "1"
}


if (count_of_drivers_relation == 1) {
  sample_matrix_driv <-  matrix(as.numeric(0) ,
                                ncol = 3,
                                nrow = length(stored_stored_leaf))
  
  colnames(sample_matrix_driv) <- c("sample", "1", "2")
  row.names(sample_matrix_driv) <- stored_stored_leaf
  
  sample_matrix_driv[, colnames(sample_matrix_driv) == "1"] <-
    as.numeric(1)
  
  sample_matrix_driv[row.names(sample_matrix_driv) %in%
                       mat_relation[mat_relation[, 4] !="1" , 2],
                     colnames(sample_matrix_driv) ==
                       "2"] <-  as.numeric(1)
  
#  mat_num <- sample_matrix_driv[, colSums(sample_matrix_driv) > 0]
  
  
}

################
## Attaching neutral mutations to cell geanology
################
## we suppose that every cell division every base as the prob neutra_mut_rate_to
# to mutate binomial prob

## associating to every branch of phylogenetic cell geanology
# a number of mutations
n_neutral_mut_for_branch <- list()



for (i in 1:nrow(mat_relation)) {
  n_neutral_mut_for_branch[[i]] <-
    rbinom(1, size = genomic_seq_length,
           prob = neutral_mut_rate)
  
}

mat_relation <- cbind(mat_relation, unlist(n_neutral_mut_for_branch))

colnames(mat_relation)[6] <- "neutral_mut_per_branch"


### generating the list of random mutation
## association labels of mut to edges
list_mut_branch <- list()

for (i in 1:nrow(mat_relation)) {
  if (as.numeric(mat_relation[i, 6]) > 0) {
    list_mut_branch[[i]] <-   stri_rand_strings(as.numeric(mat_relation[i, 6]), 5)
    
  } else{
    list_mut_branch[[i]] <- NA
  }
}




mat_relation <- cbind(mat_relation, list_mut_branch)

colnames(mat_relation)[7] <- "labels_neutral_mut_per_branch"



## genotype of samples neutral


sample_matrix_neut <- as.data.frame(matrix("",
                                           ncol = 2,
                                           nrow = length(V (
                                             cell_genealogy_graph
                                           ))),
                                    stringsAsFactors = FALSE)

colnames(sample_matrix_neut) <- c("sample", "mutations")
sample_matrix_neut$sample <- names(V (cell_genealogy_graph))


sample_matrix_neut <-
  sample_matrix_neut[sample_matrix_neut$sample != "root"
                     ,]


node_OI <- names(V(cell_genealogy_graph))


for (i in 1:length(node_OI)) {
  if (node_OI[i] != "root") {
    if (node_OI[i] %in% stored_stored_leaf) {
      nodes_of_interest <- node_OI[i]
      
    } else{
      paths_from_OI <-
        as.list(all_simple_paths(cell_genealogy_graph, from =
                                   node_OI[i],
                                 to = leaf_cell_phylo))
      
      nodes_of_interest <-  unique(names(unlist(paths_from_OI)))
    }
    
    list_of_mutation <- unlist(mat_relation[mat_relation[, 2] == node_OI[i], 7])
    
    
    if (!is.na(list_of_mutation[1])) {
      list_of_mutation <- paste(list_of_mutation, collapse = ' ')
      
      
      sample_matrix_neut$mutations[sample_matrix_neut$sample %in% nodes_of_interest] <-
        as.character(paste(
          sample_matrix_neut$mutations[sample_matrix_neut$sample %in%
                                         nodes_of_interest] ,
          as.vector(list_of_mutation),
          sep = " "
        ))
      
      
      
    }
  }
}



write.table(
  sample_matrix_neut,
  file = paste(directory, "output/sample_matrix_neut_unobs.txt", sep = "/")
  ,
  sep = "\t",
  col.names = T,
  row.names = F
)


sample_matrix_neut_leafs <-
  sample_matrix_neut[sample_matrix_neut$sample
                     %in% stored_stored_leaf , ]



write.table(
  sample_matrix_neut_leafs,
  file = paste(directory, "output/sample_matrix_neut.txt", sep = "/")
  ,
  sep = "\t",
  col.names = T,
  row.names = F
)

##################
## SINGLE CELL and VF GROUND TRUTH
#################
genotype_of_samples <- sample_matrix_neut_leafs

#binding driver and neutral genotype

for (i in ncol(mat_num):1) {
  genotype_of_samples$mutations[genotype_of_samples$sample
                                %in% labels(mat_num[, i] != 0)] <-
    paste(colnames(mat_num)[i],
          genotype_of_samples$mutations[genotype_of_samples$sample %in%
                                          labels(mat_num[, i] != 0)]
          , sep = " ")
  
}

mutation_present <- unique(unlist(str_split(
  unlist(genotype_of_samples$mutations),
  " ",
  n = Inf,
  simplify = FALSE
)))


VAF_DF <- as.data.frame(matrix(0, ncol = 2,
                               nrow = length(mutation_present)))


colnames(VAF_DF) <- c("mutation", "count")
VAF_DF$count <-
  table(unlist(str_split(
    unlist(genotype_of_samples$mutations),
    " ",
    n = Inf,
    simplify = FALSE
  )))


VAF_DF$mutation <- names(table(unlist(
  str_split(
    unlist(genotype_of_samples$mutations),
    " ",
    n = Inf,
    simplify = FALSE
  )
)))

VAF_DF$VF <- VAF_DF$count / length(stored_stored_leaf)
VAF_DF <- VAF_DF[VAF_DF$VF >= detected_vaf_thr, ]

VAF_DF <- VAF_DF[VAF_DF$mutation != "" ,]


pdf(
  file = paste(directory, "output/plots/hist_vf.pdf", sep = "/"),
  width = 10,
  height = 10
)

hist(VAF_DF$VF, breaks = 1 / 0.01)

dev.off()

#Saving files


write.table(
  mat_relation,
  file = paste(directory, "output/matrix_of_cell_geneaology.txt", sep = "/")
  ,
  sep = "\t",
  col.names = T,
  row.names = F
)

if (count_of_drivers_relation > 1) {
  write.table(
    mat_driv_relation,
    file = paste(directory, "output/mat_driv_relation.txt", sep = "/")
    ,
    sep = "\t",
    col.names = T,
    row.names = F
  )
}

write.table(
  mat_phylo_relation,
  file =
    paste(
      directory,
      "output/matrix_of_phylogenetic_relations.txt",
      sep = "/"
    )
  ,
  sep = "\t",
  col.names = T,
  row.names = F
)


write.table(
  genotype_of_samples,
  file = paste(directory, "output/genotype_of_samples.txt", sep = "/"),
  sep = "\t",
  col.names = T,
  row.names = F
)

write.table(
  VAF_DF,
  file = paste(directory, "output/vaf_mutation.txt", sep = "/"),
  sep = "\t",
  col.names = T,
  row.names = F
)
