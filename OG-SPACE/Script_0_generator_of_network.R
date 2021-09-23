### Script_0_generator_of_network.R --


### Input
###################################################

### 2 or 3d lattice integer number

dimension <- as.numeric(inputs[inputs[, 1] == "dimension", ][2])


### Number of elements for  edge of the 2/3 D lattice

N_e <- as.numeric(inputs[inputs[, 1] == "N_e", ][2])


### Distance of interaction

dist_interaction <-
  as.numeric(inputs[inputs[, 1] == "dist_interaction", ][2])

##################################################


vector_for_graph <- replicate(dimension , N_e)


### Creating the network for the dynamics

### Number of nodes

N_elements <- N_e ^ dimension


### Generating igraph object

g <- make_lattice(vector_for_graph,
                  nei = dist_interaction,
                  dim = dimension,
                  directed = FALSE,
                  mutual = FALSE,
                  circular = F
                  )


### Save the graph object

save(g, file = paste(directory, "input/graph.Rdata", sep = "/"))


#### end of file -- Script_0_generator_of_network.R
