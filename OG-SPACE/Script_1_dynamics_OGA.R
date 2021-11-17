### Script_1_dynamics_OGA.R

### Load the graph

load(paste(directory, "input/graph.Rdata", sep = "/"))


#################################################
### Saving the Nearest Neighbour of each node
#################################################

N_elements <- length(V(g))

list_nn <- list()


for (i in 1:length(V(g))) {
  list_nn[[i]] <- as.vector(neighbors(g,
                                      as.vector(V(g)[i]),
                                      mode = c("all")))
}


####################################
### Fixed parameters of the dynamics, reading them from external file


### Input
###################################################

### Type of the dynamics

simulate_process <- inputs[inputs[, 1] == "simulate_process", ][2]

if( simulate_process == "contact" ) {
  print("Starting the simulation of the dynamics: contact process")
}

if ( simulate_process == "voter") {
  print("Starting the simulation of the dynamics: voter model")
}

if ( simulate_process == "h_voter") {
  print("Starting the simulation of the dynamics: hierachical voter model")
}


### 2 or 3d lattice integer number

dimension <- as.numeric(inputs[inputs[, 1] == "dimension", ][2])


### Total time of the simulation [Arb units]

Tmax <-  as.numeric(inputs[inputs[, 1] == "Tmax", ][2])


### Dynamics paramenters of wild type cells

### Birth rate [1/ time]

alpha <- as.numeric(inputs[inputs[, 1] == "alpha", ][2])

### Death rate [1/ time]
beta <- as.numeric(inputs[inputs[, 1] == "beta", ][2])

### Driver mutation rate [1/ time] max value 1
### if =>1 every birth event create a new mutation

driv_mut <- as.numeric(inputs[inputs[, 1] == "driv_mut", ][2])


### Average birth advantage per driver mutation

driv_average_advantage <-
  as.numeric(inputs[inputs[, 1] == "driv_average_advantage", ][2])


n_events_saving <-
  as.numeric(inputs[inputs[, 1] == "n_events_saving", ][2])


##############################
### Initial condition of the network, one WT cell in a random node
###############################

random_start <- as.numeric(inputs[inputs[, 1] == "random_start", ][2])

N_starting <-  as.numeric(inputs[inputs[, 1] == "N_starting", ][2])


CI <- matrix(0, ncol = N_elements)
colnames(CI) <- seq(1:N_elements)


if (random_start == 1) {
  random <-  sample(1:N_elements, N_starting)
  CI[1, random] <- 1
} else {
  node_to_start <- as.numeric(inputs[inputs[, 1] == "node_to_start", ][2])
  CI[1, node_to_start] <- 1
}


### Inizializating the birth advantage for populations with driver
### mutations

### TO DO COUNT SUB POP

alpha_driv <- list()
alpha_driv[[1]] <- alpha
n_of_driver_for_sub_pop <- list()
n_of_driver_for_sub_pop[[1]] <- 1

#########################################
##### Inizializating lists for the dynamics
#########################################

result <- list() # list for storing the state of the network (not necessary)
events <- list() # list of the events
timing <- list() # list of the time of the event
pop_dy <- list() # list for the plot of the dynamics
time_results <- list()

result[[1]] <- CI # Dummy first event
events[[1]] <- c("CI", 0, N_starting) # Dummy first event
timing[[1]] <- 0
pop_dy[[1]] <- c(0, 1)


### Temp variable

time <- 0
contatore <- 1
new_state <- matrix(0, ncol = N_elements)
old_state <- matrix(0, ncol = N_elements)
tot_pop <- N_starting # Number of elements occupied
tot_pop_old <- N_starting # Number of elements occupied


### Setting the names of the node, is not necessary for the dynamics
### but for counting the death events in the reconsctrucing node/cell
### genealogy

colnames(new_state) <- seq(1:N_elements)
colnames(old_state) <- seq(1:N_elements)


#################
### Sampling times
#################



######## #### #### #### #### #### #### #### #### #### ####
#### Part 1 forward dynamics
#### #### #### #### #### #### #### #### #### #### #### ####

count_of_events <- 0
partial_count <- 0
count_of_events_2 <- 0


###  Dynamics, while time is less than the max and the lattice is not empty


while (time < Tmax && tot_pop_old != 0) {
  
  count_of_events_2 <- count_of_events_2 + 1
  
  ## Loading last state of the network except for t == 0
  if (time == 0) {
    old_state <- CI
    
    partial_count <- partial_count + 1
    
    time_results[[partial_count]] <- time
    result[[partial_count]] <- old_state
    
  } else {
    old_state <- new_state
  }
  
  ## Number of occupied nodes
  
  tot_pop_old <- length(which(old_state != 0))
  
  ## Labels of the present populations
  
  populations_label  <-
    unique (as.vector(old_state)[old_state != 0])
  
  ## Evaluating the number of different subpopulation (phenotypes)
  ## present in the network
  
  num_of_pop <- length(populations_label)
  
  
  pop_num <- list() # Number of nodes occupied for each subpop.
  
  list_node_occ <- list() # List of labels of occupied nodes for each subpop.
  
  ## Evaluating the nodes that are able to make a birth move (i.e.,
  ## they have at least one free neighbor)
  
  possible_nodes <- list()
  
  n_occ <- 0
  
  specific_n_occ <- replicate(num_of_pop, 0)
  
  for (i in 1:num_of_pop) {
    list_node_occ[[i]] <- which(old_state == populations_label[i])
    
    pop_num[[i]] <- length(list_node_occ[[i]])
  }
  
  possible_nodes <- unlist(list_node_occ)
  
  ## Evaluating exp rate for the next event
  
  ## Counting birth
  
  alpha_tot <- list()
  
  for (i in 1:num_of_pop) {
    alpha_tot[[i]] <- alpha_driv[[i]] * length(list_node_occ[[i]])
  }
  
  alpha_fin = sum(unlist(alpha_tot)) # Why '=' instead of '<-' ?
  
  ## Counting death
  
  beta_fin = beta * tot_pop_old
  
  
  ## Rate
  
  lambda = alpha_fin + beta_fin
  
  
  
  ## Time of the next event
  
  new_time <- rexp(1, lambda)
  
  ## Updating time
  
  time <- time + new_time
  
  ## New state
  
  new_state <- old_state
  
  ## Generating random number
  
  random_number <- runif (1, min = 0, max = 1)
  
  ## MC move
  
  
  ###### Generation of the list of moves
  
  
  ## Renormalization of prob of single events
  
  birth = unlist(alpha_tot)/ lambda
  
  death = beta_fin / lambda

  # probability vector
  prob_vet <- append(unlist(birth),death)
  prob_cum <- cumsum( prob_vet)
  
  ### Moves
  #subpop that makes a move
  target_subpop <- min(which(random_number <= prob_cum))
  
  ### the event is a death
  if(target_subpop==num_of_pop+1){
    
    count_of_events <- count_of_events + 1
    
    random_node <-
      unlist(list_node_occ)[sample(1:length(unlist(list_node_occ)), 1)]
    
    new_state[1, random_node] <- 0
    
    ## The node label must be changed
    
    colnames(new_state)[random_node] <-
      paste(labels(new_state[, random_node]), "*", sep = "")
    
    events[[count_of_events]] <-
      c("D",
        new_state[1, random_node], 
        labels(new_state[1, random_node]),
        NA,
        time)
    
    
    
  }else{
    
    ## The event is a birth 
    
    
    ## selecting a random occupied node by the target sub pop
  
    if(length(list_node_occ[[target_subpop]])==1){
      
      random_node <- list_node_occ[[target_subpop]]
      
    }else{
      
      
      
      random_node <- sample(size=1,x=list_node_occ[[target_subpop]])    
      
    }
  
      
    
    ## the neighbors of the random node
    
    NN <- unlist(list_nn[random_node])
    
    
    ## random jump
    
    
    jump <- sample(x=NN,size= 1)
    
    
    ## Changing rule in function of the selected interaction rule
    
    if ( simulate_process == "voter" ) {
      rule <- new_state[1, jump] != new_state[1, random_node]
    }
    
    if ( simulate_process == "contact") {
      rule <- new_state[1, jump] == 0
    }
    
    if ( simulate_process == "h_voter") {
      if( new_state[1, jump] == 0){
        rule <- T
      }else{    
        rule <- as.numeric(
          n_of_driver_for_sub_pop[[new_state[1, jump]]]) < 
          as.numeric(n_of_driver_for_sub_pop[[new_state[1, random_node]]])
      }
    }
    
    
    
    if (rule == T) {
      
      
      count_of_events <- count_of_events + 2
      
      ## Asymmetric generation of driver mutation. If a birth occurs
      ## we have the probability driv_mut that ONE son recive a driver
      ## mutation
      
      random_number_mut <- runif (1, min = 0, max = 1)
      
      if (random_number_mut > driv_mut) {
        ## No new driver mut
        
        new_state[1, jump] <- new_state[1, random_node]
        
        
        ## The node label must be changed
        
        colnames(new_state)[random_node] <-
          paste(labels(new_state[, random_node]), 
                ".",
                sep = "")
        
        events[[count_of_events - 1]] <-
          c("B", 
            new_state[1, random_node], 
            labels(old_state)[2][[1]][random_node], 
            labels(new_state)[2][[1]][random_node], 
            new_state[1, random_node], 
            time
          )
        
        events[[count_of_events]] <-
          c("B", 
            new_state[1, jump], 
            labels(old_state)[2][[1]][random_node], 
            labels(new_state)[2][[1]][jump], 
            new_state[1, jump], 
            time
          )
      } else {
        ## New driver mut
        
        new_state[1, jump] <- length(alpha_driv) + 1
        
        ## Creating the new driver alpha gaussian advantage
        
        alpha_driv[[ new_state[1, jump]]] <-
          alpha_driv[[target_subpop]] +
          rnorm(1, mean = driv_average_advantage, 
                sd = driv_average_advantage / 2)
        ## 
        ## counting the nuber of mutations of new pop
        
        n_of_driver_for_sub_pop[[length(alpha_driv)]] <-
          as.numeric(n_of_driver_for_sub_pop[[new_state[1, random_node]]]) + 1
        
        ##         The node label must be changed
        
        colnames(new_state)[random_node] <-
          paste(labels(new_state[, random_node]), 
                "^",
                sep = "")
        
        
        events[[count_of_events]] <-
          c("MUT", 
            new_state[1, random_node] , 
            labels(old_state)[2][[1]][random_node], 
            labels(new_state)[2][[1]][jump], 
            new_state[1, jump], 
            time
          )
        
        events[[count_of_events - 1]] <-
          c("B", 
            new_state[1, random_node] , 
            labels(old_state)[2][[1]][random_node], 
            labels(new_state)[2][[1]][random_node], 
            new_state[1, random_node], 
            time
          )
      }
    }
    
  }
  

  ## Storing state
  ## Necessary for part 2
  
  timing[[count_of_events_2]] <- as.numeric(time)
  
  
  ## Counting the number of occupied nodes for each pop
  tot_pop == length(which(new_state != 0))
  
  
  if (count_of_events_2 %% n_events_saving == 0) {
    
    
    partial_count <- partial_count + 1
    
    for (i in 1:length(alpha_driv)) {
      list_node_occ[[i]] <- which(new_state == i)
      
      pop_num[[i]] <- length(list_node_occ[[i]])
    }
    
    pop_dy[[partial_count]] <- c(unlist(pop_num))
    
    time_results[[partial_count]] <- time
    result[[partial_count]] <- new_state
  }
  
  
  if (beta == 0 & tot_pop == N_elements) {
    time <- Tmax + 1
  }
} # End while dynamics


### Saving last state of dynamics

save(file = paste(directory, "output/last_state.Rdata", sep = "/"), 
     new_state)
save(file = paste(directory, "output/events.Rdata", sep = "/"),
     events)
save(file = paste(directory, "output/result.Rdata", sep = "/"),
     result)
save(file = paste(directory, "output/timing.Rdata", sep = "/"),
     timing)
save(file = paste(directory, "output/pop_dy.Rdata", sep = "/"),
     pop_dy)
save(file = paste(directory, "output/time_results.Rdata", sep = "/"), 
     time_results)

### end of file -- Script_1_dynamics_OGA.R
