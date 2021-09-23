### Run_OG-SPACE.R

### Libraries

library(igraph)
library(gtools)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(stringi)
library(stringr)
library(shiny)
library(manipulateWidget)
library(rgl)

#directory<-"C:/Users/fabri/Desktop/V01_simulator_contact_process/V08"



FILE <- file.choose()

directory  <- dirname(FILE)


inputs <- read.table(paste(directory,
                           "input/Parameters.txt",
                           sep = "/"),
                     colClasses = "character")

generate_lattice <- as.numeric(inputs[inputs[, 1] == "generate_lattice", ][2])

seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])
  
  set.seed(seed)

if (generate_lattice == 1) {
  print("Generating the lattice")
  source(paste(directory, "Script_0_generator_of_network.R", sep = "/"))
}

rm(list = setdiff(ls(), c("inputs", "directory")))

simulate_process <- inputs[inputs[, 1] == "simulate_process", ][2]

seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])

set.seed(seed)

source(paste(directory, "Script_1_dynamics_OGA.R", sep = "/"))

rm(list = setdiff(ls(), c("inputs", "directory")))


seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])

set.seed(seed)

print("sampling phylogentic relation and genotype")


source(paste(directory,
             "Script_2_sampling_phylogentic_relation_and_genotype.R",
             sep = "/"))

rm(list = setdiff(ls(), c("inputs", "directory")))


seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])

set.seed(seed)

do_spatial_dyn_plot <-
  inputs[inputs[, 1] == "do_spatial_dyn_plot", ][2]

if ( do_spatial_dyn_plot == 1  ) {  
  print("Plotting the dynamics")
  source(paste(directory, "Script_3_saving_dynamics_plot.R", sep = "/"))
}

rm(list = setdiff(ls(), c("inputs", "directory")))



seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])

set.seed(seed)

do_bulk_exp <- inputs[inputs[, 1] == "do_bulk_exp", ][2]
do_sc_exp <- inputs[inputs[, 1] == "do_sc_exp", ][2]


if (do_sc_exp == 1 || do_bulk_exp == 1 ) {
  print("Simulating the experimental noise")
  source(paste(directory, "Script_4_simulating_bulk_sc_exp.R", sep = "/"))
}

rm(list = setdiff(ls(), c("inputs", "directory")))


seed <- as.numeric(inputs[inputs[, 1] == "set_seed", ][2])

set.seed(seed)

do_geneaology_tree <- inputs[inputs[, 1] == "do_geneaology_tree", ][2]
do_phylo_tree <- inputs[inputs[, 1] == "do_phylo_tree", ][2]


if (do_geneaology_tree == 1 || do_phylo_tree == 1 ) {
  print("Plotting trees")
  source(paste(directory, "Script_5_plotting_trees.R", sep = "/"))
}

rm(list = ls())


#### end of file -- Run_OG-SPACE.R
