# OG-SPACE
# Introduction 
Optimized Gillespie algorithm for simulating  Stochastic sPAtial models of Cancer Evolution (OG-SPACE) is a computational framework to simulate the spatial evolution of cancer cells and the experimental procedure of bulk and Single-cell DNA-seq experiments.
OG-SPACE relies on an optimized Gillespie algorithm for a large number of cells able to handle a variety of Birth-Death processes on a lattice and an efficient procedure to reconstruct the phylogenetic tree and the genotype of the sampled cells. 


# REQUIRED  SOFTWARE AND PACKAGE

- R (tested on version 4.0) https://cran.r-project.org
- The following R libraries:
  - igraph
  - gtools
  - ggplot2
  - gridExtra
  - reshape2
  - stringi
  - stringr
  - shiny
  - manipulateWidget
  - rgl

# RUN OG-SPACE
- Download the folder OG-SPACE.
- use the following command "Rscript.exe  my_path\Run_OG-SPACE.R". "my_path" is the path to the folder containing the OG-SPACE scripts.
- When the pop-up window appears, select the file "Run_OG-SPACE.R" in the working folder.
Alternatively, you can launch OG-SPACE, with software like RStudio. In this case, simply run the script "Run_OG-SPACE.R" and when the pop-up window appears, select the file "Run_OG-SPACE.R" in the working folder.
 
# PARAMETERS OF OG-SPACE
Most of the parameters of OG-SPACE could be modified by editing with a text editor the file "input/Parameters.txt". Here a brief description of each parameters.
- simulate_process            three values "contact","voter" and "h_voter". This parameter selects which model simulate with OG-SPACE.
- generate_lattice = if 1 OG-SPACE generate a regular lattice for the dynamics. If 0 OG-SPACE takes an Igraph object named "g.Rdata" in the folder "input".
- dimension = an integer number, the dimensionality of the generated regular lattice.
- N_e  			   =  an integer number, number of elements of the edge of the generated regular lattice.
- dist_interaction = an integer number, the distance of interaction between nodes of the lattice.
- simulate_experiments       = if 1 OG-SPACE generates bulk and sc-DNA seq experiments data. If 0, no. 
- do_bulk_exp		     = if 1 OG-SPACE generates  bulk seq experiment data . If 0,  no
- do_sc_exp		     = if 1 OG-SPACE generates sc-DNA seq experiments data . If 0, no
- to_do_plots_of_trees         = if 1 OG-SPACE generates the plots of the trees . If 0, no. 
- do_pop_dyn_plot              = if 1 OG-SPACE generates the plots of the dynamics . If 0, no. 
- do_spatial_dyn_plot         = if 1 OG-SPACE generates the plots of the spatial dynamics . If 0, no.  
- do_geneaology_tree           = if 1 OG-SPACE generates the plots of the cell genealogy trees . If 0,  no. 
- do_phylo_tree              =  if 1 OG-SPACE generate the plots of the phylogenetic trees . If 0  no. 
- size_of_points_lattice       = an integer number, size of the points in the plot of spatial dynamics.
- size_of_points_trees         = an integer number, size of the points in the plot of trees.
- set_seed                    = the random seed of the computation.
- Tmax  = maximum time of the computation [arb. units] . 
- alpha  = birth rate of the first subpopulation [1/time].
- beta     = death rate of the first subpopulation [1/time].
- driv_mut      =   probability of developing a driver mutation (between 0 and 1).
- driv_average_advantadge = average birth rate advantage per driver [1/time].
- random_start               = if 1 OG-SPACE select randomly the spatial position of the first cell . If 0 it use the variable "node_to_start" .
- node_to_start              = if  random_start=0 OG-SPACE, the variable should be setted to the label of the node of starting.
- N_starting                 = Number of starting cells.  Works only with  random_start=1.
- n_events_saving             = integer number, frequency of the number of events when saving the dynamics for the plot.
- do_random_sampling           = if 2 OG-SPACE samples randomly the cells.
- -n_sample		    = integer number of the number of sampled cell. Ignored if do_random_sampling  = 0
- dist_sampling                = The radius of the spatial sampled region. Ignored if do_random_sampling  = 1
- genomic_seq_length          = number of bases of the genome under study.
- neutral_mut_rate            = neutral mutational rate per base [1/time].
- n_time_sample               = integer number, number of the plots of the dynamics.
- detected_vaf_thr 	    = VAF threshold. If a VAF is lesser than this number is considered not observed.
- sequencing_depth_bulk        = integer number, the sequencing depth of bulk sequencing.
- prob_reads_bulk		    = number between 0 and 1, 1- the prob of a false negative in bulk read
- mean_coverage_cell_sc      = integer number, mean number of read per cells
- fn_rate_sc_exp               = number between 0 and 1, 1- the prob of a false negative in sc read
- fp_rate_sc_exp              = number between 0 and 1, 1- the prob of a false positive in sc read
- minimum_reads_for_cell      = integer number, the minimum number of reads per cell in order to call a mutation
- detection_thr_sc             = ratio of successful reads necessary to call  a mutation

# OUTPUTS OF OG-SPACE

In the folder "output", you will find all the .txt data files of the output. Note that the trees are returned as edge list matrices.
The files will contain:
- The state of the lattice, with the position of each cell.
- The Ground Truth (GT) genotype of the sampled cells.
- The  GT  Variant Allele Frequency (VAF) spectrum of the sampled cells.
- The GT genealogy tree  of the sampled cells.
- The GT phylogenetic  tree  of the sampled cells.
- The mutational tree of the driver mutations appeared during the simulation of the dynamics.
-  The  genotype of the sampled cells after simulating a sc-DNA-seq experiment (if required).
- The VAF spectrum of the sampled cells after simulating a  bulk DNA-seq experiment (if required).
      
In the folder "output/plots", you will find all required plots.
