# OG-SPACE
# Introduction 
Optimized Gillespie algorithm for simulating  Stochastic sPAtial models of Cancer Evolution (OG-SPACE) is a computational framework to simulate the spatial evolution of cancer cells and the experimental procedure of bulk and Single-cell DNA-seq experiments.
OG-SPACE relies on an optimized Gillespie algorithm for a large number of cells able to handle a variety of Birth-Death processes on a lattice and an efficient procedure to reconstruct the phylogenetic tree and the genotype of the sampled cells. 
More in detail, the dynamics is modeled by a stochastic multi-type Birth-Death (BD) process over an arbitrary lattice, with distinct possible interaction rules, and models the spatio-temporal dynamics of a tumor within a 2D or a 3D environment. 
 All cells can acquire and accumulate random mutations over time, which can either be passengers -- with no functional effects --, or drivers -- thus enhancing the birth rate of all descendants. 
After simulating the dynamics,
OG-SPACE mimics sequencing and variant calling of a portion of cells (e.g., a biopsy), at either the the bulk or the single-cell resolution, with the possibility of including experiment-specific errors in the simulated output data. 
 
OG-SPACE provides as outputs: 
- the state of the lattice at any time of the simulation;
- the Ground Truth (GT) genotype of the sampled cells;
- the  GT  Variant Allele Frequency (VAF) spectrum of the sampled cells;
- the GT phylogenetic  tree  of the sampled cells, in which the leaves represent the cells, whereas the internal nodes represent the most recent common ancestors;
- the mutational tree of the driver mutations (if present), where the nodes represent the mutations and edges model the accumulation temporal direction;
-  the noisy  genotype of the sampled cells obtained by simulating the errors of a sc-DNA-seq experiment;
-  the noisy VAF spectrum of the sampled cells obtained by simulating the errors of bulk DNA-seq experiment.



The inputs of  OG-SPACE are described below. 
 
# Implementation 
OG-SPACE employs the efficient strategy to decouple the simulation of the Birth-death process on a lattice and the genetic evolution of the sequence of the cells.
Here we present a brief description of the algorithms used to simulate such a system.

# OGA
To simulate the spatial dynamics of a multi-type BD process on a lattice via an OGA the following lists are needed:<img src="https://render.githubusercontent.com/render/math?math=\mathcal{V}_i"> the list of nodes occupied by the <img src="https://render.githubusercontent.com/render/math?math=i^\text{th}"> type of particle (i.e., subpopulation) present in the lattice and <img src="https://render.githubusercontent.com/render/math?math=\mathcal{N}_l"> the list of the neighbours of each node in the network.

Then, to evaluate the dynamics OG-SPACE employs the Algorithm 1.
![image](https://user-images.githubusercontent.com/43064628/138439742-cbd216ca-c88c-4f3b-ba66-3c74c7c030a3.png)


Where,<img src="https://render.githubusercontent.com/render/math?math=Nu()">  indicated the number of elements of a set ,<img src="https://render.githubusercontent.com/render/math?math=Z^d"> is the lattice, <img src="https://render.githubusercontent.com/render/math?math=T_{max}"> the time of the simulation, <img src="https://render.githubusercontent.com/render/math?math=\alpha">   the birth rate of the wild-type cells,  <img src="https://render.githubusercontent.com/render/math?math=\beta">  the death rate of the cells,  <img src="https://render.githubusercontent.com/render/math?math=\mu_\text{dri}"> the probabilty to acquire a new driver mutation, and  <img src="https://render.githubusercontent.com/render/math?math=\bar{\alpha}_\text{dri}"> and  <img src="https://render.githubusercontent.com/render/math?math=\sigma^2"> are the parameters of the distribution of the birth adavantage given by a driver mutation.   

The step "Evaluate if the event is a phantom event" is different for every contact rule included in the current implementation of OG-SPACE.
For the contact process, the algorithm checks if y is empty; in the voter model, if the state of x is different from the state of y; in the hierarchical voter model, if x$bears more driver mutations respect to y. If one of these conditions is true, then the event is flagged as not phantom.


# Generating the Phylogeny and the genotype of the sampled cells and the 

After the simulation,OG-SPACE generates the list <img src="https://render.githubusercontent.com/render/math?math=\mathcal{S}_{n_{\text{fin}}}"> of sampled cells composed by a user-selected number of randomly distributed cells or cells that falls in a circular (2D scenario)/spherical (3D scenario) region with a user-selected radius.
Then, OG-SPACE reconstructs the phylogenetic tree of the sampled cells and their genotypes by computing the \emph{tree of the genealogy of the sampled cells} <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}=(V,E)">. 
To reconstruct <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}=(V,E)">, OG-SPACE saves the following lists: <img src="https://render.githubusercontent.com/render/math?math=\mathcal{PA}_m">, i.e.,  the label of the parental cell of the <img src="https://render.githubusercontent.com/render/math?math=m^{\text{th}}">  event,  <img src="https://render.githubusercontent.com/render/math?math=\mathcal{DA}_m"> the list of the labels of the two nodes occupied, and  <img src="https://render.githubusercontent.com/render/math?math=\mathcal{T}_m">  the time.
OG-SPACE then applies the algorithm presented in Table Algorithm 2 to finally obtain <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}=(V,E)">

![image](https://user-images.githubusercontent.com/43064628/138438185-13ba6f62-9d6f-4d37-9b2f-29daa13639fd.png)


Therefore, by deleting all the nodes with degree equal to 2 of <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}"> and redrawing the edges between the remaining node coherently, OG-SPACE obtains the phylogenetic tree of the sampled cells <img src="https://render.githubusercontent.com/render/math?math=\mathcal{S}_{n_{\text{fin}}}">. 



It is important to note that, since the nodes of <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}"> with degree 3 are coalescent events,  they also represent the internal nodes of a standard phylogenetic tree, whereas the nodes with degree 1 are either the root or the leaves of such tree. 
Therefore, by deleting all the nodes with degree equal to 2 and redrawing the edges between the remaining node coherently, it is easy to obtain the phylogenetic tree of the sampled cells . 

As specified in the Background section, the large majority of mutations that can hit a given cell during its lifetime have no functional effect (i.e., they are passengers). 
From the computational perspective, it would be fallacious to explicitly consider the accumulation of such mutations during the simulations. 
However, the VAF spectrum generated by considering drivers only would be insufficient ad unrealistic. 
For these reasons, OG-SPACE implements an a posteriori attachment of neutral mutations to the cells of a simulation. 
Supposing that the Infinite Site Assumption holds, OG-SPACE assigns to each edge of <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}"> a number of mutations via a Bernoulli trial, where the number of trials is the length of the genome and the probability of success (i.e., the emergence of a new neutral mutation ) is defined via
a user-selected parameter <img src="https://render.githubusercontent.com/render/math?math=\mu_\text{neut}">. 
By associating a unique label for each mutation, it is then possible to retrieve the genotype of each sample by first enumerating the edges of the paths between the root and a leaf of the tree (i.e., a sampled cell) and then associating all the mutations present on the edges of such path to the cell.

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



