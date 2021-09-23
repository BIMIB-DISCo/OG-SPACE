

load(file = paste(directory, "output/last_state.Rdata", sep = '/'))

load(paste(directory, "output/pop_dy.Rdata", sep = '/'))

load(paste(directory, "output/timing.Rdata", sep = '/'))

load(file = paste(directory, "output/time_results.Rdata", sep = '/'))

load(paste(directory, "input/graph.Rdata", sep = '/'))

load(paste(directory, "output/result.Rdata", sep = '/'))










##################
# input

n_time_sample = min(as.numeric(inputs[inputs[, 1] == "n_time_sample", ][2]), length(result))
dimension <- as.numeric(inputs[inputs[, 1] == "dimension", ][2])
do_spatial_dyn_plot <- as.numeric(inputs[inputs[, 1] == "do_spatial_dyn_plot", ][2])
size_of_points_lattice     <- as.numeric(inputs[inputs[, 1] == "size_of_points_lattice", ][2])  

simulate_process <-  as.character(inputs[inputs[, 1] == "simulate_process", ][2])  

############### Plot of the graph evolution
if (do_spatial_dyn_plot == 1) {
  
  
  
  if(dimension == 2){
  if (length(result) > 2) {
    seq_of_time_sample <- seq(floor(length(result) / n_time_sample),
                              from = 1,
                              to = length(result) - 1)
    
    
    
    for (i in 1:length(seq_of_time_sample)) {

      temp <- as.numeric(result[[seq_of_time_sample[i]]])
      temp[temp!=0] <- temp[temp!=0] + 20
      temp[temp==0] <-temp[temp==0] + 1
      palette_node <-  colors()[temp]
      V(g)$color <- palette_node
      
      pdf(
        file = paste(directory,paste(  "output/plots/dynamics_at_time",
          as.character(time_results[[seq_of_time_sample[i]]]),"_model_",simulate_process,
          ".pdf",
          sep = ""
        ), sep="/" ),
        width = 10,
        height = 10
      )
      
      plot(
        g,
        vertex.label = NA,
        edge.arrow.size = 0.02,
        vertex.size = size_of_points_lattice,
        vertex.color = V(g)$color,
        layout = layout_on_grid(g,dim=dimension)
      )
      
      dev.off()
      
      
    }
  }
  

  temp <- as.numeric(new_state)
  temp[temp!=0] <- temp[temp!=0] + 20
  temp[temp==0] <-temp[temp==0] + 1
  palette_node <-  colors()[temp]
  V(g)$color <- palette_node
  
  
  pdf(
    file = paste(directory,paste("output/plots/dynamics_at_final_time.pdf","_model_",simulate_process,".pdf",sep=""), sep =
                   "/"),
    width = 10,
    height = 10
  )
  
  plot(
    g,
    vertex.label = NA,
    edge.arrow.size = 0.02,
    vertex.size = size_of_points_lattice,
    vertex.color = V(g)$color,
    layout = layout_on_grid(g,dim=dimension)
  )
  
  dev.off()
  
  
  
  
  temp <- as.numeric(result[[1]])
  temp[temp!=0] <- temp[temp!=0] + 20
  temp[temp==0] <-temp[temp==0] + 1
  palette_node <-  colors()[temp]
  V(g)$color <- palette_node  
  
  pdf(
    file = paste(directory, paste("output/plots/dynamics_intial","_model_",simulate_process,".pdf",sep=""), sep = "/"),
    width = 10,
    height = 10
  )
  
  plot(
    g,
    vertex.label = NA,
    edge.arrow.size = 0.02,
    vertex.size =   size_of_points_lattice,
    vertex.color = V(g)$color,
    layout = layout_on_grid(g,dim=dimension)
  )
  
  dev.off()
  }
  
  
  
if(dimension == 3){
  
  #size of the window of the plot
  r3dDefaults$windowRect <- c(0,50, 800, 800) 
  
  if (length(result) > 2) {
    seq_of_time_sample <- seq(floor(length(result) / n_time_sample),
                              from = 1,
                              to = length(result) - 1)
    
    
    
    for (i in 1:length(seq_of_time_sample)) {
      
     
      temp <- as.numeric( result[[seq_of_time_sample[i]]])
      temp[temp!=0] <- temp[temp!=0] + 20
      temp[temp==0] <-temp[temp==0] + 1
      palette_node <-  colors()[temp]
      V(g)$color <- palette_node
      
      
      rgl.viewpoint(theta = 35,phi = 45)
      
 
      plot3d(layout_on_grid(g,dim=3),col=as.character(V(g)$color),xlab = NULL,ylab = NULL,zlab = NULL,size = 7)
      
      rgl.postscript(paste(directory,paste(  "output/plots/dynamics_at_time",
                                             as.character(time_results[[seq_of_time_sample[i]]])
                                             ,
                                             ".pdf",
                                             sep = ""
      ), sep="/" ), fmt = 'pdf', drawText = F)
      
      
    }
  }
  
  
  temp <- as.numeric(new_state)
  
  temp[temp!=0] <- temp[temp!=0] + 20
  temp[temp==0] <-temp[temp==0] + 1
  palette_node <-  colors()[temp]
  V(g)$color <- palette_node
  

  
  rgl.viewpoint(theta = 35,phi = 45)
  
  plot3d(layout_on_grid(g,dim=3),col=as.character(V(g)$color),xlab = NULL,ylab = NULL,zlab = NULL,size = 7.5)
  rgl.postscript(paste(directory, "output/plots/dynamics_at_final_time.pdf", sep =
                         "/"), fmt = 'pdf', drawText = F)
  
  
  
  
  temp <- as.numeric(result[[1]])
  temp[temp!=0] <- temp[temp!=0] + 20
  temp[temp==0] <-temp[temp==0] + 1
  palette_node <-  colors()[temp]
  V(g)$color <- palette_node
  
  rgl.viewpoint(theta = 35,phi = 45)
  
  plot3d(layout_on_grid(g,dim=3),col=as.character(V(g)$color),xlab = NULL,ylab = NULL,zlab = NULL,size = 7.5 ) 
  rgl.postscript(paste(directory, "output/plots/dynamics_initial_time.pdf", sep =
                         "/"), fmt = 'pdf', drawText = F)
  
}
  
  
  }
