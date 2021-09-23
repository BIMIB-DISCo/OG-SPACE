







do_bulk_exp <- as.numeric(inputs[inputs[, 1] == "do_bulk_exp", ][2])

do_sc_exp <- as.numeric(inputs[inputs[, 1] == "do_sc_exp", ][2])

sequencing_depth_bulk  <-
  as.numeric(inputs[inputs[, 1] == "sequencing_depth_bulk"
                    , ][2])

mean_coverage_cell_sc <-
  as.numeric(inputs[inputs[, 1] == "mean_coverage_cell_sc"
                    , ][2])

prob_reads_bulk <-
  as.numeric(inputs[inputs[, 1] == "prob_reads_bulk", ][2])

fn_rate_sc_exp   <-
  as.numeric(inputs[inputs[, 1] == "fn_rate_sc_exp", ][2])

fp_rate_sc_exp   <-
  as.numeric(inputs[inputs[, 1] == "fp_rate_sc_exp", ][2])

detection_thr_sc <-
  as.numeric(inputs[inputs[, 1] == "detection_thr_sc", ][2])

minimum_reads_for_cell <-
  as.numeric(inputs[inputs[, 1] == "minimum_reads_for_cell"
                    , ][2])

detected_vaf_thr <-
  as.numeric(inputs[inputs[, 1] == "detected_vaf_thr", ][2])
### Bulk experiment From :"Spatially constrained tumour growth affects
#the patterns of clonal selection and neutral drift in cancer genomic data"

if (do_bulk_exp == 1) {
  ## generating gt vaf
  gt_VAF <-
    as.data.frame(read.table(paste(
      directory, "output/vaf_mutation.txt", sep = "/"
    ),
    sep = "\t", colClasses = "character"))
  
  colnames(gt_VAF) <- as.character(gt_VAF[1, ])
  
  gt_VAF <- gt_VAF[-1, ]
  
  ## generating coverage for site
  coverage_sites <- rpois(dim(gt_VAF)[1] , sequencing_depth_bulk)
  gt_VAF$error_VF <- NA
  
  for (i in 1:dim(gt_VAF)[1]) {
    gt_VAF$error_VF[i] <- rbinom(1,
                                 min(as.numeric(gt_VAF$count[i]),
                                     coverage_sites[i]),
                                 prob_reads_bulk *
                                   as.numeric(gt_VAF$VF[i])) /
      min(max(as.numeric(gt_VAF$count)), coverage_sites[i])
    
  }
  
  
  pdf(
    file = paste(
      directory,
      "output/plots/hist_detected_vaf_with_errors.pdf",
      sep = "/"
    ),
    width = 10,
    height = 10
  )
  
  hist(as.numeric(gt_VAF$error_VF)[as.numeric(gt_VAF$error_VF) > detected_vaf_thr],
       breaks = 1 / 0.01)
  dev.off()
  
  write.table(
    gt_VAF,
    file = paste(directory, "output/detected_vaf_with_errors.txt", sep = "/"),
    sep = "\t",
    col.names = T,
    row.names = T
  )
  
  
  
  
  
} else{
  print("no bulk experiment required")
  
}



if (do_sc_exp == 1) {
  matrix_of_cell_genotype <-
    read.table(paste(directory, "output/genotype_of_samples.txt", sep =
                       "/"),
               sep = "\t", colClasses = "character")
  
  
  colnames(matrix_of_cell_genotype) <-  matrix_of_cell_genotype[1, ]
  matrix_of_cell_genotype <-  matrix_of_cell_genotype[-1, ]
  
  matrix_of_cell_genotype$mutations_with_errors <- ""
  # generating number of read for each cells
  
  list_of_mutations_1 <- paste(matrix_of_cell_genotype$mutations,
                               collapse = ' ')
  
  list_of_mutations <- unique(unlist(str_split(
    list_of_mutations_1, " ", n = Inf,
    simplify = FALSE
  )))
  
  list_of_mutation <- list_of_mutations[list_of_mutations != ""]
  ### Sc experiment
  
  
  ###loop over sampled cells
  for (j in 1:nrow(matrix_of_cell_genotype)) {
    ## generating reads for site
    
    gt_cell <- paste(matrix_of_cell_genotype$mutations[j],
                     collapse = ' ')
    
    gt_cell <- unique(unlist(str_split(
      gt_cell, " ", n = Inf,
      simplify = FALSE
    )))
    
    gt_cell <- gt_cell[gt_cell != ""]
    
    number_of_read_site <- rpois(length(list_of_mutation) ,
                                 mean_coverage_cell_sc)
    
    
    
    
    
    for (k in 1:length(list_of_mutation)) {
      if (number_of_read_site[k] > minimum_reads_for_cell) {
        ### FN
        if (list_of_mutation[k] %in% gt_cell) {
          if (detection_thr_sc < rbinom(1, number_of_read_site[k],
                                        1 - fn_rate_sc_exp) / number_of_read_site[k]) {
            matrix_of_cell_genotype$mutations_with_errors[j] <- paste(
              matrix_of_cell_genotype$mutations_with_errors[j],
              list_of_mutation[k],
              sep = " "
            )
            
            
          }
          
          
        } else{
          ### FP
          if (detection_thr_sc < rbinom(1, number_of_read_site,
                                        fp_rate_sc_exp) / number_of_read_site[k]) {
            matrix_of_cell_genotype$mutations_with_errors[j] <- paste(
              matrix_of_cell_genotype$mutations_with_errors[j],
              list_of_mutation[k],
              sep = " "
            )
            
            
            
          }
          
          
        }
        
        
        
        
        
        
      } else{
        # is NA
        
        matrix_of_cell_genotype$mutations_with_errors[j] <- paste(
          matrix_of_cell_genotype$mutations_with_errors[j],
          paste(list_of_mutation[k], "NA", sep = "_"),
          sep = " "
        )
      }
      
    }
    
  }
  
  
  
  
  
  
  
  write.table(
    matrix_of_cell_genotype[, c(1, 3)],
    file = paste(
      directory,
      "output/matrix_of_cell_genotype_with_errors.txt",
      sep = "/"
    ),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  
  
  
  
  
  
  
} else{
  print("no SC experiment required")
}
