# Efficiency model functions 
# Write a function to calculate alpha values
calculate.alpha <- function(strain, prot_expr_data, sequence_data, id_wt_seq) {
  id_strain <- which(prot_expr_data$strain == strain)
  id_wt_strain <- intersect(id_strain, id_wt)
  avg_wt_prot_expr_strain <- mean(prot_expr_data$value[id_wt_strain])
  alpha <- avg_wt_prot_expr_strain * sequence_data$elongation_times_sequences[id_wt_seq]
  return(alpha)
}

# Accuracy model functions
add.mean.fluo.col <- function(data, seq_col = "X...sequence_name...", strain_col = "strain...", fluo_col = "fluorescence.value..AU.") {
  # keep only columns seq_col, strain_col and fluo_col
  data <- data[,c(seq_col, strain_col, fluo_col)]

  # Compute mean fluorescence value for each (sequence, strain) group
  mean_fluo <- aggregate(data[[fluo_col]] ~ data[[seq_col]] + data[[strain_col]], data = data, FUN = mean)
  names(mean_fluo) <- c(seq_col, strain_col, "mean_fluorescence")
  
  # Add the mean fluorescence value for each row
  data <- merge(data, mean_fluo, by = c(seq_col, strain_col))

  # Rename columns in dataframe for better understanding
  names(data) <- c("sequence_name", "strain", "value", "mean_fluorescence")

  # Remove rows where sequence_name is pET28b_empty (they're controls and carry no info)
  data <- data[data$sequence_name != "pET28b_empty",]

  return(data)
}

# This function outputs the total nonsense error rate of an mRNA sequence 
# from the input codons of the sequence (argument: sequence) 
# and the nonsense error rates for each individual codon of said sequence.
# Nonsense error rates were computed by the FONSE model (Gilchrist, 2015) (file: "nonsense_error_rates_fonse.csv").
# All three stop codons are excluded.

translation.success.rate.of.mRNA <- function(seq_cand, nonsense_error_rates_ecoli, verbose = FALSE) {
  codons_of_sequence <- seq.string.to.cod.string(seq_cand)
  L <- length(codons_of_sequence)
  
  if (verbose) {
    cat("length candidate sequence is:", L, '\n')
  }
  
  # Exclude the last codon (stop codon)
  codons <- codons_of_sequence[1:(L-1)]
  
  # Identify positions of codons in the nonsense error rate table
  positions_in_table <- match(codons, names(nonsense_error_rates_ecoli))

  # Retrieve nonsense error rates for codons
  error_rate_mRNA <- nonsense_error_rates_ecoli[positions_in_table]

  if (verbose) {
    cat("Elongation times of the codons:", error_rate_mRNA, '\n')
  }
  
  # Calculate total translation success rate of the mRNA sequence
  total_success_rate <- prod(1-error_rate_mRNA)
  return(total_success_rate)
}

# This function adds the total success rate column to dataframe
add.success.rate.col <- function(data, nonsense_error_rates_ecoli){
  data$total_success_rate <- sapply(data$Sequence, translation.success.rate.of.mRNA, nonsense_error_rates_ecoli = nonsense_error_rates_ecoli)
  return(data)
}

# This function calculates the r value for a given strain by multiplying the wild type success rate its corresponding mean fluorescence value
calculate.r <- function(strain, data, success_rate_wt) {
  mean_fluo_wt <- data[data$sequence_name == "V015-wildtype" & data$strain == strain, "mean_fluorescence"][1]
  r <- success_rate_wt * mean_fluo_wt
  return(r)
}

add.r.column <- function(data, success_rate_wt) {
  r_k12 <- calculate.r("K12", data, success_rate_wt)
  r_bl21 <- calculate.r("BL21DE3", data, success_rate_wt)
  data$r <- ifelse(data$strain == "K12", r_k12, r_bl21)
  return(data)
}

energy.expenditure.of.mRNA <- function(seq_cand, roc_table, elongation_times_table) {
    codons_of_sequence <- seq.string.to.cod.string(seq_cand)
    L <- length(codons_of_sequence)
    
    # Exclude the last codon (stop codon)
    codons <- codons_of_sequence[1:(L-1)]
    
    # Identify positions of codons in the elongation times table
    positions_in_table <- match(codons, names(elongation_times_table))
    
    # Retrieve elongation times for codons
    elong_times <- elongation_times_table[positions_in_table]
    
    if (verbose) {
    cat("Elongation times of the codons:", elong_times, '\n')
    }
  
    # Calculate total elongation time
    total_elong_time <- sum(elong_times)

    # Identify positions of codons in the ROC table
    positions_in_table <- match(codons, names(roc_table))

    # Retrieve ROC for codons
    roc_codon <- roc_table[positions_in_table]

    # Calculate total roc
    total_roc <- sum(roc_codon)
    
    # Calculate energy expenditure of the mRNA sequence
    energy_expenditure <- total_elong_time * roc_codon
    
    return(energy_expenditure)
}

# function to compute the accuracy part of the equation 
product.of.q1s <- function(seq_cand, nonsense_error_rates_ecoli, verbose = FALSE) {
  codons_of_sequence <- seq.string.to.cod.string(seq_cand)
  L <- length(codons_of_sequence)
  
  if (verbose) {
    cat("length candidate sequence is:", L, '\n')
  }
  
  # Exclude the last codon (stop codon)
  codons <- codons_of_sequence[1:(L-1)]
  
  # Identify positions of codons in the nonsense error rate table
  positions_in_table <- match(codons, names(nonsense_error_rates_ecoli))

  # Retrieve nonsense error rates for codons
  error_rate_mRNA <- nonsense_error_rates_ecoli[positions_in_table]
  
  # Calculate total translation success rate of the mRNA sequence
  prod_q1 <- prod(1-error_rate_mRNA)^-1
  return(prod_q1)
}