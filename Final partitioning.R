source('clean up and first partition.r')
source('second partition.r')



# Return matrix of scores. Rows representing cases, columns representing partition. Padraic code.
get_scores_matrix <- function(data, partition, weights){
  vapply(partition, function(part) get_score(data, part, weights), numeric(nrow(data)))
}



# EM Algorithm
# Helper functions - Padraic code

# Returns list with scores for each part as vector elements
get_scores_list <- function(data, partition, weights){
  lapply(partition, function(x) get_score(data, x, weights))
}


# Util partition functions  -Padraic code
random_partition <- function(num_parts, num_snps){
  split(1:num_snps, sample(1:num_parts, num_snps, replace = TRUE))
}




# Update partition functions   -Padraic code
update_snp_partition <- function(cases, controls, case_partition){
  # Starting with a partition of the cases.
  # Calculate the frequencies of snps on each case partition:
  freqs <- sapply(case_partition, function(part){
    get_allele_freqs(cases[part,])
  })
  snps <- 1:ncol(cases)
  # Put the snps into the case partition on which they have the highest frequency.
  return(split(snps, max.col(freqs)))
}

get_fit2 <- function(cases, controls, case_partition, snp_partition, weights){
  case_scores <- data.frame(get_scores_list(cases, snp_partition, weights))
  control_scores <- get_scores_list(controls, snp_partition, weights)
  # TODO: mean vs median?
  # control_medians <- sapply(control_scores, median)
  #  relative_case_scores <- sweep(case_scores, 2, control_medians, '/')
  # fit <- mean(diag(sapply(case_partition, function(part){
  #  colSums(relative_case_scores[part,])
  #  })))
  #return(fit)
}

get_allele_freqs <- function(data){
  colSums(data, na.rm = TRUE)/(2*colSums(!is.na(data)))
}

# Update partition function   -Padraic code
update_case_partition <- function(cases, controls, weights, snp_partition){
  # Re-partition cases based on SNP partition
  case_scores <- get_scores_matrix(cases, snp_partition, weights)
  control_scores <- get_scores_matrix(controls, snp_partition, weights)
  # TODO: mean vs median?
  control_medians <- apply(control_scores, 2, median)
  relative_case_scores <- sweep(case_scores, 2, control_medians, '/')
  # (A) Add a "don't know" box:
  # know <- abs(relative_case_scores[,1] - relative_case_scores[,2]) > 0.1 # Only for 2 parts
  know <- apply(relative_case_scores, 1, max) > 1
  case_part_vector <- max.col(relative_case_scores[know,])
  case_partition <- split((1:nrow(cases))[know], case_part_vector)
  # (B) Without a don't know box:
  case_part_vector <- max.col(relative_case_scores)
  case_partition <- split(1:nrow(cases), case_part_vector)
  return(case_partition)
}

# EM Algorithm - random partition - Muireann code
#MF maybe not necessary to initialize snps.... :  snp_partition <- random_partition(rank, num_snps)
# Initial partition
# Partition 1 = high_cases
# Partition 2 = other_partition

initial_partition = list('1'=as.integer(high_case[,1]),'2'=as.integer(other_partition))
initial_partition

cases<- as.matrix(cases)
controls<-as.matrix(controls)

#### EM function
em_partition <- function(cases, controls, rank, weights, max_iter = 2000, nrun = 20){
  num_snps <- ncol(cases)
  output <- lapply(1:nrun, function(i){
    num_iter = 0
    # Initialise
    case_partition <- initial_partition
    # Iterate until convergence or until max_iter is reached.
    old_parts <- NULL
    while(!(identical(old_parts, case_partition) | (num_iter > max_iter))){
      num_iter <- num_iter + 1
      old_parts <- case_partition
      # Re-partition SNPs based on cases partition
      snp_partition <- update_snp_partition(cases, controls, case_partition)
      # Re-partition cases based on SNP partition
      case_partition <- update_case_partition(cases, controls, weights, snp_partition)
    }
    fit <- get_fit2(cases, controls, case_partition, snp_partition, weights)
    return(list(fit = fit, snp_partition = snp_partition))
  })
  
  
}


part_1 <-em_partition(cases , controls , 2 , weights, max_iter=3000, nrun=20)
part_1


#snp partitions (Muireann code)
snp_partition_1 <- part_1$snp_partition$`1`
snp_partition_1<-data.frame(snp_partition_1)

snp_partition_2 <-part_1[[2]]$snp_partition$`2`
snp_partition_2<-data.frame(snp_partition_2)

#Create a matrix with the column numbers representing the snp columns numbers generated from the algorithm
#This matrix contains 28 snp columns numbers
snp<-as.matrix(c("1",  "4",  "7",  "8",  "9", "10", "13", "15", "16", "17", "21", "22", "23", "25", "26", "29", "30", "32","33", "34", "35", "36", "40", "41", "50", "64", "65", "72"))


# Pull partitioned snps from cleaned data (this contains 86 snps). Muireann code
snps_clean <- (cases[,snp])
snps_clean


