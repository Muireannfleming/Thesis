library(ggplot2)
library(dplyr)


t <- read.table('WTCCC2_SCZ_PGC2_128SNPs_SZnarrow.ped')
map <- read.table('WTCCC2_SCZ_PGC2_128SNPs_SZnarrow.map')
lookup <- read.table('scz2.rep.128.txt', header = TRUE)

#PLINK format:
# Column 1: Family ID ('FID')
# Column 2: Within-family ID ('IID'; cannot be '0')
# Column 3: Within-family ID of father ('0' if father isn't in dataset)
# Column 4: Within-family ID of mother ('0' if mother isn't in dataset)
# Column 5: Sex code ('1' = male, '2' = female, '0' = unknown)
# Column 6: Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing_columns data if case/control)
# Column 7,8: first snp
# Column 9,10: second snp
# ...and so on...


# Prepping the data - Padraic code
rownames(t) <- t$V2 #  Individual ID
# Ignore family structure and sex
t <- t[, -c(1:5)]
# Remove missing_columns phenotypes (encoded as -9 or a 0)
t <- t[t$V6 != -9,]
t <- t[t$V6 != 0,]
# can now set 0 to NA
t[t == 0] <- NA
# The default is to encode an affected individual as 2 and unaffected as a 1.
# 2 => TRUE, #1 => FALSE
t$case <- t$V6 == 2
t$V6 <- NULL
# Convert factors back to strings for analysis
t <- data.frame(lapply(t, as.character), stringsAsFactors=FALSE)
# Calculate the number of SNPS by ignoring the phenotype column and divide by 2:
num_snps <- (ncol(t) - 1)/2
rownames(lookup) <- lookup$snpid
# I want to add lines to ordered_lookup from lookup
# I need to first set columns names
ordered_lookup <- lookup[0,]

snp_names <- as.character(map$V2)
snp_locations <- map$V4
# NB: The following is a bunch of hacks that might not work on a different dataset (not reusable!)
new_data <- data.frame(t$case)
for(i in 1:length(snp_names)){
  # Attempt to lookup the snp by name
  line <- lookup[snp_names[i],]
  if(is.na(line[1])){
    # Failed ID lookup. snp name is different in lookup than in map.
    # To fix this: use the snp location instead.
    line <- lookup[lookup$bp == snp_locations[i],]
  }
  l <- as.character(line$a1a2)
  ordered_lookup[nrow(ordered_lookup) +1,] <- line
  pair <- unlist(strsplit(l, ''))
  if(length(pair) == 2){
    # This means the SNP is a substitution
    # Select two columns that encode one snp
    col1 <- t[,2*i-1]
    col2 <- t[,2*i]
    converted <- rep(0, length(col1))
    # Let 1 encode a2
    converted <- converted + (col1 == pair[2])
    converted <- converted + (col2 == pair[2])
  }else{
    # Probably an insertion/deletion
    ref <- pair[1]
    # Find how long each snp value is
    col1 <- sapply(t[,2*i-1], function(x) nchar(x, keepNA = TRUE))
    col2 <- sapply(t[,2*i], function(x) nchar(x, keepNA = TRUE))
    converted <- rep(0, length(col1))
    if(ref == 'I'){
      # Insertion is the reference (a1).
      # 0 for the longer snp
      # 1 for the shorter snp
      # Let 1 encode a2  
      # NB it will not always be length 1
      converted <- converted + (col1 == 1)
      converted <- converted + (col2 == 1)
    }else{
      # Deletion is the reference (a1).
      # So 0 for the deletion (a1)
      # 1 for the insertion (longer snp)
      # NB: the following is terrible code:
      # NB: it will not always be length 3
      converted <- converted + (col1 == 3)
      converted <- converted + (col2 == 3)
    }
  }
  new_data[ncol(new_data) +1] <- converted
}
#Padraic code
# Clean up:
new_data$case <- new_data$t.case
new_data$t.case <- NULL
colnames(new_data) <- c(1:num_snps, 'case')
d <- new_data

# Clean up data:
reference <- ordered_lookup[,c('a1a2', 'frqa','frqu','p','or')]
# Lets make 1 be the risk allele (odds ratio > 1)
# Currently 1 is a2
# So this is correct if or < 1 for a1(reference)
# If not we need to invert (also invert frqa, frqu and take inverse of or column in reference)
wrong <- reference$or > 1 # Reference allele is risk allele (not what we want)
# 2 => 0, 1 => 1, 0 => 2:
d[,c(wrong,FALSE)] <- 2 - d[,c(wrong,FALSE)]
# Invert or's:
reference$or[wrong] <- 1/reference$or[wrong]
reference[wrong,c('frqa','frqu')] <- 1 - reference[wrong,c('frqa','frqu')]
risk <- data.frame(1-reference$frqa, 1-reference$frqu, 1/reference$or,reference$p)
colnames(risk) <- c('frqa','frqu','or','p')

write.table(d, file = 'real.d', quote=FALSE, row.names = FALSE)
cases <- subset(d, case == TRUE)[,1:num_snps]
controls <- subset(d, case == FALSE)[,1:num_snps]

# Unclean cases and controls
unclean_cases <- subset(d, case == TRUE)[,1:num_snps]
unclean_controls <- subset(d, case == FALSE)[,1:num_snps]


#Padraic code
# Deal with missing values:
missing_columns <- sort(colSums(is.na(cases)), decreasing = TRUE)
missing_rows <- sort(rowSums(is.na(cases)), decreasing = TRUE)
# plot(missing_columns, main = 'number of NA\'s per SNP in cases only')
# plot(missing_rows, main = 'number of NA\'s per row in cases only')
# Remove columns with more than 30 missing values
to_keep <- names(missing_columns[missing_columns < 30])
cases <- cases[, colnames(cases) %in% to_keep]
controls <- controls[, colnames(controls) %in% to_keep]
risk <- risk[rownames(risk) %in% to_keep,]
# Replace missing SNP values with columns mean for cases
replace_na_mean <- function(x){
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

# (A) Replace missing values with mean:
cases <- data.frame(apply(cases, 2, replace_na_mean))
controls <- data.frame(apply(controls, 2, replace_na_mean))

# (B) Remove rows with missing values
# cases <- na.omit(cases)
# controls <- na.omit(controls)

# Use matrices for efficiency.
# Muireann code
cases <- as.matrix(cases)
controls <- as.matrix(controls)
#colnames(cases)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86")
#colnames(controls)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86")
cases
controls

# Allele Frequencies of SNPs (cases)
# Muireann code
all_freq_cases<-apply(cases,2,mean)
all_freq_cases
all_freq_case<-data.frame(all_freq_cases)

#Filter frequencies
high_cases_freqs<-filter(all_freq_case, all_freq_cases > 1.5)
high_cases_freqs
low_cases_freqs<-filter(all_freq_case, all_freq_cases <1.5)
low_cases_freqs

#Get allele frequencies
get_allele_freqs <- function(data){
  colSums(data, na.rm = TRUE)/(2*colSums(!is.na(data)))
}



M <- length(to_keep)
# Clean up:
rm(col1, col2, converted, i, l, num_snps,
   pair, ref, snp_locations, snp_names, wrong, line,
   map, new_data, t, lookup, ordered_lookup, reference, d)

# odds_ratio <- get_odds_ratio(cases, controls) #  Calculate from the data directly. Padraic code.
odds_ratio <- risk$or #  Use odds ratios calculated from the full GWAS.
odds_ratio

get_allele_freqs <- function(data){
  colSums(data, na.rm = TRUE)/(2*colSums(!is.na(data)))
}
get_odds_ratio <- function(cases, controls){
  risk_case <- get_allele_freqs(cases)
  reference_case <- 1 - risk_case
  risk_controls <- get_allele_freqs(controls)
  reference_controls <- 1 - risk_controls
  return((risk_case * reference_controls) / (risk_controls * reference_case))
}

odd_ratio<-get_odds_ratio(cases, controls)
odd_ratio

# Convert to matrix. Muireann code.
OR_matrix <- as.matrix(odd_ratio)
OR_matrix

# Get the natural log (needed for weights). Muireann code. 
log_OR <- log(OR_matrix)
log_OR

# Calculate the weights. Muireann code.
weights<-as.matrix(log_OR[1:86,])
weights

# Helper functions are all Padraics code. 
get_score <- function(data, part, weights){
  B <- as.matrix(weights[part])
  X <- data[,part]
  return(as.vector(X %*% B))
  
}

#Get and print PRS of cases and controls. Padraic code.
score_cases <- get_score(cases, 1:86, weights)
score_cases<-data.frame(score_cases)

score_controls <- get_score(controls, 1:86, weights)
score_controls<-data.frame(score_controls)

#Convert into dataframe -  MF code.
controls_score<-data.frame(score_controls)

cases_score<-data.frame(score_cases)




#Plot graphs
plot_cases<-ggplot(data=score_cases, aes(x=score_cases))+geom_histogram(fill="grey", colour="black")
plot_cases
plot_cases+theme_classic()+geom_vline(xintercept=8.5, color="Red", size=.75)+ggtitle("PRS of cases")+ylab("Individuals")+xlab("PRS")                   
#plot_cases+theme_classic()+geom_hline(yintercept=8.5, color="Red", size=.75)+ggtitle("PRS of cases")+ylab("PRS")+xlab("Individuals")
plot_dist_cases<-ggplot(cases_score, aes(score_cases))+geom_density()+theme_classic()+geom_vline(xintercept=8.5, color="Red", size=.5)+xlab("PRS")+ggtitle("Density of PRS of Cases")
plot_dist_cases

plot_controls<-ggplot(data=score_controls, aes(score_controls))+geom_histogram(fill="grey", colour="black")
plot_controls+theme_classic()+ggtitle("PRS of controls")+ylab("Individuals")+xlab("PRS")
plot_dist_controls<-ggplot(controls_score, aes(score_controls))+geom_density()+theme_classic()+xlab("PRS")+ggtitle("Density of PRS of Controls")
plot_dist_controls

# Order cases
order_cases<-cases_dataframe[order(cases_dataframe$score_cases),]
order_cases




# Partition cases based on score (=>8.5) (Need dplyr for this!). Muireann code.
high_cases<-filter(order_cases, score_cases > 8.5)
high_cases<-data.frame(high_cases)
low_cases<-filter(order_cases, score_cases <8.5)
low_cases<-data.frame(low_cases)

# # Partition cases based on score (=>8.5) (Need dplyr for this!) - using unordered dataframe
# Muireann code
high_case<-filter(cases_dataframe, score_cases >= 8.5)
# high_case is used in to initialize the em algorithm, as the high scoring cases
high_case

