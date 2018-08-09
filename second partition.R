# Second partition 

source('clean up and first partition.r')


# Removing high. freq SNPs from top 10% of PRS in cases. Muireann code. 
high<-cases[high_case$X1.1574,]
high
low<-cases[-high_case$X1.1574,]
low


# Allele frequencies in high/low cases. Muireann code.
freq_high<- apply(high,2,mean)
freq_high

freq_low<-apply(low, 2, mean)
freq_low

# Quantiles - Using quantiles of the freq_high, as removal of all freq_high results in just 12 remaining snps
q_frequencies<-quantile(freq_high)
q_frequencies

quantile(freq_low - freq_high,0.9)
quantile_threshold = quantile(freq_low - freq_high, 0.5)

# New cases - Remove SNPs that have f.high > f.low in the low partition. MF code.
new_cases<-low[,freq_low - freq_high >= quantile_threshold]
#new_cases<-as.integer(new_cases)



# Re-calculate PRS 
# natural log or OR from above
weights<-as.matrix(log_OR[1:86,])
weights

colnames(new_cases)<-c("10","13","14","16","17","19","27","29","33","34","38","39","41","46","52","55","56","57","59","60","61","62","68","70","75","77","82","84","96","98","99","100","101","102","104","105","108","110","114","115","119","120","123")
#new_cases<-as.integer(new_cases)

part<- as.integer(colnames(new_cases))
data<-low     

# Helper functions are all Padraics code. 
get_score <- function(data, part, weights){
  B <- as.matrix(weights[part])
  X <- data[,part]
  return(as.vector(X %*% B))
  
}



#= used to pull the subset of snps needed.  
new_cases_score<-get_score(low, part, weights)
#new_cases_score<-get_score(low, 1:43, weights)
new_cases_score <-data.frame(new_cases_score)
new_cases_score

# Plot distribution
plot_dist<-ggplot(new_cases_score, aes(new_cases_score))+geom_density()+theme_classic()+xlab("PRS")+ggtitle("Density of PRS of Newly Partitioned Cases")
plot_dist

# Plot scores 
plot_cases<-ggplot(data=new_cases_score, aes(x=new_cases_score))+geom_histogram(fill="grey", colour="black")
plot_cases+theme_classic()+ggtitle("PRS of Newly Partitioned Cases")+ylab("Individuals")+xlab("PRS")

# Convert to a vector
#new_cases_score <-as.vector(new_cases_score) 
#low<-as.vector(low)

# Other partition. Muireann code.
q_threshold2 = 0.9
new_cases_score<-as.matrix(new_cases_score)

other_partition = rownames(low)[new_cases_score > quantile(new_cases_score,q_threshold2)]
other_partition 

