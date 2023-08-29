library(dplyr)
library(poLCA)
library(ggplot2)

###########user part#################
user = read.csv("./Data_Survey_user_0813.csv")


convert_to_integer <- function(data) {
  mapping_df <- data.frame()
  
  for (col in colnames(data)) {
    #if (is.character(data[[col]])) {
      unique_strings <- unique(data[[col]])
      unique_values <- seq_along(unique_strings)
      
      mapping_df <- rbind(mapping_df, data.frame(Var = col,String = unique_strings, Value = unique_values))
      
      data[[col]] <- as.integer(factor(data[[col]], levels = unique_strings))
    #}
  }
  
  return(list(data = data, mapping = mapping_df))
}

convert_to_original <- function(converted_data,mapping) {
  data <- converted_data
  mapping_df <- mapping
  data$Var2=as.character(data$Var2)
  
  for (i in seq_len(nrow(data))) {
    col <- data$Var2[i]
    col = as.integer(gsub("Pr\\((\\d+)\\)", "\\1", col))
    value <- data$L2[i]
    
    original_value <- mapping_df$String[match(paste(value, col, sep = "-"), paste(mapping_df$Var,mapping_df$Value, sep = "-"))]
    #cat(original_value)
    data$Var2[i] <- original_value
  }
  
  return(data)
}


user_int = convert_to_integer(user)$data
user_int = user_int[,-c(1)]

f = as.matrix(user_int)~1

max_II <- -100000
min_bic <- 100000
for(i in 2:10){
  lc <- poLCA(f, user_int, nclass=2, maxiter=3000, 
              tol=1e-5, na.rm=FALSE,  
              nrep=10, verbose=TRUE, calc.se=TRUE)
  if(lc$bic < min_bic){
    min_bic <- lc$bic
    LCA_best_model<-lc
  }
}   

LCA_best_model

#AIC(2): 7619.417
#BIC(2): 8256.157
#G^2(2): 6103.901 (Likelihood ratio/deviance statistic) 
#X^2(2): 2.747466e+17 (Chi-square goodness of fit) 


#model <- poLCA(f, data=user_int, nclass = 2)

model = LCA_best_model

entropy<-function (p) sum(-p*log(p))

error_prior<-entropy(model$P)
error_post<-mean(apply(model$posterior,1, entropy),na.rm = TRUE)
round(((error_prior-error_post) / error_prior),3)
#entropy:0.957


#Calculate distribution difference with kl divergence

class_probs = model$probs

#kl_divergence <- function(p, q) {
#  sum(p * log(p / q))
#
# Define a function to calculate the KL divergence, adding Laplace smoothing
kl_divergence <- function(p, q, smooth=0.0001) {
  p_smoothed <- (p + smooth) / (1 + length(p) * smooth)
  q_smoothed <- (q + smooth) / (1 + length(q) * smooth)
  sum(p_smoothed * log(p_smoothed / q_smoothed))
}

kl_divergence_results <- c()

for (i in colnames(user_int)) {
  # Extract the probs of the current variable
   class1_probs <- class_probs[[i]][1]
  class2_probs <- class_probs[[i]][2]
  
  kl_divergence_result <- kl_divergence(class1_probs, class2_probs)
  
  kl_divergence_results <- c(kl_divergence_results, kl_divergence_result)
}

kl_divergence_results = abs(kl_divergence_results)
names(kl_divergence_results) = colnames(user_int)

#The ones with large differences in extraction distribution 
#indicate that they are characteristic in this class.
kl_divergence_results[kl_divergence_results>0.2]

plot(kl_divergence_results[kl_divergence_results>0.2])

# Get a list of unique variables in L2

lcmodel <- reshape2::melt(model$probs, level=2)
lcmodel = convert_to_original(lcmodel,convert_to_integer(user)$mapping)
unique_vars <- unique(lcmodel$L2)

# Create a list to store drawing objects
plots_list <- lapply(unique_vars, function(var) {
  subset_data <- subset(lcmodel, L2 == var)
  
  # Calculate the percentage for each category
  subset_data$percentage <- subset_data$value* 100
  
  # Create a bar graph
  p <- ggplot(subset_data, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(x = "Antwortkategorien", y = "Prozentsatz") +
    ggtitle(paste("Variable:", var)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5))
  filename <- paste("plot/", var, ".png", sep = "")
  ggsave(filename, p, width = 10, height = 6, dpi = 300)
  
  return(p)
})

for (plot in plots_list) {
  print(plot)
}

###########no user part#################

no_user = read.csv("./downloads/new_non_user.csv")

no_user_int = convert_to_integer(no_user)$data
no_user_int = no_user_int[,-c(1)]

f = as.matrix(no_user_int)~1

max_II <- -100000
min_bic <- 100000
for(i in 2:10){
  lc <- poLCA(f, no_user_int, nclass=2, maxiter=3000, 
              tol=1e-5, na.rm=FALSE,  
              nrep=10, verbose=TRUE, calc.se=TRUE)
  if(lc$bic < min_bic){
    min_bic <- lc$bic
    LCA_best_model<-lc
  }
}   



LCA_best_model

#AIC(2): 9694.716
#BIC(2): 10244.61
#G^2(2): 7156.97 (Likelihood ratio/deviance statistic) 
#X^2(2): 6.13879e+13 (Chi-square goodness of fit) 

model = LCA_best_model

entropy<-function (p) sum(-p*log(p))

error_prior<-entropy(model$P)
error_post<-mean(apply(model$posterior,1, entropy),na.rm = TRUE)
round(((error_prior-error_post) / error_prior),3)
#entropy:0.953

class_probs = model$probs

kl_divergence <- function(p, q, smooth=0.0001) {
  p_smoothed <- (p + smooth) / (1 + length(p) * smooth)
  q_smoothed <- (q + smooth) / (1 + length(q) * smooth)
  sum(p_smoothed * log(p_smoothed / q_smoothed))
}

kl_divergence_results <- c()

for (i in colnames(no_user_int)) {
  class1_probs <- class_probs[[i]][1]
  class2_probs <- class_probs[[i]][2]
  
  kl_divergence_result <- kl_divergence(class1_probs, class2_probs)
  
  kl_divergence_results <- c(kl_divergence_results, kl_divergence_result)
}

kl_divergence_results = abs(kl_divergence_results)
names(kl_divergence_results) = colnames(no_user_int)

kl_divergence_results[kl_divergence_results>0.2]

lcmodel <- reshape2::melt(model$probs, level=2)
lcmodel = convert_to_original(lcmodel,convert_to_integer(no_user)$mapping)
unique_vars <- unique(lcmodel$L2)

#plot
plots_list <- lapply(unique_vars, function(var) {
  subset_data <- subset(lcmodel, L2 == var)
  print(var)
  subset_data$percentage <- subset_data$value* 100
  
  p <- ggplot(subset_data, aes(x = Var2, y = percentage, fill = Var2)) +
    geom_bar(stat = "identity") +
    facet_grid(Var1 ~ .) +
    labs(x = "Antwortkategorien", y = "Prozentsatz") +
    ggtitle(paste("Variable:", var)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5))
  filename <- paste("plot/", var, ".png", sep = "")
  ggsave(filename, p, width = 10, height = 6, dpi = 300)
  
  return(p)
})

for (plot in plots_list) {
  print(plot)
}




