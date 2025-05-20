#Topological Measures for Identifying and Predicting the Spread of Complex Contagions
#Guilbeault & Centola, 2021, Nature Communications 
#July 2021 

#Overview: 
#This script calculates average complex path length for all nodes in a given graph, g
#which provides a measure of complex centrality; methods w/wo parallel processing are provided. 
#A pipeline is provided for replicating the main analyses from Guilbeault & Centola (2021)
#using both simulated graphs and empirical Addhealth networks (the latter requires downloading
#the addhealth network topologies from the following github: https://github.com/drguilbe/complexpaths). 

#Note: 
#For the purposes of algorithmic efficiency, this code makes a simplifying assumption in the calculation
#of complex path length - i.e., the complex path length between seed node i and target node j is calculated on the 
#final subgraph of activated nodes based on seeding with N[i], rather than at every time step in the diffusion process, 
#thereby significantly reducing the number of calculations required. This introduces a small amount of noise when 
#estimating complex path length, but this does not qualitatively or significantly alter the ability for this code 
#to replicate our main results when comparing seeding strategies. The full code base with both versions of the procedure
#will be made available in future releases. 

# El unico cambio que hice fue length(seeds) por length(unique(seeds)) 

#Load Libraries
#rm(list=ls());gc() # clean environment
library(dplyr)
library(tidyr)
library(influential)
library(fastnet)
library(igraph)
library(doParallel)
library(parallel)

####### Model Functions -------------------------------------------

# Variable between 0 and 1
min_max_norm<-function(x){(x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))}

# Select the neighborhoods of a seed.
clustered_seeding <- function(seeds, g, seeds_needed) {
  print("Using clustered_seeding()")
  
  # Find the neighbors of each seed in the graph and combine them into a single vector
  possible_seeds <- unique(unlist(sapply(seeds, function(x) neighbors(g, x, mode = "total"))))
  
  # Filter the neighbors that are not already in the seed set
  seeds_to_add <- possible_seeds[!possible_seeds %in% seeds]
  
  # Check if more seeds are needed than can be added
  need_seeds <- seeds_needed > length(seeds_to_add)
  
  # If more seeds are needed than available, return all possible neighbors
  if (need_seeds) {
    return(possible_seeds)
  } else {
    # If not, randomly select the needed seeds and combine them with the initial seeds
    final_seeds <- c(seeds, sample(seeds_to_add, seeds_needed))
    return(final_seeds)
  }
}

# Data frame with centralities
get_simple_centralities <- function(g) {
  
  # Create a data frame with node centralities
  centrality_df <- data.frame(
    seed = as.numeric(V(g)), 
    degree = as.numeric(degree(g)), 
    betweenness = as.numeric(betweenness(g)), 
    eigen = as.numeric(eigen_centrality(g)$vector)
  )
  
  # Calculate the collective influence of the nodes
  # https://cran.r-project.org/web/packages/influential/vignettes/Vignettes.html#collective-influence
  percolation <- collective.influence(graph = g, vertices = V(g), mode = "all", d = 3) #3 steps of influence
  
  # Add the collective influence to the centrality data frame
  centrality_df$percolation <- as.numeric(percolation)
  
  # Return the data frame with all centralities
  return(centrality_df)
}

# Simulation of complex spread
get_complex <- function(seed, N, g, gmat, thresholds, num_seeds_to_add, model_output_list) {
  # Initialize a simulation matrix with zeros
  gmat_simulation <- matrix(nrow = N, ncol = N, 0)
  
  # Get the number of seeds to add for the given seed
  num_seeds_i <- num_seeds_to_add[seed]
  seeds_to_add <- numeric(num_seeds_i)
  
  if (num_seeds_i > 0) {
    # Get the neighbors of the seed in the graph
    seeds <- as.numeric(neighbors(g, seed, mode = "total"))
    
    # Calculate wich nodes will be seeds
    num_seeds_needed <- num_seeds_i - length(unique(seeds)) # x(-1) ¿? (esto será más como los nodos que me sobran de mi vecindario)
    need_seeds <- num_seeds_needed > 0
    
    if (need_seeds) { # Initiate clustered seeding
      seeds <- clustered_seeding(seeds, g, num_seeds_needed)
      num_seeds_needed <- num_seeds_i - length(unique(seeds))
      need_seeds <- num_seeds_needed > 0
    } else if (length(seeds) > 1) {
      seeds <- c(seed, sample(seeds, num_seeds_i))
    } else {
      seeds <- c(seed, seeds)
    }
  }
  
  # Initialize the activated vector and update the simulation matrix
  activated <- logical(N)
  activated[seeds] <- TRUE
  gmat_simulation[seeds, ] <- 1
  gmat_simulation[, seeds] <- 1
  gmat_simulation_run <- gmat * gmat_simulation
  
  spread = TRUE
  while (spread) {
    # Calculate the influence of each node
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    
    # Update the activated nodes
    t <- which(activated)
    t_1 <- which(influence_activated)
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t) # there is a moment when spread stops
    activated[t_full] <- TRUE
    
    # Update the simulation matrix for the new adopters
    adopters <- which(activated)
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
  }
  
  # Once stopped the spread, relevant quantities are calculated
  
  # Calculate the number of adopters
  num_adopters <- sum(activated)
  
  # Create a sub-graph from the adjacency matrix of the simulation
  complex_g <- graph_from_adjacency_matrix(gmat_simulation_run)
  
  # Calculate the distances from the seed to all other nodes
  all_distances <- distances(complex_g, seed, V(complex_g), mode = "out")
  all_distances[all_distances == "Inf"] <- 0
  
  # Calculate the average path length
  PLci <- mean(all_distances)
  
  # Update the model output list
  model_output_list <- c(seed, N, num_seeds_i, num_adopters, PLci)
  
  # Return the model output list
  return(model_output_list)
}


#####Replicate Results With Parallel Methods ### -----------------------------

# Detect the number of cores available for parallel processing
detectCores(logical = TRUE)

# Create a cluster with 8 cores
cluster <- makeCluster(8)

# Register the cluster for parallel processing
registerDoParallel(cluster)

# Export necessary functions and objects to the cluster with 8 cores
clusterExport(cluster, list('neighbors', 'graph_from_adjacency_matrix', 'distances', 'V', 'clustered_seeding'))

# Replicate with SIMULATED graphs
num_graphs <- 20
N <- 200
# Initialize a list to store the graphs
graphs <- c()

# Generate 20 simulated graphs using the Holme & Kim Scale-Free algorithm
for (i in 1:num_graphs) {
  #graphs[[i]] <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)
  graphs[[i]] <- graph
}

# Initialize a list to store the adjacency matrices of the graphs
gmats <- c()
gmats[[1]] <- as.matrix(as_adjacency_matrix(graph))

# Convert each graph to its adjacency matrix and store it in the list
for (i in 1:num_graphs) {
  gmats[[i]] <- as.matrix(as_adjacency_matrix(graphs[[i]]))
}

# Replicate with EMPIRICAL Addhealth networks 
# addhealth_mats <- list.files("C:/Users/dougl/Desktop/AddHealth_Networks_Largest_Components/") # Download zipped folder from github and set directory 
# gmats <- c(); for(i in 1:length(addhealth_mats)){gmats[[i]] <- read.csv(paste("C:/Users/dougl/Desktop/AddHealth_Networks_Largest_Components/", addhealth_mats[i], sep=""))}
# for(i in 1:length(gmats)){colnames(gmats[[i]]) <- 1:length(gmats[[i]]); rownames(gmats[[i]]) <- 1:length(gmats[[i]])}  
# for(i in 1:length(gmats)){gmats[[i]] <- as.matrix(gmats[[i]])}  
# graphs <- c(); for(i in 1:length(addhealth_mats)){graphs[[i]] <- graph_from_adjacency_matrix(gmats[[i]])} 

# Set threshold type
T_dist <- c("homo", "hetero")[1]  # Set the threshold distribution type to "homo" (homogeneous)
T_type <- c("abs", "frac")[1]  # Set the threshold type to "abs" (absolute)

# Initialize an empty data frame to store comparison results
comparison_df <- data.frame()

# Loop over each graph in the list of graphs
for(i in 1:length(graphs)) { 
  print(paste("Graph: ", i, sep=""))  # Print the current graph index
  g <- graphs[[i]]  # Get the current graph
  gmat <- gmats[[i]]  # Get the adjacency matrix of the current graph
  N <- nrow(gmat)  # Get the number of nodes in the graph
  
  # Check if the graph has more than 10 nodes to avoid errors with empty networks
  if(length(V(g)) > 10) { 
    
    # Set thresholds for each node (default to 3)
    thresholds <- replicate(N, 3) 
    #thresholds <- replicate(N, runif(1, 0.1,0.5)) #for Fractional Thresholds
    # If the threshold type is fractional, adjust the thresholds based on the number of neighbors
    if(T_type == "frac") {
      thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(g, x, mode = "total")))))
    }
    
    # Calculate the number of seeds to add for each node
    num_seeds_to_add <- thresholds - 1
    num_seeds_to_add[num_seeds_to_add < 0] = 0  # Ensure no negative values
    thresholds[thresholds <= 0] = 1  # Ensure thresholds are at least 1
    
    # Get simple centralities for the current graph
    simple_centralities_df <- get_simple_centralities(g)
    
    # Get complex centrality 
    # Initialize a list to store model output for each node
    model_output_list <- vector(mode = "list", length = N)
    start_time <- proc.time()  # Record the start time of the model run
    
    # Run the model in parallel for each node. Return: (seed, N, num_seeds_i, num_adopters, PLci)
    model_output_list <- clusterMap(cluster, get_complex, seed = 1:N, MoreArgs = list(N = N, g = g, gmat = gmat, thresholds = thresholds, num_seeds_to_add = num_seeds_to_add, model_output_list = model_output_list)) 
    print(proc.time() - start_time)  # Print the runtime of the model
    
    # Combine the model output into a data frame
    model_output_list_df <- as.data.frame(do.call('rbind', model_output_list))
    colnames(model_output_list_df) <- c("seed", "N", "num_neigh_seeds", "num_adopters", "PLci")
    
    # Normalize the path length centrality index (PLci)
    model_output_list_df$PLci_norm <- min_max_norm(model_output_list_df$PLci)
    #### ----
    
    # Merge the model output with the simple centralities data frame
    model_output_parallel_full <- merge(model_output_list_df, simple_centralities_df, by = "seed") 
    model_output_parallel_full$threshold <- thresholds  # Add thresholds to the data frame
    model_output_parallel_full$T_dist <- T_dist  # Add threshold distribution type to the data frame
    model_output_parallel_full$T_type <- T_type  # Add threshold type to the data frame
    
    # Select the best node 
    # normalized PLci
    top_plc <- model_output_parallel_full %>% top_n(1, wt = PLci_norm) %>% mutate(top = "Complex") %>% sample_n(1)
    # betweenness centrality
    top_btwn <- model_output_parallel_full %>% top_n(1, wt = betweenness) %>% mutate(top = "betweenness") %>% sample_n(1)
    # degree centrality
    top_degree <- model_output_parallel_full %>% top_n(1, wt = degree) %>% mutate(top = "degree") %>% sample_n(1)
    # eigenvector centrality
    top_eigen <- model_output_parallel_full %>% top_n(1, wt = eigen) %>% mutate(top = "eigen") %>% sample_n(1)
    # percolation centrality
    top_percolation <- model_output_parallel_full %>% top_n(1, wt = percolation) %>% mutate(top = "percolation") %>% sample_n(1)
    
    # Combine the top nodes into a single data frame
    top_final <- rbind(top_plc, top_btwn, top_degree, top_eigen, top_percolation)
    top_final$graph <- i  # Add the graph index to the data frame
    
    # Append the results to the comparison data frame
    comparison_df <- rbind(comparison_df, top_final)
  }
}

stopCluster(cluster)

# Perform pairwise Wilcoxon tests to compare the number of adopters between different centrality measures
pairwise.wilcox.test(comparison_df$num_adopters, comparison_df$top, p.adjust.method = "none", paired = TRUE)

# Calculate the mean number of adopters for
# "Complex" centrality 
mean(subset(comparison_df, top == "Complex")$num_adopters)
# "betweenness" centrality 
mean(subset(comparison_df, top == "betweenness")$num_adopters)
#"degree" centrality 
mean(subset(comparison_df, top == "degree")$num_adopters)
# "eigen" centrality 
mean(subset(comparison_df, top == "eigen")$num_adopters)
# "percolation" centrality 
mean(subset(comparison_df, top == "percolation")$num_adopters)

#############
## Model without Parallel Processing ## --------------------------------------

N<-200
g<-net.holme.kim(N,4,0.5)
g<-to.igraph(g)
gmat<-as.matrix(as_adjacency_matrix(g))
T_dist<-c("homo", "hetero")[1]
T_type<-c("abs","frac")[1]
thresholds<-replicate(N, 3)
#thresholds<-replicate(N, runif(1, 0.1,0.5)) #fractional Ts distributed uniformly at random. 
if(T_type == "frac"){thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(g, x, mode = "total")))))}
num_seeds_to_add<-thresholds-1
num_seeds_to_add[num_seeds_to_add<0]=0
thresholds[thresholds<=0]=1
simple_centralities_df<-get_simple_centralities(g)

#Run Model
model_output_list <- vector(mode = "list", length = N)
start_time <- proc.time()
model_output_list<-lapply(1:N, function(x) get_complex(x, N, g, gmat, thresholds, num_seeds_to_add, model_output_list)) #Run model; change 1:N to any subset of specific seeds to narrow search
print(proc.time() - start_time) #view runtime 
model_output_df <- as.data.frame(Reduce(rbind, model_output_list))
colnames(model_output_df)<-c("seed","N","num_neigh_seeds", "num_adopters","PLci")
model_output_df$PLci_norm<-min_max_norm(model_output_df$PLci)
model_output_full<-merge(model_output_df, simple_centralities_df, by="seed") #map to simple centralities of seed
