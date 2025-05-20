library(igraph)

set.seed(123)
N <- 324
density <- 0.029
num_edges <- round(density * N * (N - 1) / 2)

# Regular network
k <- round(2 * num_edges / N) + 1
regular_network <- sample_smallworld(dim = 1, size = N, nei = k/2, p = 0.0)

# Small-world network
p <- 0.1
synthetic_sm_network <- sample_smallworld(dim = 1, size = N, nei = k/2, p = p)

# Load real network

edge_data <- read.csv("Talaga-homophily-network/talaga-homophily-network-edges-N324_10.csv") #edge_data <- read.csv("Talaga-homophily-network/talaga-homophily-network-edges.csv")
real_sm_network <- graph_from_data_frame(edge_data, directed = FALSE, vertices = NULL)

# Scale-free network
# 256 m=4, 324 m=5, 400 m=6
scale_free_network <- barabasi.game(N, m = 5, directed = FALSE)

# Erdos-Renyi network
erdos_renyi_network <- erdos.renyi.game(N, num_edges, type = "gnm")

# Function to calculate network properties
calculate_properties <- function(network) {
  list(
    mean_degree = mean(degree(network)),
    sd_degree = sd(degree(network)),
    density = edge_density(network),
    clustering = transitivity(network, type = "global"),
    assortativity = assortativity_degree(network),
    avg_path_length = average.path.length(network)
  )
}

# Calculate properties for each network
properties_regular <- calculate_properties(regular_network)
properties_synthetic_sm <- calculate_properties(synthetic_sm_network)
properties_real_sm <- calculate_properties(real_sm_network)
properties_scale_free <- calculate_properties(scale_free_network)
properties_erdos_renyi <- calculate_properties(erdos_renyi_network)

# Mean degree
cat("Regular Network:", properties_regular$mean_degree, "\n")       # Regular:     12
cat("Small-world:", properties_synthetic_sm$mean_degree, "\n")      # Sintético:   12
cat("Real Small-world:", properties_real_sm$mean_degree, "\n")      # Real:        12.255
cat("Scale-free:", properties_scale_free$mean_degree, "\n")         # Scale-free:  11.895
cat("Erdos-Renyi:", properties_erdos_renyi$mean_degree, "\n\n")     # Erdos-Renyi: 11.57

# Sd degree
cat("Regular Network:", properties_regular$sd_degree, "\n")         # Regular:     0.0
cat("Small-world:", properties_synthetic_sm$sd_degree, "\n")        # Sintético:   1.639121 
cat("Real Small-world:", properties_real_sm$sd_degree, "\n")        # Real:        5.674546
cat("Scale-free:", properties_scale_free$sd_degree, "\n")           # Scale-free:  10.31601
cat("Erdos-Renyi:", properties_erdos_renyi$sd_degree, "\n\n")       # Erdos-Renyi: 3.237484

# Density
cat("Regular Network:", properties_regular$density, "\n")           # Regular:     0.03007519 
cat("Small-world:", properties_synthetic_sm$density, "\n")          # Sintético:   0.03007519 
cat("Real Small-world:", properties_real_sm$density, "\n")          # Real:        0.03071429 
cat("Scale-free:", properties_scale_free$density, "\n")             # Scale-free:  0.02981203 
cat("Erdos-Renyi:", properties_erdos_renyi$density, "\n\n")         # Erdos-Renyi: 0.02899749 

# Average Path Length
cat("Regular Network:", properties_regular$avg_path_length, "\n")   # Regular:     17.12782
cat("Small-world:", properties_synthetic_sm$avg_path_length, "\n")  # Sintético:   3.020852 
cat("Real Small-world:", properties_real_sm$avg_path_length, "\n")  # Real:        3.058208
cat("Scale-free:", properties_scale_free$avg_path_length, "\n")     # Scale-free:  2.583371
cat("Erdos-Renyi:", properties_erdos_renyi$avg_path_length, "\n\n") # Erdos-Renyi: 2.705113

# Clustering Coefficient
cat("Regular Network:", properties_regular$clustering, "\n")        # Regular:     0.6818182 
cat("Small-world:", properties_synthetic_sm$clustering, "\n")       # Sintético:   0.3228764 
cat("Real Small-world:", properties_real_sm$clustering, "\n")       # Real:        0.226257
cat("Scale-free:", properties_scale_free$clustering, "\n")          # Scale-free:  0.0766702 
cat("Erdos-Renyi:", properties_erdos_renyi$clustering, "\n\n")      # Erdos-Renyi: 0.02903955 

# Assortativity Coefficient
cat("Regular Network:", properties_regular$assortativity, "\n")     # Regular:     NaN
cat("Small-world:", properties_synthetic_sm$assortativity, "\n")    # Sintético:  -0.01604385
cat("Real Small-world:", properties_real_sm$assortativity, "\n")    # Real:        0.2331087
cat("Scale-free:", properties_scale_free$assortativity, "\n")       # Scale-free: -0.04929535
cat("Erdos-Renyi:", properties_erdos_renyi$assortativity, "\n\n")   # Erdos-Renyi:-0.003567295