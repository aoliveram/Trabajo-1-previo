library(igraph)
library(doParallel)
library(animation)
library(magick)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(parallel)
library(influential)

plot <- FALSE

########################## Graph generation ----------------------------------------
rows <- 16
cols <- 16
N <- rows * cols

######## Small-world-ish

set.seed(123)
graph <- make_lattice(dimvector = c(rows, cols), nei = 1, directed = FALSE)
add_diagonal_edges <- function(graph, rows, cols) {
  for (i in 1:rows) {
    for (j in 1:cols) {
      node <- (i - 1) * cols + j
      if (i < rows && j < cols) {
        graph <- add_edges(graph, c(node, node + cols + 1))
      }
      if (i < rows && j > 1) {
        graph <- add_edges(graph, c(node, node + cols - 1))
      }
    }
  }
  return(graph)
}
graph <- add_diagonal_edges(graph, rows, cols)

rewire_edges <- function(graph, p) {
  edges <- as_edgelist(graph)
  num_edges <- nrow(edges)
  num_rewire <- round(p * num_edges)
  for (i in 1:num_rewire) {
    edge_to_remove <- sample(1:num_edges, 1)
    graph <- delete_edges(graph, E(graph)[edge_to_remove])
    repeat {
      new_edge <- sample(V(graph), 2)
      if (!are_adjacent(graph, new_edge[1], new_edge[2])) {
        graph <- add_edges(graph, new_edge)
        break
      }
    }
  }
  return(graph)
}
graph <- rewire_edges(graph, 0.1)

######### Small-world (0.2)

set.seed(123)
density <- 0.029
num_edges <- round(density * N * (N - 1) / 2)
k <- round(2 * num_edges / N) + 1 # Number of neighbors
p <- 0.2
graph <- sample_smallworld(dim = 1, size = N, nei = k/2, p = p)

######### Scale-free

#set.seed(123)
#graph <- barabasi.game(N, m = 4, directed = FALSE)

######### Erdos-Renyi

set.seed(123)
density <- 0.029
num_edges <- round(density * N * (N - 1) / 2)
graph <- erdos.renyi.game(N, num_edges, type = "gnm")

######### Plot

print(edge_density(graph))

if (plot) {
  layout <- layout_on_grid(graph)
  plot(graph, layout = layout, vertex.size = 3, edge.width = 0.5, vertex.label.cex = 0.5, 
       edge.color = "gray", edge.arrow.size = 0)
}

########################## Functions -----------------------------------------

State <- function(value = "State1") {
  if (!value %in% c("State1", "State2")) {
    stop("Estado no válido. Use 'State1' o 'State2'.")
  }
  return(value)
}

clustered_seeding <- function(seeds, g, num_seeds_needed) {
  new_seeds <- sample(setdiff(V(g), seeds), num_seeds_needed)
  return(c(seeds, new_seeds))
}

plot_network <- function(g, activated, N, layout, step, file_dir, alpha_1, homoph) {
  vertex_color <- ifelse(activated == "State2", "red", "grey")
  
  edge_color <- rep("gray", ecount(g))
  E(g)$color <- edge_color
  
  for (i in 1:N) {
    if (activated[i] == "State2") {
      neighbors_i <- neighbors(g, i, mode = "total")
      for (j in neighbors_i) {
        E(g, P = c(i, j))$color <- "red"
      }
    }
  }
  
  png(paste0(file_dir, "/step_", sprintf("%03d", step), ".png"), width = 1200, height = 900)
  plot(g, layout = layout, vertex.color = vertex_color,
       vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.arrow.size = 0, 
       main = sprintf("Homophily:  %.2f | Step: %d", homoph, step))
  dev.off()
}

get_complex_plot <- function(seed, N, g, thresholds, num_seeds_to_add, plot, file_dir, alpha_1, homoph) {
  
  set.seed(123)
  
  # Selección de semillas adicionales para la infección
  seeds <- as.numeric(neighbors(g, seed, mode = "total"))
  num_seeds_i <- num_seeds_to_add[seed]
  num_seeds_needed <- num_seeds_i - length(seeds)
  
  if (num_seeds_needed > 0) {
    seeds <- clustered_seeding(seeds, g, num_seeds_needed)
  } else if (length(seeds) > 1) {
    seeds <- c(seed, sample(seeds, num_seeds_i))
  } else {
    seeds <- c(seed, seeds)
  }
  
  # Inicializar estados
  activated <- rep("State1", N)
  activated[seeds] <- "State2"
  
  gmat_simulation <- matrix(0, nrow = N, ncol = N)
  gmat_simulation[seeds, ] <- 1
  gmat_simulation[, seeds] <- 1
  gmat <- as.matrix(as_adjacency_matrix(g))
  gmat_simulation_run <- gmat * gmat_simulation
  spread <- TRUE
  step <- 1
  
  # Fijar el layout de la red
  layout <- layout_on_grid(g)
  
  while (spread) {
    #print(step)
    if (plot) {
      plot_network(g, activated, N, layout, step, file_dir, alpha_1, homoph)
    }
    
    # Determinar nuevos infectados
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    t <- which(activated == "State2")
    t_1 <- which(influence_activated & activated == "State1" & V(g)$alpha < alpha_1)
    
    # Nueva regla: permitir infección si alpha es similar al de los vecinos
    t_2 <- which(influence_activated & activated == "State1" & V(g)$alpha >= alpha_1)
    for (i in t_2) {
      vecinos <- neighbors(g, i)
      similar_alpha_neighbors <- sum(abs(V(g)$alpha[i] - V(g)$alpha[vecinos]) < homoph)
      if (similar_alpha_neighbors >= thresholds[i]) {
        t_1 <- c(t_1, i)
      }
    }
    
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t)
    
    activated[t_full] <- "State2"
    
    adopters <- which(activated == "State2")
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
    
    all_infected <- all(activated != "State1")
    if (all_infected) {
      spread <- FALSE
    }
    
    step <- step + 1
  }
  
  if (plot) {
    plot_network(g, activated, N, layout, step, file_dir, alpha_1, homoph)
  }
  
  num_adopters <- sum(activated == "State2")
  return(list(alpha_1 = alpha_1 ,homophily = homoph, seed = seed, num_adopters = num_adopters, num_steps = step))
}

sweep_homoph_parameter <- function(seed, homoph_values, N, graph, thresholds, num_seeds_to_add, plot, file_dir, alpha_values, graph_file, centrality) {
  
  # Define the path of the main folder
  by_homoph <- diff(homoph_values)[1]
  num_decimals_homoph <- nchar(strsplit(sub('0+$', '', as.character(by_homoph)), ".", fixed = TRUE)[[1]][[2]])
  format_homoph <- sprintf("%%.%df", num_decimals_homoph)
  homoph_min <- sprintf(format_homoph, min(homoph_values))
  homoph_max <- sprintf(format_homoph, max(homoph_values))
  
  by_alpha <- diff(alpha_values)[1]
  num_decimals_alpha <- nchar(strsplit(sub('0+$', '', as.character(by_alpha)), ".", fixed = TRUE)[[1]][[2]])
  format_alpha <- sprintf("%%.%df", num_decimals_alpha)
  alpha_min <- sprintf(format_alpha, min(alpha_values))
  alpha_max <- sprintf(format_alpha, max(alpha_values))
  
  parameter_sweep_homoph_folder <- sprintf("Simulation_Stats_Simple/%s/1inf_hom_%s_%s_by_%s_alpha_%s_%s_by_%s", 
                                           graph_file, homoph_min, homoph_max, by_homoph, alpha_min, alpha_max, by_alpha)
  
  # Create the main folder if it does not exist
  if (!dir.exists(parameter_sweep_homoph_folder)) dir.create(parameter_sweep_homoph_folder, recursive = TRUE)
  
  # Initialize an empty list to store data frames
  list_of_dfs <- list()
  
  # Create a cluster
  cluster <- makeForkCluster(10)
  
  # Define the function to be run in parallel
  run_simulation <- function(alpha_1) {
    sub_list_of_dfs <- list()
    for (j in seq_along(homoph_values)) {
      homoph <- homoph_values[j]
      results <- get_complex_plot(seed, N, graph, thresholds, num_seeds_to_add, plot, file_dir, alpha_1, homoph)
      results_df <- as.data.frame(results)
      sub_list_of_dfs <- append(sub_list_of_dfs, list(results_df))
      
      if (plot) {
        images <- list()
        for (k in 1:results$num_steps) {
          img <- image_read(sprintf("%s/step_%03d.png", file_dir, k))
          images <- c(images, img)
        }
        image_write(image_animate(image_join(images), fps = 1), sprintf("%s/1inf_%02d.gif", parameter_sweep_homoph_folder, which(alpha_1 == alpha_values)))
      }
    }
    return(sub_list_of_dfs)
  }
  
  # Use clusterMap to run the simulations in parallel
  results_list <- clusterMap(cluster, run_simulation, alpha_values)
  
  # Stop the cluster
  stopCluster(cluster)
  
  # Combine all data frames into one and assign to a dynamically named variable
  for (res in results_list) {
    list_of_dfs <- append(list_of_dfs, res)
  }
  
  combined_results_df <- do.call(rbind, list_of_dfs)
  combined_results_name <- sprintf("results_%s_%s", homoph_min, homoph_max)
  assign(combined_results_name, combined_results_df, envir = .GlobalEnv)
  
  # Save the combined data frame in the main folder
  write.csv(get(combined_results_name), sprintf("%s/results_1inf_hom_%s_%s_by_%s_alpha_%s_%s_by_%s_%s.csv", parameter_sweep_homoph_folder, homoph_min, homoph_max, by_homoph, alpha_min, alpha_max, by_alpha, centrality), row.names = FALSE)
  
  return("Done!")
}

# Función para ploteo de la red en una infección simple
plot_network_single <- function(g, activated, N, layout, step, file_dir, homoph) {
  vertex_color <- ifelse(activated, "red", "grey")
  
  edge_color <- rep("gray", ecount(g))
  E(g)$color <- edge_color
  
  for (i in 1:N) {
    if (activated[i]) {
      neighbors_i <- neighbors(g, i, mode = "total")
      for (j in neighbors_i) {
        E(g, P = c(i, j))$color <- "red"
      }
    }
  }
  
  png(paste0(file_dir, "/step_", sprintf("%03d", step), ".png"), width = 1200, height = 900)
  plot(g, layout = layout, vertex.color = vertex_color,
       vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.arrow.size = 0, 
       main = sprintf("Homophily:  %.2f | Step: %d", homoph, step))
  dev.off()
}

# To get centralities

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
  percolation <- collective.influence(graph = g, vertices = V(g), mode = "all", d = 3) # 3 steps of influence
  
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
    
    # Calculate which nodes will be seeds
    num_seeds_needed <- num_seeds_i - length(unique(seeds))
    need_seeds <- num_seeds_needed > 0
    
    if (need_seeds) {
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
  
  spread <- TRUE
  while (spread) {
    # Calculate the influence of each node
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    
    # Update the activated nodes
    t <- which(activated)
    t_1 <- which(influence_activated)
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t)
    activated[t_full] <- TRUE
    
    # Update the simulation matrix for the new adopters
    adopters <- which(activated)
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
  }
  
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

calculate_top_nodes <- function(graph, T_dist, T_type) {
  N <- vcount(graph)
  gmat <- as.matrix(as_adjacency_matrix(graph))  # Get the adjacency matrix of the current graph
  
  thresholds <- replicate(N, 2) 
  # If the threshold type is fractional, adjust the thresholds based on the number of neighbors
  if (T_type == "frac") {
    thresholds <- round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(graph, x, mode = "total")))))
  }
  
  # Calculate the number of seeds to add for each node
  num_seeds_to_add <- thresholds - 1
  num_seeds_to_add[num_seeds_to_add < 0] <- 0  # Ensure no negative values
  thresholds[thresholds <= 0] <- 1  # Ensure thresholds are at least 1
  
  # Get simple centralities for the current graph
  simple_centralities_df <- get_simple_centralities(graph)
  
  # Run Model
  model_output_list <- vector(mode = "list", length = N)
  model_output_list <- lapply(1:N, function(x) get_complex(x, N, graph, gmat, thresholds, num_seeds_to_add, model_output_list))
  model_output_df <- as.data.frame(Reduce(rbind, model_output_list))
  colnames(model_output_df) <- c("seed", "N", "num_neigh_seeds", "num_adopters", "PLci")
  model_output_full <- merge(model_output_df, simple_centralities_df, by = "seed")
  
  return(model_output_full)
}

########################## Complex simulation (2D) ----------------

file_dir <- "Simulation_Stats_Simple/plot_data"
if (!dir.exists(file_dir)) {
  dir.create(file_dir)
}

# Crear un valor de alpha arbitrario entre 0 y 1 para cada nodo
set.seed(123)  # Para reproducibilidad
V(graph)$alpha <- runif(N, 0, 1)

# Definir umbrales para la infección
thresholds <- replicate(N, 1)
num_seeds_to_add <- thresholds - 1
num_seeds_to_add[num_seeds_to_add < 0] <- 0
thresholds[thresholds <= 0] <- 1

# Seleccionar semillas basadas en el valor de alpha
alpha_1 <- 0.5
seed_candidates <- V(graph)[alpha < alpha_1]

set.seed(123)
seed <- as.integer(sample(seed_candidates, 1))

# Configurar clúster para computación en paralelo
cores <- detectCores() - 1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
clusterExport(cluster, list('plot_network', 'neighbors', 'graph_from_adjacency_matrix',
                            'distances', 'V', 'clustered_seeding', 'get_complex_plot', 'sweep_homoph_parameter'))

homoph_values <- seq(0.0, 0.5, by = 0.005)

plot <- FALSE

# Ejecutar el modelo
results <- sweep_homoph_parameter(seed, homoph_values, N, graph, thresholds, num_seeds_to_add, plot, file_dir, alpha_1)

# Finalizar clúster
stopCluster(cluster)

data_homophily_sweep <- read.csv("Simulation_Stats_Simple/homoph_sweep_single_00_03/results_00_03.csv")

ggplot(data_homophily_sweep, aes(x = homophily, y = num_adopters)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(title = "Relationship between Homophily and Number of Adopters",
       x = "Level of Homophily",
       y = "Number of Adopters") +
  theme_minimal()

########################## Complex simulation (3D) ---------------------------

graph_file <- 'regular-graph-ish'

# Crear un valor de alpha arbitrario entre 0 y 1 para cada nodo
set.seed(123)  # Para reproducibilidad
V(graph)$alpha <- runif(N, 0, 1)

# Definir umbrales para la infección
thresholds <- replicate(N, 1)
num_seeds_to_add <- thresholds - 1
num_seeds_to_add[num_seeds_to_add < 0] <- 0
thresholds[thresholds <= 0] <- 1

# Set threshold type
T_dist <- c("homo", "hetero")[1]
T_type <- c("abs", "frac")[1]
result <- calculate_top_nodes(graph, T_dist, T_type)

# Top nodes centralities
top_PLci <- order(result$PLci, decreasing = TRUE)[1:10]
top_degree <- order(result$degree, decreasing = TRUE)[1:10]
top_betweenness <- order(result$betweenness, decreasing = TRUE)[1:10]
top_eigen <- order(result$eigen, decreasing = TRUE)[1:10]
top_percolation <- order(result$percolation, decreasing = TRUE)[1:10]

centralities <- list('PLci', 'degree', 'betweenness', 'eigen', 'percolation')
top_nodes_list <- list(top_PLci, top_degree, top_betweenness, top_eigen, top_percolation)

homoph_values <- seq(0.0, 0.5, by = 0.01)
alpha_values <- seq(0.0, 1.0, by = 0.01)

plot <- FALSE
centralities <- list('PLci', 'degree', 'betweenness', 'eigen', 'percolation')
execution_times <- list()

for (i in 1:length(top_nodes_list)) {
  seed <- top_nodes_list[[i]][1]
  centrality <- centralities[[i]]
  
  start_time <- Sys.time()
  sweep_results <- sweep_homoph_parameter(seed, homoph_values, N, graph, thresholds,
                                          num_seeds_to_add, plot, file_dir, alpha_values, graph_file, centrality)
  end_time <- Sys.time()
  time <- end_time - start_time
  print(time)
  execution_times[[i]] <- time
}

########################## Plot 3D -------------------------------------------

centralities <- list('PLci', 'degree', 'betweenness', 'eigen', 'percolation')

data_dir <- sprintf("Simulation_Stats_Simple/%s/1inf_hom_0.00_0.50_by_0.01_alpha_0.00_1.00_by_0.01", graph_file)
file_list <- list.files(path = data_dir, pattern = "*.csv", full.names = TRUE)

for (i in 1:length(centralities)) {
  
  centrality <- centralities[[i]]
  file_path <- file_list[grep(centrality, file_list)]
  data_sims <- read.csv(file_path)
  num_adopters_matrix_data <- matrix(data_sims$num_adopters, nrow = length(unique(data_sims$alpha_1)), byrow = TRUE)
  
  fig <- plot_ly(
    x = unique(data_sims$homophily),
    y = unique(data_sims$alpha_1),
    z = num_adopters_matrix_data,
    type = 'surface'
  ) %>%
    layout(
      title = sprintf("Relationship between Homophily, Alpha_1, and Number of Adopters (%s)", centrality),
      scene = list(
        xaxis = list(title = 'Hom. Window'),
        yaxis = list(title = 'Alpha'),
        zaxis = list(title = 'Number of Adopters')
      )
    )
  
  output_file <- sprintf("Simulation_Stats_Simple/%s/1inf_hom_0.00_0.50_by_0.01_alpha_0.00_1.00_by_0.01_%s.html", graph_file, centrality)
  saveWidget(fig, output_file)
}
