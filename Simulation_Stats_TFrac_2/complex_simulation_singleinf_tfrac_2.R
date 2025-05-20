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
real_sm_net_homoph <- TRUE

########################## Graph generation ----------------------------------------

######### Small-world SDA

if (real_sm_net_homoph) {
  networks_dir <- "Talaga-homophily-network/"  # Cambia esto a la ruta correcta
  real_small_world_networks <- list()
  attributes <- list()
  
  for (i in 1:5) {
    edge_file <- paste0(networks_dir, "talaga-homophily-network-edges-N324_", i, ".csv")
    edges <- read.csv(edge_file)
    attribute_file <- paste0(networks_dir, "talaga-homophily-network-node-attributes-N324_", i, ".csv")
    node_attributes <- read.csv(attribute_file)
    
    g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = node_attributes)
    
    relevant_dim <- node_attributes$relevant_dim
    
    V(g)$alpha <- relevant_dim
    
    graphs[[i]] <- g
    attributes[[i]] <- node_attributes
  }
}
#graph <- graphs[[1]]
#V(graph)$alpha

######### Small-world (0.2)
#mean(degree(graph))256 =8, mean(degree(graph))400 =12
#density256 =0.0313, density256 =0.030
generate_small_world_networks <- function(N, num_edges, p) {
  small_world_graphs <- list()
  for (i in 1:5) {
    set.seed(i)  # Semilla diferente para cada iteración
    k <- round(2 * num_edges / N) + 1
    graph <- sample_smallworld(dim = 1, size = N, nei = k / 2, p = p)
    V(graph)$alpha <- runif(N, 0, 1)
    small_world_graphs[[i]] <- graph
  }
  return(small_world_graphs)
}

######### Scale-free
#mean(degree(graph))256 =7.92 [m=4], mean(degree(graph))400 =11.895 [m=6]
#density256 =0.031 [m=4], density400 =0.0298 [m=6]
generate_scale_free_networks <- function(N, m) {
  scale_free_graphs <- list()
  for (i in 1:5) {
    set.seed(i)  # Semilla diferente para cada iteración
    graph <- barabasi.game(N, m = m, directed = FALSE)
    V(graph)$alpha <- runif(N, 0, 1)
    scale_free_graphs[[i]] <- graph
  }
  return(scale_free_graphs)
}

######### Erdos-Renyi
#mean(degree(graph))256 =7.39, mean(degree(graph))400 =11.57
#density400 =0.0290, density400 =0.02899}
generate_erdos_renyi_networks <- function(N, num_edges) {
  erdos_renyi_graphs <- list()
  for (i in 1:5) {
    set.seed(i)  # Semilla diferente para cada iteración
    graph <- erdos.renyi.game(N, num_edges, type = "gnm")
    V(graph)$alpha <- runif(N, 0, 1)
    erdos_renyi_graphs[[i]] <- graph
  }
  return(erdos_renyi_graphs)
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

sweep_homoph_parameter <- function(seed, homoph_values, N, graph, thresholds, num_seeds_to_add, plot, file_dir, alpha_values, graph_file, centrality, T_str, i) {
  
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
  
  parameter_sweep_homoph_folder <- sprintf("Simulation_Stats_TFrac_2/%s/1inf_N%s_hom_%s_%s_by_%s_alpha_%s_%s_by_%s_T%s",  
                                           graph_file, N, homoph_min, homoph_max, by_homoph, alpha_min, alpha_max, by_alpha, T_str)
  
  # Create the main folder if it does not exist
  if (!dir.exists(parameter_sweep_homoph_folder)) dir.create(parameter_sweep_homoph_folder, recursive = TRUE)
  
  # Initialize an empty list to store data frames
  list_of_dfs <- list()
  
  # Create a cluster
  cluster <- makeForkCluster(24)
  
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
          img <- image_read(sprintf("Simulation_Stats_TFrac_2/%s/step_%03d.png", file_dir, k))
          images <- c(images, img)
        }
        image_write(image_animate(image_join(images), fps = 1), sprintf("Simulation_Stats_TFrac_2/%s/1inf_%02d.gif", parameter_sweep_homoph_folder, which(alpha_1 == alpha_values)))
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
  write.csv(get(combined_results_name), sprintf("%s/1inf_hom_%s_%s_by_%s_alpha_%s_%s_by_%s_T%s_graph_%s_seed%s.csv", parameter_sweep_homoph_folder, homoph_min, homoph_max, by_homoph, alpha_min, alpha_max, by_alpha, T_str, i, seed), row.names = FALSE)
  
  return(combined_results_df)
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

calculate_top_nodes <- function(graph, T_dist, T_type, thresholds) {
  N <- vcount(graph)
  gmat <- as.matrix(as_adjacency_matrix(graph))  # Get the adjacency matrix of the current graph
  # If the threshold type is fractional, adjust the thresholds based on the number of neighbors
  #if (T_type == "frac") {
  #  thresholds_local <- round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(graph, x, mode = "total")))))
  #}
  # Calculate the number of seeds to add for each node
  #num_seeds_to_add <- thresholds_local - 1
  #num_seeds_to_add[num_seeds_to_add < 0] <- 0  # Ensure no negative values
  #thresholds_local[thresholds_local <= 0] <- 1  # Ensure thresholds_local are at least 1
  
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

########################## Complex simulation --------------------------------





Ts = 4
N_redes = 5
topolog_redes = 4
homoph_val = 6
alpha_val = 25
n_PLci_nodes = 5

Ts * N_redes * topolog_redes * homoph_val * alpha_val * n_PLci_nodes / topolog_redes






# Ejemplo de uso
N <- 324
density <- 0.029
num_edges <- round(density * N * (N - 1) / 2)
m <- 5
p <- 0.2

#graphs <- generate_small_world_networks(N, num_edges, p) #OK
#graphs <- generate_scale_free_networks(N, m)             #OK
graphs <- generate_erdos_renyi_networks(N, num_edges)

graph_file <- 'erdos-renyi'

# Parameter space
homoph_values <- seq(0.0, 0.3, length.out = 6)#seq(0.0, 0.5, by = 0.02)
alpha_values <- seq(0.0, 1.0, by = 0.02)
num_adopters_min <- 0.1

# Set threshold type
T_dist <- c("homo", "hetero")[1]
T_type <- c("abs", "frac")[2]

threshold_list <- c(0.10, 0.15, 0.20, 0.25)

for (threshold in threshold_list){
  
  T_str <- as.character(threshold)
  
  for (i in seq_along(graphs)) {
    graph <- graphs[[i]]
    thresholds <- replicate(N, threshold) #thresholds <- replicate(N, runif(1, 0.1,0.5))
    if(T_type == "frac"){thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(graph, x, mode = "total")))))}
    num_seeds_to_add <- thresholds - 1
    num_seeds_to_add[num_seeds_to_add < 0] <- 0
    thresholds[thresholds <= 0] <- 1
    
    gmat <- as.matrix(as_adjacency_matrix(graph))
    complex_centrality <- vector(mode = "list", length = N)
    complex_centrality <- lapply(1:N, function(x) get_complex(x, N, graph, gmat, thresholds, num_seeds_to_add, model_output_list))
    complex_centrality <- sapply(complex_centrality, function(x) tail(x, n = 1))
    top_PLci <- order(complex_centrality, decreasing = TRUE)[1:20]
    
    plot <- FALSE
    
    j <- 0  # Inicializar j como 0
    failed_nodes <- c()
    
    execution_times <- list()
    start_time <- Sys.time()
    for (seed_node in top_PLci) {
      if (j >= 5) {  # Verificar si ya se han encontrado 5 simulaciones exitosas
        break
      }
      
      temp_results <- sweep_homoph_parameter(seed_node, homoph_values, N, graph, thresholds,
                                             num_seeds_to_add, plot, file_dir, alpha_values, graph_file,
                                             'PLci', T_str, i)
      if (any(temp_results$num_adopters / N > num_adopters_min)) {
        results <- temp_results 
        cat("Successful infection.\n")
        j <- j + 1
      } else {
        failed_nodes <- c(failed_nodes, seed_node)
        cat("Unsuccessful infection.\n")
      }
    }
    
    end_time <- Sys.time()
    simulation_time <- end_time - start_time
    cat("Simulation time for graph", i, ":", simulation_time, "\n")
    
    # Save failed nodes
    failed_nodes[[i]] <- failed_nodes
  }
}

#saveRDS(failed_nodes, file = sprintf("Simulation_Stats_TFrac_2/%s/failed_nodes.rds", graph_file))

########################## Plot  -----------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

threshold <- 0.25
T_str <- as.character(threshold)

# Directorio de datos
data_dir <- sprintf("Simulation_Stats_TFrac_2/%s/1inf_N%s_hom_0.00_0.30_by_0.06_alpha_0.00_1.00_by_0.02_T%s", graph_file, N, T_str)
file_list <- list.files(path = data_dir, pattern = "*.csv", full.names = TRUE)

# Leer y combinar los datos
all_data <- lapply(file_list, read_csv)
combined_data <- bind_rows(all_data)

# Calcular la media de num_adopters/N para cada combinación de alpha_1 y homophily
summary_data <- combined_data %>%
  group_by(alpha_1, homophily) %>%
  summarize(mean_adopters = mean(num_adopters / N))

fig <- ggplot() +
  geom_line(data = combined_data, aes(x = alpha_1, y = num_adopters / N, group = interaction(seed, homophily), color = as.factor(homophily)), size = 0.5, alpha = 0.4) +
  #geom_line(data = summary_data, aes(x = alpha_1, y = mean_adopters, color = as.factor(homophily)), linewidth = 1.5) +
  geom_point(data = combined_data, aes(x = alpha_1, y = num_adopters / N, color = as.factor(homophily)), size = 1.5, alpha = 0.4) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.9, direction = -1, 
                        labels = c('0.00', '0.06', '0.12', '0.18', '0.24', '0.30')) +
  labs(x = expression(Gamma),  # Usar LaTeX para la etiqueta del eje-x
       y = "Adopters",
       color = "Homophily") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.33),  # Posición de la leyenda dentro del gráfico
        legend.direction = "vertical",  # Orientación vertical de la leyenda
        legend.box = "vertical",
        plot.title = element_blank()) +  # Eliminar el título
  guides(shape = "none")  # Ocultar la leyenda de las formas

print(fig)

# Guardar la figura como PDF
output_file_pdf <- sprintf("Simulation_Stats_TFrac_2/%s/plot_T%s.pdf", graph_file, T_str)
ggsave(output_file_pdf, plot = fig, device = "pdf", width = 4, height = 3, dpi = 175)
