# Seeds selected by alpha's value, and requirement of alpha_1, alpha_2.
# With alphas(activated_neighboors).

library(igraph)
library(doParallel)
library(animation)
library(magick)

###################### Graph ---------------------------------------------
rows <- 16
cols <- 16
N <- rows * cols
set.seed(123)

#library(fastnet)
#graph <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)(density=0.03)
#graph <- erdos.renyi.game(N, p=0.05, directed=FALSE)

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
graph <- rewire_edges(graph, 0.1) #(density=0.0285)

edge_density(graph)

layout <- layout_on_grid(graph)
#layout = layout_with_fr(graph)

plot(graph, layout = layout, vertex.size = 3, edge.width = 0.5, vertex.label.cex = 0.5, 
     edge.color = "gray", edge.arrow.size = 0)

###################### Functions --------------------------------------------

# Definir la clase State
setClass("State",
         slots = list(value = "character"),
         prototype = list(value = "State1"))

# Función para crear instancias de la clase State
State <- function(value = "State1") {
  if (!value %in% c("State1", "State2", "State3")) {
    stop("Estado no válido. Use 'State1', 'State2' o 'State3'.")
  }
  new("State", value = value)
}

# Función para ploteo de la red
plot_network <- function(g, activated, N, layout, step, plot_dir, alpha_1, alpha_2, homoph) {
  vertex_color <- sapply(activated, function(x) {
    if (x@value == "State2") "red" else if (x@value == "State3") "blue" else "grey"
  })
  
  edge_color <- rep("gray", ecount(g))
  E(g)$color <- edge_color
  
  for (i in 1:N) {
    if (activated[[i]]@value == "State2") {
      neighbors_i <- neighbors(g, i, mode = "total")
      for (j in neighbors_i) {
        E(g, P = c(i, j))$color <- "red"
      }
    } else if (activated[[i]]@value == "State3") {
      neighbors_i <- neighbors(g, i, mode = "total")
      for (j in neighbors_i) {
        E(g, P = c(i, j))$color <- "blue"
      }
    }
  }
  
  png(paste0(plot_dir, "/step_", sprintf("%03d", step), ".png"), width = 1200, height = 900)
  plot(g, layout = layout, vertex.color = vertex_color,
       vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.arrow.size = 0, 
       main = sprintf("Alpha Seed 1: %.2f | Alpha Seed 2: %.2f | Homophily:  %.2f | Step: %d", alpha_1, alpha_2, homoph, step))
  dev.off()
}

# Función para simular la difusión compleja y guardar los plots
get_complex_plot <- function(seed1, seed2, N, g, thresholds1, thresholds2, num_seeds_to_add1, num_seeds_to_add2, plot, plot_dir, alpha_1, alpha_2, homoph) {
  
  set.seed(123)
  
  # Selección de semillas adicionales para la infección 1
  seeds1 <- as.numeric(neighbors(g, seed1, mode = "total"))
  num_seeds_i1 <- num_seeds_to_add1[seed1]
  num_seeds_needed1 <- num_seeds_i1 - length(seeds1)
  
  if (num_seeds_needed1 > 0) {
    seeds1 <- clustered_seeding(seeds1, g, num_seeds_needed1)
  } else if (length(seeds1) > 1) {
    seeds1 <- c(seed1, sample(seeds1, num_seeds_i1))
  } else {
    seeds1 <- c(seed1, seeds1)
  }
  
  # Selección de semillas adicionales para la infección 2
  seeds2 <- as.numeric(neighbors(g, seed2, mode = "total"))
  num_seeds_i2 <- num_seeds_to_add2[seed2]
  num_seeds_needed2 <- num_seeds_i2 - length(seeds2)
  
  if (num_seeds_needed2 > 0) {
    seeds2 <- clustered_seeding(seeds2, g, num_seeds_needed2)
  } else if (length(seeds2) > 1) {
    seeds2 <- c(seed2, sample(seeds2, num_seeds_i2))
  } else {
    seeds2 <- c(seed2, seeds2)
  }
  
  # Inicializar estados
  activated <- vector("list", N)
  for (i in 1:N) {
    activated[[i]] <- State("State1")
  }
  
  for (i in seeds1) {
    activated[[i]] <- State("State2")
  }
  
  for (i in seeds2) {
    activated[[i]] <- State("State3")
  }
  
  gmat_simulation_1 <- matrix(0, nrow = N, ncol = N)
  gmat_simulation_2 <- matrix(0, nrow = N, ncol = N)
  
  gmat_simulation_1[seeds1, ] <- 1
  gmat_simulation_1[, seeds1] <- 1
  gmat_simulation_2[seeds2, ] <- 1
  gmat_simulation_2[, seeds2] <- 1
  gmat <- as.matrix(as_adjacency_matrix(g))
  gmat_simulation_run_1 <- gmat * gmat_simulation_1
  gmat_simulation_run_2 <- gmat * gmat_simulation_2
  spread <- TRUE
  step <- 1 # Variable para contar los pasos
  
  # Fijar el layout de la red
  layout <- layout_on_grid(g)
  
  while (spread) {
    print(step)
    if (plot) {plot_network(g, activated, N, layout, step, plot_dir, alpha_1, alpha_2, homoph)}
    
    # Determinar nuevos infectados para la infección 1
    influence1 <- colSums(gmat_simulation_run_1)
    influence_activated_1 <- influence1 >= thresholds1
    t1 <- which(sapply(activated, function(x) x@value) == "State2")
    t1_1 <- which(influence_activated_1 & sapply(activated, function(x) x@value) == "State1" & V(graph)$alpha < alpha_1)
    
    # Nueva regla: permitir infección si alpha es similar al de los vecinos
    t1_2 <- which(influence_activated_1 & sapply(activated, function(x) x@value) == "State1" & V(graph)$alpha >= alpha_1)
    for (i in t1_2) {
      vecinos <- neighbors(g, i)
      similar_alpha_neighbors <- sum(abs(V(graph)$alpha[i] - V(graph)$alpha[vecinos]) < homoph)
      if (similar_alpha_neighbors >= thresholds1[i]) {
        t1_1 <- c(t1_1, i)
      }
    }
    
    t1_full <- union(t1, t1_1)
    spread1 <- length(t1_full) > length(t1)
    
    # Determinar nuevos infectados para la infección 2
    influence2 <- colSums(gmat_simulation_run_2)
    influence_activated_2 <- influence2 >= thresholds2
    t2 <- which(sapply(activated, function(x) x@value) == "State3")
    t2_1 <- which(influence_activated_2 & sapply(activated, function(x) x@value) == "State1" & V(graph)$alpha > alpha_2)
    
    # Nueva regla: permitir infección si alpha es similar al de los vecinos
    t2_2 <- which(influence_activated_2 & sapply(activated, function(x) x@value) == "State1" & V(graph)$alpha <= alpha_2)
    for (i in t2_2) {
      vecinos <- neighbors(g, i)
      similar_alpha_neighbors <- sum(abs(V(graph)$alpha[i] - V(graph)$alpha[vecinos]) < homoph)
      if (similar_alpha_neighbors >= thresholds2[i]) {
        t2_1 <- c(t2_1, i)
      }
    }
    
    t2_full <- union(t2, t2_1)
    spread2 <- length(t2_full) > length(t2)
    
    spread <- spread1 || spread2
    
    for (i in t1_full) {
      activated[[i]] <- State("State2")
    }
    
    for (i in t2_full) {
      activated[[i]] <- State("State3")
    }
    
    adopters_1 <- which(sapply(activated, function(x) x@value == "State2"))
    adopters_2 <- which(sapply(activated, function(x) x@value == "State3"))
    gmat_simulation_1[adopters_1, ] <- 1
    gmat_simulation_1[, adopters_1] <- 1
    gmat_simulation_2[adopters_2, ] <- 1
    gmat_simulation_2[, adopters_2] <- 1
    gmat_simulation_run_1 <- gmat * gmat_simulation_1
    gmat_simulation_run_2 <- gmat * gmat_simulation_2
    
    # Verificar si todos los nodos están infectados
    all_infected <- all(sapply(activated, function(x) x@value) != "State1")
    if (all_infected) {
      spread <- FALSE
    }
    
    step <- step + 1 # Incrementar el número de pasos
  }
  
  # Graficar el estado final de la red
  if (plot) {plot_network(g, activated, N, layout, step, plot_dir, alpha_1, alpha_2, homoph)}
  
  num_adopters_1 <- sum(sapply(activated, function(x) x@value == "State2"))
  num_adopters_2 <- sum(sapply(activated, function(x) x@value == "State3"))
  return(list(homophily = homoph, seed1 = seed1, num_adopters_1 = num_adopters_1, seed2 = seed2, num_adopters_2 = num_adopters_2, num_steps = step))
}

# Función auxiliar para seleccionar semillas adicionales (asumiendo que existe)
clustered_seeding <- function(seeds, g, num_seeds_needed) {
  new_seeds <- sample(setdiff(V(g), seeds), num_seeds_needed)
  return(c(seeds, new_seeds))
}

###################### Complex simulation --------------------------------

# Configurar clúster para computación en paralelo
cores <- detectCores() - 1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
clusterExport(cluster, list('plot_network', 'neighbors', 'graph_from_adjacency_matrix',
                            'distances', 'V', 'clustered_seeding'))

# Directorio para almacenar los plots temporales
plot_dir <- "plots_TEST"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Crear un valor de alpha arbitrario entre 0 y 1 para cada nodo
set.seed(123)  # Para reproducibilidad
alpha <- runif(N, 0, 1)
V(graph)$alpha <- alpha

# Definir umbrales para las infecciones
thresholds1 <- replicate(N, 2)
thresholds2 <- replicate(N, 2)

num_seeds_to_add1 <- thresholds1 - 1
num_seeds_to_add1[num_seeds_to_add1 < 0] <- 0
thresholds1[thresholds1 <= 0] <- 1

num_seeds_to_add2 <- thresholds2 - 1
num_seeds_to_add2[num_seeds_to_add2 < 0] <- 0
thresholds2[thresholds2 <= 0] <- 1

# Seleccionar semillas basadas en el valor de alpha
alpha_1 <- 0.5
alpha_2 <- 0.5
seed1_candidates <- V(graph)[alpha < alpha_1]
seed2_candidates <- V(graph)[alpha > alpha_2]

#set.seed(123)
seed1 <- as.integer(sample(seed1_candidates, 1))
seed2 <- as.integer(sample(seed2_candidates, 1))

# Run the model
plot <- TRUE

# Define the values of homoph to sweep
homoph_values <- seq(0.0, 0.5, by = 0.1)

# Define the path of the main folder
homoph_min <- sprintf("%02d", min(homoph_values) * 10)
homoph_max <- sprintf("%02d", max(homoph_values) * 10)
parameter_sweep_homoph_folder <- sprintf("parameter_sweep_homoph_%s_%s", homoph_min, homoph_max)

# Create the main folder if it does not exist
if (!dir.exists(parameter_sweep_homoph_folder)) dir.create(parameter_sweep_homoph_folder)

# Initialize an empty list to store data frames
list_of_dfs <- list()

# Iterate over each value of homoph
for (i in seq_along(homoph_values)) {
  homoph <- homoph_values[i]
  print(paste("Running simulation for homophily: ", homoph))
  results <- get_complex_plot(seed1, seed2, N, graph, thresholds1, thresholds2, 
                              num_seeds_to_add1, num_seeds_to_add2, plot, plot_dir, alpha_1, alpha_2, homoph)
  results_df <- as.data.frame(results)
  print(results_df)
  write.csv(results_df, sprintf("%s/results_%02d.csv", parameter_sweep_homoph_folder, i), row.names = FALSE)
  list_of_dfs <- append(list_of_dfs, list(results_df))
  images <- list()
  for (j in 1:results$num_steps) {
    img <- image_read(sprintf("%s/step_%03d.png", plot_dir, j))
    images <- c(images, img)
  }
  image_write(image_animate(image_join(images), fps = 1), sprintf("%s/complex_base_barrido_homoph_%02d.gif", parameter_sweep_homoph_folder, i))
}

# Combine all data frames into one and assign to a dynamically named variable
assign(sprintf("results_%s_%s", homoph_min, homoph_max), do.call(rbind, list_of_dfs))

# Save the combined data frame in the main folder
write.csv(get(sprintf("results_%s_%s", homoph_min, homoph_max)), sprintf("%s/results_%s_%s.csv", parameter_sweep_homoph_folder, homoph_min, homoph_max), row.names = FALSE)

# Print the combined data frame
print(get(sprintf("results_%s_%s", homoph_min, homoph_max)))

################################## TEST ---------------------------------

sweep_homoph_parameter <- function(seed1, seed2, N, graph, thresholds1, thresholds2, 
                                   num_seeds_to_add1, num_seeds_to_add2, plot, plot_dir, alpha_1, alpha_2) {
  # Define the values of homoph to sweep
  homoph_values <- seq(0.0, 0.5, by = 0.1)
  
  # Define the path of the main folder
  homoph_min <- sprintf("%02d", min(homoph_values) * 10)
  homoph_max <- sprintf("%02d", max(homoph_values) * 10)
  parameter_sweep_homoph_folder <- sprintf("parameter_sweep_homoph_%s_%s", homoph_min, homoph_max)
  
  # Create the main folder if it does not exist
  if (!dir.exists(parameter_sweep_homoph_folder)) dir.create(parameter_sweep_homoph_folder)
  
  # Initialize an empty list to store data frames
  list_of_dfs <- list()
  
  # Iterate over each value of homoph
  for (i in seq_along(homoph_values)) {
    homoph <- homoph_values[i]
    print(paste("Running simulation for homophily: ", homoph))
    results <- get_complex_plot(seed1, seed2, N, graph, thresholds1, thresholds2, 
                                num_seeds_to_add1, num_seeds_to_add2, plot, plot_dir, alpha_1, alpha_2, homoph)
    results_df <- as.data.frame(results)
    print(results_df)
    write.csv(results_df, sprintf("%s/results_%02d.csv", parameter_sweep_homoph_folder, i), row.names = FALSE)
    list_of_dfs <- append(list_of_dfs, list(results_df))
    if (plot) {
      images <- list()
      for (j in 1:results$num_steps) {
        img <- image_read(sprintf("%s/step_%03d.png", plot_dir, j))
        images <- c(images, img)
      }
      image_write(image_animate(image_join(images), fps = 1), sprintf("%s/complex_base_barrido_homoph_%02d.gif", parameter_sweep_homoph_folder, i))
    }
  }
  
  # Combine all data frames into one and assign to a dynamically named variable
  combined_results_df <- do.call(rbind, list_of_dfs)
  combined_results_name <- sprintf("results_%s_%s", homoph_min, homoph_max)
  assign(combined_results_name, combined_results_df, envir = .GlobalEnv)
  
  # Save the combined data frame in the main folder
  write.csv(get(combined_results_name), sprintf("%s/results_%s_%s.csv", parameter_sweep_homoph_folder, homoph_min, homoph_max), row.names = FALSE)
}

# Example usage
sweep_homoph_parameter(seed1, seed2, N, graph, thresholds1, thresholds2, num_seeds_to_add1, num_seeds_to_add2, plot, plot_dir, alpha_1, alpha_2)
