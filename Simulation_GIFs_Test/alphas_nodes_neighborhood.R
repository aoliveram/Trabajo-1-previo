############################# Network construction ------------------------

# Definir el número de nodos y el grado de conectividad
N <- 250  # Número de nodos 
z <- 4   # Grado de conectividad (número de vecinos por nodo)

# Crear una lista de aristas para una red regular con z = 4
edges <- c()

for (i in 1:N) {
  edges <- c(edges, i, (i - 2 - 1) %% N + 1)
  #print(edges)
  edges <- c(edges, i, (i - 1 - 1) %% N + 1)
  #print(edges)
  edges <- c(edges, i, (i + 1 - 1) %% N + 1)
  #print(edges)
  edges <- c(edges, i, (i + 2 - 1) %% N + 1)
  #print(edges)
}

# Crear el grafo a partir de la lista de aristas
graph <- graph(edges = edges, directed = FALSE)


######### Small World graph

# Ejemplo de uso con un grafo de mundo pequeño (nei == numero de vecinos originales, p=0.1)
graph <- sample_smallworld(dim=1, N, nei = 4, p = 0.15, loops = FALSE, multiple = FALSE)

?sample_smallworld
######### Scale Free graph

library(fastnet)
graph <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)

plot(graph, 
     layout = layout_with_fr(graph), 
     #layout = layout_in_circle(graph),
     vertex.size = 3, 
     edge.width = 0.5, 
     vertex.label.cex = 0.4, 
     edge.color = "gray", 
     edge.arrow.size = 0)

############################# Betas definition -------------------------------


library(fastnet)
N <- 250
graph <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)

# Generar valores de beta_o siguiendo una distribución normal truncada N(0.5, 0.5/3)
set.seed(123)
#beta_o <- pmin(pmax(rnorm(N, mean = 0.5, sd = 0.5/3), 0), 1) - 0.5 # Normal distr.
beta_o <- replicate(N, runif(1, 0, 1))
hist(beta_o, breaks = N, main = "Histograma de beta_o", xlab = "beta_o", ylab = "Frecuencia", col = "lightblue")
plot(beta_o)

# Asignar colores basados en los valores de beta_o
color_palette <- colorRampPalette(c("white", "blue"))
node_colors_o <- color_palette(100)[as.numeric(cut(beta_o, breaks = 100))]


layout <- layout_with_fr(graph)
# Graficar la red en forma circular
plot(graph, layout = layout, main = "Red beta_o",
     vertex.size = 4.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.color = "lightgray",
     edge.arrow.size = 0, vertex.color = node_colors_o)

# Calcular beta_neighborhood para cada nodo
beta_neighborhood <- sapply(1:N, function(i) {
  neighbors <- neighbors(graph, i)
  neighbor_indices <- as.integer(neighbors)  # Asegurarse de que los índices sean enteros
  (beta_o[i] + sum(beta_o[neighbor_indices])) / (length(neighbor_indices) + 1)
})

hist(beta_neighborhood, breaks = N, main = "Histograma de beta_o", xlab = "beta_neighborhood", ylab = "Frecuencia", col = "lightblue")

plot(beta_o - beta_neighborhood)

node_colors_neighborhood <- color_palette(100)[as.numeric(cut(beta_neighborhood, breaks = 100))]

# Graficar la red 
plot(graph, layout = layout, main = "Red beta_neighborhood",
     vertex.size = 4.5,  edge.width = 0.3, vertex.label.cex = 0.5, edge.color = "lightgray",
     edge.arrow.size = 0, vertex.color = node_colors_neighborhood)

############################# Network atributte ------------------------------

# Asignar beta_o como atributo de los nodos
V(graph)$beta_o <- beta_o
V(graph)$beta_neighborhood <- beta_neighborhood

min(beta_o); max(beta_o)
min(beta_neighborhood); max(beta_neighborhood)

# Hay vecindarios con beta_neighborhood < 0.2 ?

sum(beta_neighborhood < 0.25)

# Nodos con beta_neighborhood < 0.25
allow_infection_1 <- which(beta_neighborhood < 0.25)

# Nodos con beta_neighborhood entre 0.25 y 0.5
allow_infection_2 <- which(beta_neighborhood >= 0.25 & beta_neighborhood < 0.5)

########### vecindarios disponibles

# Encontrar los vecinos de los nodos rojos
neighbors <- unique(unlist(neighborhood(graph, order = 1, nodes = allow_infection_1)))

# Combinar los nodos rojos y sus vecinos
neighborhoods_allow_infection_1 <- unique(c(allow_infection_1, neighbors))

# Inicializar todos los nodos como grises
color_nodes <- rep("grey", vcount(g))

# Pintar de rojo los nodos definidos y sus vecinos
color_nodes[neighborhoods_allow_infection_1] <- "red"

# Graficar el grafo con los nodos pintados de rojo
plot(graph, layout = layout, vertex.color = color_nodes,
     vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.color = "gray",
     edge.arrow.size = 0)

############################### Simulation funtions ---------------------------

library(doParallel)
library(animation)
library(magick)


# Función para simular la difusión compleja y guardar los plots
get_complex_plot <- function(seed, N, g, gmat, thresholds, num_seeds_to_add, plot_dir, layout) {
  gmat_simulation <- matrix(0, nrow = N, ncol = N)
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
  
  activated <- logical(N)
  activated[seeds] <- TRUE
  gmat_simulation[seeds, ] <- 1
  gmat_simulation[, seeds] <- 1
  gmat_simulation_run <- gmat * gmat_simulation
  
  spread <- TRUE
  step <- 1  # Variable para contar los pasos
  
  # Crear directorio para almacenar los plots si no existe
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  # Fijar el layout de la red
  #layout <- layout_with_fr(g)
  #layout <- layout_in_circle(graph)
  
  while (spread) {
    # Graficar la red y guardar como archivo PNG
    png(paste0(plot_dir, "/step_", sprintf("%03d", step), ".png"), width = 1200, height = 900)
    plot(g, layout = layout, vertex.color = ifelse(activated, "red", "grey"),
         vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.color = "gray",
         edge.arrow.size = 0, main = paste("Step", step))
    dev.off()
    
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    t <- which(activated)
    t_1 <- which(influence_activated)
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t)
    activated[t_full] <- TRUE
    adopters <- which(activated)
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
    
    step <- step + 1  # Incrementar el número de pasos
  }
  
  num_adopters <- sum(activated)
  complex_g <- graph_from_adjacency_matrix(gmat_simulation_run)
  all_distances <- distances(complex_g, seed, V(complex_g), mode = "out")
  all_distances[all_distances == "Inf"] <- 0
  PLci <- mean(all_distances)
  
  return(list(seed = seed, num_adopters = num_adopters, PLci = PLci, num_steps = step - 1))
}

# Función auxiliar para seleccionar semillas adicionales (asumiendo que existe)
clustered_seeding <- function(seeds, g, num_seeds_needed) {
  new_seeds <- sample(setdiff(V(g), seeds), num_seeds_needed)
  return(c(seeds, new_seeds))
}

############################### Simulation -----------------------------------

# Configurar clúster para computación en paralelo
cores <- detectCores() - 1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
clusterExport(cluster, list('neighbors', 'graph_from_adjacency_matrix', 'distances', 'V', 'clustered_seeding'))

# Directorio para almacenar los plots temporales
plot_dir <- "plots_betas_neighborhood"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

gmat <- as.matrix(as_adjacency_matrix(graph))
thresholds <- replicate(N, 2)
#thresholds <- replicate(N, runif(1, 0.1,0.5)) #for Fractional Thresholds
T_dist <- c("homo", "hetero")[1]  # Set the threshold distribution type to "homo" (homogeneous)
T_type <- c("abs", "frac")[1]  # Set the threshold type to "abs" (absolute)

# If the threshold type is fractional, adjust the thresholds based on the number of neighbors
if(T_type == "frac") {
  thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(graph, x, mode = "total")))))
} # mean 4.7; mean[0:200] 3.05

num_seeds_to_add <- thresholds - 1
num_seeds_to_add[num_seeds_to_add < 0] <- 0
thresholds[thresholds <= 0] <- 1

# Semilla para la difusión compleja (solo un nodo)
seed <- sample(allow_infection_1, 1)

# Run the model
results <- get_complex_plot(seed, N, graph, gmat, thresholds, num_seeds_to_add, plot_dir, layout)

# Detener el clúster después de la computación en paralelo
#stopCluster(cluster)

# Convertir los resultados a un data frame
results_df <- as.data.frame(results)
print(results_df)

# Crear el GIF animado
images <- list()
for (i in 1:results$num_steps) {
  filename <- sprintf("%s/step_%03d.png", plot_dir, i)
  img <- image_read(filename)
  images <- c(images, img)
}
animation <- image_animate(image_join(images), fps = 1)
image_write(animation, "simulation GIF/betas_neighborhood.gif")
