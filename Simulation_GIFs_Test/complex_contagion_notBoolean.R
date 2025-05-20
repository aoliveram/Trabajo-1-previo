library(igraph)
library(doParallel)
library(animation)
library(magick)

######### Scale Free graph

library(fastnet)
N <- 250
graph <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)

plot(graph, 
     layout = layout_in_circle(graph),
     vertex.size = 3, 
     edge.width = 0.5, 
     vertex.label.cex = 0.4, 
     edge.color = "gray", 
     edge.arrow.size = 0)

###################### Complex simulation --------------------------------

# Definir la clase State
setClass("State",
         slots = list(value = "character"),
         prototype = list(value = "State1"))

# Función para crear instancias de la clase State
State <- function(value = "State1") {
  if (!value %in% c("State1", "State2")) {
    stop("Estado no válido. Use 'State1' o 'State2'.")
  }
  new("State", value = value)
}

# Función para simular la difusión compleja y guardar los plots
get_complex_plot <- function(seed, N, g, gmat, thresholds, num_seeds_to_add, plot_dir) {
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
  
  # Crear una lista de instancias de la clase State
  activated <- vector("list", N)
  for (i in 1:N) {
    activated[[i]] <- State("State1")
  }
  for (i in seeds) {
    activated[[i]] <- State("State2")
  }
  
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
  layout <- layout_with_fr(g)
  #layout <- layout_in_circle(graph)
  
  while (spread) {
    print(step)
    # Graficar la red y guardar como archivo PNG
    png(paste0(plot_dir, "/step_", sprintf("%03d", step), ".png"), width = 1200, height = 900)
    plot(g, layout = layout, vertex.color = ifelse(sapply(activated, function(x) x@value) == "State2", "red", "grey"),
         vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.color = "gray",
         edge.arrow.size = 0, main = paste("Step", step))
    dev.off()
    
    influence <- colSums(gmat_simulation_run)
    influence_activated <- influence >= thresholds
    t <- which(sapply(activated, function(x) x@value) == "State2")
    t_1 <- which(influence_activated)
    t_full <- union(t, t_1)
    spread <- length(t_full) > length(t)
    for (i in t_full) {
      activated[[i]] <- State("State2")
    }
    adopters <- which(sapply(activated, function(x) x@value) == "State2")
    gmat_simulation[adopters, ] <- 1
    gmat_simulation[, adopters] <- 1
    gmat_simulation_run <- gmat * gmat_simulation
    
    step <- step + 1  # Incrementar el número de pasos
  }
  
  num_adopters <- sum(sapply(activated, function(x) x@value) == "State2")
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

# Configurar clúster para computación en paralelo
cores <- detectCores() - 1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
clusterExport(cluster, list('neighbors', 'graph_from_adjacency_matrix', 'distances', 'V', 'clustered_seeding'))

# Directorio para almacenar los plots temporales
plot_dir <- "plots"
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
seed <- 2

# Run the model
results <- get_complex_plot(seed, N, graph, gmat, thresholds, num_seeds_to_add, plot_dir)

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
image_write(animation, "complex_notBoolean_02.gif")
