# Cargar librerías necesarias
library(igraph)
library(magick)

# Definir parámetros
num_nodes <- 250  # Número de nodos
N <- num_nodes
prob <- 0.1  # Probabilidad de que exista una arista entre dos nodos
num_steps <- 10  # Número de pasos de la simulación

library(fastnet)
network <- to.igraph(net.holme.kim(N, 4, 0.5)) # m=4, pt=0.5 (triad prob.)

# Inicializar estados de nodos (0: susceptible, 1: infectado por contagio 1, 2: infectado por contagio 2)
V(network)$state1 <- 0
V(network)$state2 <- 0

# Definir nodos iniciales infectados para ambos contagios
initial_infected1 <- sample(V(network), 5)
initial_infected2 <- sample(V(network), 5)

V(network)$state1[initial_infected1] <- 1
V(network)$state2[initial_infected2] <- 2

# Función para simular un paso de contagio complejo
simulate_step <- function(graph) {
  new_state1 <- V(graph)$state1
  new_state2 <- V(graph)$state2
  
  for (v in V(graph)) {
    if (V(graph)$state1[v] == 0 && V(graph)$state2[v] == 0) {
      neighbors <- neighbors(graph, v)
      infected_neighbors1 <- sum(V(graph)$state1[neighbors] == 1)
      infected_neighbors2 <- sum(V(graph)$state2[neighbors] == 2)
      
      if (infected_neighbors1 >= 2) {
        new_state1[v] <- 1
      } else if (infected_neighbors2 >= 2) {
        new_state2[v] <- 2
      }
    }
  }
  
  V(graph)$state1 <- new_state1
  V(graph)$state2 <- new_state2
  return(graph)
}

# Crear una lista para guardar las imágenes de cada paso
images <- list()

# Calcular el layout una vez
layout <- layout.fruchterman.reingold(network)

# Simular múltiples pasos
for (step in 1:num_steps) {
  network <- simulate_step(network)
  
  # Guardar la imagen del estado de la red en cada paso
  img <- image_graph(width = 600, height = 600, res = 96)
  plot(network, layout = layout, vertex.color = ifelse(V(network)$state1 == 1, "red",
                                                       ifelse(V(network)$state2 == 2, "blue", "white")),
       main = paste("Paso", step),
       vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.arrow.size = 0)
  dev.off()
  images[[step]] <- img
}

# Guardar la imagen del estado final de la red
img <- image_graph(width = 600, height = 600, res = 96)
plot(network, layout = layout, vertex.color = ifelse(V(network)$state1 == 1, "red",
                                                     ifelse(V(network)$state2 == 2, "blue", "white")),
     main = "Estado Final",
     vertex.size = 2.5, edge.width = 0.3, vertex.label.cex = 0.5, edge.arrow.size = 0)
dev.off()
images[[num_steps + 1]] <- img

# Crear el GIF
gif <- image_animate(image_join(images), fps = 1)
image_write(gif, "diffusion_process.gif")
