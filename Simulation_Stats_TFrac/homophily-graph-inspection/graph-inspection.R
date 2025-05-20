library(igraph)

set.seed(123)
N <- 400
density <- 0.029
num_edges <- round(density * N * (N - 1) / 2)
k <- round(2 * num_edges / N) + 1 # Number of neighbors
p <- 0.3
synthetic_network <- sample_smallworld(dim = 1, size = N, nei = k/2, p = p)

edge_data <- read.csv("Talaga-homophily-network/talaga-homophily-network-edges.csv")
real_network <- graph_from_data_frame(edge_data, directed = FALSE, vertices = NULL)

clustering_synthetic <- transitivity(synthetic_network, type = "global")
clustering_real <- transitivity(real_network, type = "global")

assortativity_synthetic <- assortativity_degree(synthetic_network)
assortativity_real <- assortativity_degree(real_network)

avg_path_length_synthetic <- average.path.length(synthetic_network)
avg_path_length_real <- average.path.length(real_network)

# Imprimir los resultados
cat("Clustering:\n")
cat("Sintético:", clustering_synthetic, "\n")     #Sintético: 0.09828563 
cat("Real:", clustering_real, "\n\n")             #Real:      0.226257 

cat("Assortatividad de grado:\n")       
cat("Sintético:", assortativity_synthetic, "\n")  #Sintético: -0.02153456 
cat("Real:", assortativity_real, "\n\n")          #Real: 0.2331087 

cat("Longitud promedio de los caminos:\n")        
cat("Sintético:", avg_path_length_synthetic, "\n")#Sintético: 2.739148 
cat("Real:", avg_path_length_real, "\n")          #Real: 3.058208 

# Degree plot
synthetic_df <- data.frame(degree = degree(synthetic_network), network = 'Synthetic')
real_df <- data.frame(degree = degree(real_network), network = 'Real')

combined_df <- rbind(synthetic_df, real_df)

mean_synthetic <- mean(degree(synthetic_network))
std_synthetic <- sd(degree(synthetic_network))
mean_real <- mean(degree(real_network))
std_real <- sd(degree(real_network))

library(ggplot2)

fig <- ggplot(combined_df, aes(x = degree, fill = network)) +
  geom_histogram(binwidth = 1, color = "black", position = "identity", alpha = 0.5) +
  geom_vline(aes(xintercept = mean_synthetic), color = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = mean_real), color = "orange", linetype = "dashed") +
  labs(title = "Degree distribution",
       x = "Grado",
       y = "Frecuencia") +
  scale_fill_manual(values = c("blue", "orange")) +
  theme_minimal() +
  annotate("text", x = 25.4, y = 55, label = paste("Synthetic Mean:", round(mean_synthetic, 2), "±", round(std_synthetic, 2)), color = "orange") +
  annotate("text", x = 25.1, y = 45, label = paste("Real Mean:", round(mean_real, 2), "±", round(std_real, 2)), color = "blue")

output_file_pdf <- "Simulation_Stats_TFrac/homophily-graph-inspection/degree-distribution-comparison.pdf"
#save_image(fig, output_file_pdf, format = "pdf", width = 4*150, height = 3*150)
ggsave(output_file_pdf, plot = fig, device = "pdf", width = 16/2.9, height = 10/2.9)
