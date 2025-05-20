rows <- 20; cols <- 20; N <- rows * cols

library(igraph)
library(ggplot2)

# Crear una función para calcular C(p)/C(0) y L(p)/L(0)
calculate_cl <- function(p) {
  graph <- sample_smallworld(dim = 1, size = N, nei = k/2, p = p)
  c <- transitivity(graph, type = "global")
  l <- mean_distance(graph, directed = FALSE)
  return(c(c, l))
}

# Valores de p
p_values <- 10^seq(-4, 0, length.out = 50)
cl_values <- sapply(p_values, calculate_cl)

# Normalizar C(p)/C(0) y L(p)/L(0)
c0 <- cl_values[1, 1]
l0 <- cl_values[2, 1]
cl_values[1, ] <- cl_values[1, ] / c0
cl_values[2, ] <- cl_values[2, ] / l0

# Crear el dataframe para ggplot
df <- data.frame(
  p = rep(p_values, 2),
  value = c(cl_values[1, ], cl_values[2, ]),
  type = rep(c("C(p)/C(0)", "L(p)/L(0)"), each = length(p_values))
)

# Graficar
ggplot(df, aes(x = p, y = value, shape = type, color = type)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Rewiring Probability (p)", y = "Normalized Value") +
  geom_vline(xintercept = c(0.05, 0.1, 0.2), linetype = "dashed", color = "black") +
  theme_minimal()
