# Instalar y cargar las librerías necesarias
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("patchwork")) install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Parámetros
N <- 324  # Ajusta este valor según tu caso
graph_files <- c('erdos-renyi', 'scale-free', 'small-world', 'small-world-SDA')
threshold_list <- c(0.10, 0.15, 0.20)
threshold_labels <- c('T=0.10', 'T=0.15', 'T=0.20')
graph_labels <- c('Erdös–Rényi', 'Scale-Free', 'Small-World', 'Small-World SDA')

# Crear una lista para almacenar los gráficos
plots <- list()

# Generar gráficos para cada combinación de graph_file y threshold
for (graph_file in graph_files) {
  for (threshold in threshold_list) {
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
    
    # Crear el gráfico
    fig <- ggplot() +
      geom_line(data = combined_data, aes(x = alpha_1, y = num_adopters / N, group = interaction(seed, homophily), color = as.factor(homophily)), size = 0.5, alpha = 0.4) +
      geom_point(data = combined_data, aes(x = alpha_1, y = num_adopters / N, color = as.factor(homophily)), size = 1.5, alpha = 0.4) +
      scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.9, direction = -1, 
                            labels = c('0.00', '0.06', '0.12', '0.18', '0.24', '0.30')) +
      labs(x = NULL, y = NULL, color = "Social Distance") +
      theme_minimal() +
      theme(legend.position = "none", plot.title = element_blank())
    
    # Añadir el gráfico a la lista
    plots[[paste(graph_file, T_str, sep = "_")]] <- fig
  }
}

# Asignar títulos y etiquetas según la posición
for (i in 1:3) {
  # Primera fila: añadir título con el valor de threshold
  plots[[i]] <- plots[[i]] + ggtitle(threshold_labels[i])
}

for (j in 1:4) {
  # Primera columna: añadir etiqueta del eje-y con el tipo de red
  plots[[1 + (j - 1) * 3]] <- plots[[1 + (j - 1) * 3]] + ylab(paste(graph_labels[j], "\nAdopters"))
}

for (i in 1:3) {
  plots[[i + 9]] <- plots[[i + 9]] + xlab(bquote(Gamma ~ "\n" ~ .(threshold_labels[i])))
}

# Combinar los gráficos en un arreglo 4x3
combined_plot <- (plots[[1]] | plots[[5]] | plots[[9]]) /
  (plots[[2]] | plots[[6]] | plots[[10]]) /
  (plots[[3]] | plots[[7]] | plots[[11]]) /
  (plots[[4]] | plots[[8]] | plots[[12]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Mostrar el gráfico combinado
print(combined_plot)

# Guardar el gráfico combinado como PDF
output_file_pdf <- "Simulation_Stats_TFrac_2/combined_plot_3.pdf"
ggsave(output_file_pdf, plot = combined_plot, device = "pdf", width = 7, height = 8, dpi = 175)
