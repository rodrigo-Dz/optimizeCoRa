library(dbscan)
library(plotly)
library(dplyr)
library(readr)
library(scales)
library(geometry)
library(viridisLite)  
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(viridis)
library(patchwork)
library(ggplot2)


# ---- Explore ----

curves <- read_tsv("./Output/OUT_ExplCoRa_ATFv1_p02_Fig1_mY_mY.txt", col_names = TRUE)

# set scale [0, 1]
# curves$robustness[1] = 1

cont_colors <- viridis(100)[as.numeric(cut(curves$robustness, breaks = 100))]
final_colors <- ifelse(curves$oscilations != 0, "purple", cont_colors)
curves$color <- final_colors

curves <- curves %>% filter(negative_sol < 5 & other_errors < 5)

p <- plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(
               size = 10, 
               color = ~robustness,  # Usar robustness para la barra de color
               colorscale = "Viridis",
               colorbar = list(title = "Robustness"),
               showscale = TRUE,
               # Mapear los colores personalizados
               cmin = min(curves$robustness),
               cmax = max(curves$robustness)
             ),
             type = "scatter3d", 
             mode = "markers") %>%
  # Sobrescribir colores con tu esquema personalizado
  add_markers(color = ~I(final_colors)) %>%
  layout(scene = list(
    xaxis = list(title = "mU", type = "log"),
    yaxis = list(title = "mW", type = "log"),
    zaxis = list(title = "eP", type = "log")
  ))

p



# ---- Optimize ----

optimization_history <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_OptCoRa_ATFv1_p01_Fig1_mY_mY.txt", col_names = TRUE)

hist(optimization_history$robustness)
hist(optimization_history$`min(CoRA)`)

p <- plot_ly(optimization_history, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(
               size = 10, 
               color = ~robustness,  # Usar robustness para la barra de color
               colorscale = "Viridis",
               colorbar = list(title = "Robustness"),
               showscale = TRUE,
               # Mapear los colores personalizados
               cmin = min(optimization_history$robustness),
               cmax = max(optimization_history$robustness)
             ),
             type = "scatter3d", 
             mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = "mU", type = "log"),
    yaxis = list(title = "mW", type = "log"),
    zaxis = list(title = "eP", type = "log")
  ))

p


optimization_history[2:8] <- log10(optimization_history[2:8])

param_combinations <- combn(names(optimization_history)[2:8], 2, simplify = FALSE)

for(i in seq_along(param_combinations)) {
  # Crear el gráfico
  p <- ggplot(optimization_history) +
    geom_point(aes(x = .data[[param_combinations[[i]][1]]], 
                   y = .data[[param_combinations[[i]][2]]], 
                   color = robustness)) +
    scale_color_viridis_c() +
    guides(color = "none") +
    xlim(-3, 3) +  
    ylim(-3, 3) +
    labs(title = paste(param_combinations[[i]][1], "vs", param_combinations[[i]][2]))
  
  # Crear nombre del archivo
  file_name <- paste0("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_", param_combinations[[i]][1], "_vs_", 
                      param_combinations[[i]][2], ".png")
  
  # Guardar el gráfico
  ggsave(filename = file_name, 
         plot = p, 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  # Mensaje de confirmación
  cat("Guardado:", file_name, "\n")
}



