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

curves <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_ExplCoRa_ATFv1_p01_Fig1_mY_mY.txt", col_names = TRUE)

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
# Generar gradiente continuo con 256 colores
color_palette <- colorRampPalette(c("#215EA1", "#22BAB9", "#FF5A5A",  "#FFF970"))(256)
color_palette <- colorRampPalette(c("#215EA1", "#22BAB9", "#94E2EA", "#FFC7D8", "#F15E8C", "#FF5A5A"))(256)

p <- plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(
               size = 10, 
               color = ~robustness,
               colorscale = setNames(data.frame(
                 seq(0, 1, length.out = 256),
                 color_palette
               ), NULL),
               colorbar = list(title = "Robustness"),
               showscale = TRUE,
               cmin = min(curves$robustness),
               cmax = max(curves$robustness)
             ),
             type = "scatter3d", 
             mode = "markers")  %>%
  layout(scene = list(
    xaxis = list(title = "mU", type = "log"),
    yaxis = list(title = "mW", type = "log"),
    zaxis = list(title = "eP", type = "log")
  ))

p


curves <- curves %>% filter(robustness>0.9)

curves[1:3] <- log10(curves[1:3])

hist(curves$mU)
param_combinations <- combn(names(curves)[1:3], 2, simplify = FALSE)



for (i in seq_along(param_combinations)) {
  # Crear el gráfico
  p <- ggplot(curves) +
    geom_point(aes(
      x = .data[[param_combinations[[i]][1]]],
      y = .data[[param_combinations[[i]][2]]],
      color = log10(steady_state)
    )) +
    scale_color_viridis_c(name = expression(log[10](steady_state))) +  # Escala con leyenda visible y título
    xlim(-3, 3) +
    ylim(-3, 3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(param_combinations[[i]][1], "vs", param_combinations[[i]][2]),
      x = param_combinations[[i]][1],
      y = param_combinations[[i]][2]
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"  # Muestra la barra de color a la derecha
    )
  
  # Crear nombre del archivo
  file_name <- paste0(
    "/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_",
    param_combinations[[i]][1], "_vs_",
    param_combinations[[i]][2], ".png"
  )
  
  # Guardar el gráfico
  ggsave(
    filename = file_name,
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  # Mensaje de confirmación
  cat("Guardado:", file_name, "\n")
}





# ---- Optimize ----

optimization_history <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p01_Fig2_mY_mY.txt", col_names = TRUE)

hist(optimization_history$robustness)
hist(optimization_history$`min(CoRA)`)
hist(optimization_history$mU, xlim=c(-3,3))

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

optimization_history <- optimization_history %>% distinct(mU, mW, .keep_all = TRUE)


optimization_history[2:4] <- log10(optimization_history[2:4])

param_combinations <- combn(names(optimization_history)[2:4], 2, simplify = FALSE)


for(i in seq_along(param_combinations)) {
  # Crear el gráfico
  p <- ggplot(optimization_history) +
    geom_point(aes(x = .data[[param_combinations[[i]][1]]], 
                   y = .data[[param_combinations[[i]][2]]], 
                   color = robustness)) +
    scale_color_gradientn(colors = c("#215EA1", "#22BAB9", "#94E2EA", "#FFC7D8", "#F15E8C", "#FF5A5A")) +
    guides(color = "none") +
    xlim(-3, 3) +  
    ylim(-3, 3) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "#E5ECF6"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_blank()  # Elimina las etiquetas de los ejes
    ) 
  # Crear nombre del archivo
  file_name <- paste0("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_", param_combinations[[i]][1], "_vs_", 
                      param_combinations[[i]][2], ".png")
  
  # Guardar el gráfico
  ggsave(filename = file_name, 
         plot = p, 
         width = 8, 
         height = 5, 
         dpi = 300)
  
  # Mensaje de confirmación
  cat("Guardado:", file_name, "\n")
}



##### Seady state
c("#FFF970", "#B2D51B", "#8AD51B")
c("#0D7F6D", "#21B15A", "#8AD51B", "#FFF970")
c("#FF5A5A", "#F15E8C", "#FFC7D8", "#DBFCFF", "#94E2EA", "#22BAB9", "#215EA1")
for(i in seq_along(param_combinations)) {
  # Crear el gráfico
  p <- ggplot(optimization_history) +
    geom_point(aes(x = .data[[param_combinations[[i]][2]]], 
                   y = .data[[param_combinations[[i]][1]]], 
                   color = log10(SS))) +
    scale_color_gradientn(colors = c("#215EA1", "#22BAB9", "#FF5A5A",  "#FFF970")
) +
    guides(color = "none") +
    xlim(-3, 3) +  
    ylim(-3, 3) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "#E5ECF6"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      #axis.title = element_blank()  # Elimina las etiquetas de los ejes
    ) 
  # Crear nombre del archivo
  file_name <- paste0("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_", param_combinations[[i]][2], "_vs_", 
                      param_combinations[[i]][1], ".png")
  
  # Guardar el gráfico
  ggsave(filename = file_name, 
         plot = p, 
         width = 8, 
         height = 5, 
         dpi = 300)
  
  # Mensaje de confirmación
  cat("Guardado:", file_name, "\n")
}









for(i in seq_along(param_combinations)) {
  # Crear el gráfico
  p <- ggplot(optimization_history) +
    geom_point(aes(x = .data[[param_combinations[[i]][2]]], 
                   y = .data[[param_combinations[[i]][1]]], 
                   color = log10(SS))) +
    scale_color_gradientn(
      colors = c("#215EA1", "#22BAB9", "#FF5A5A",  "#FFF970"),
      name = "log(SS)"  # Etiqueta del gradiente
    ) +
    xlim(-3, 3) +  
    ylim(-3, 3) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "#E5ECF6"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_blank(),  # Elimina las etiquetas de los ejes
      legend.position = "right",     # Posición de la leyenda
      legend.title = element_text(size = 10),  # Tamaño del título de la leyenda
      legend.text = element_text(size = 8)     # Tamaño del texto de la leyenda
    ) 
  
  # Crear nombre del archivo
  file_name <- paste0("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_", param_combinations[[i]][2], "_vs_", 
                      param_combinations[[i]][1], ".png")
  
  # Guardar el gráfico
  ggsave(filename = file_name, 
         plot = p, 
         width = 8, 
         height = 5, 
         dpi = 300)
  
  # Mensaje de confirmación
  cat("Guardado:", file_name, "\n")
}






######  Density

# Crear el gráfico de densidad
p_density <- ggplot(optimization_history, aes(x = mU)) +
  geom_density(fill = "#215EA1", alpha = 0.6, color = NA) +
  xlim(-3, 3) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank()
  )

# Crear nombre del archivo para la densidad
file_name_density <- "/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_robustness_density.png"

# Guardar el gráfico
ggsave(filename = file_name_density, 
       plot = p_density, 
       width = 8, 
       height = 5, 
       dpi = 300)

# Mensaje de confirmación
cat("Guardado:", file_name_density, "\n")


#### Iteraciones 
optimization_history <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p03_Fig2_mY_mY.txt", col_names = TRUE)
optimization_history[2:4] <- log10(optimization_history[2:4])
optim <- optimization_history[1:4]



library(ggplot2)
library(tidyr)
library(patchwork)

# Agregar valor de la función objetivo
historial_df <- optimization_history[1:501, ]

# Gráfico dual: parámetros y valor de función
p1 <- ggplot(historial_df %>% pivot_longer(c(mU, mW, eP)), 
             aes(x = Iteration, y = value, color = name)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "mU" = "#215EA1", 
      "mW" = "#22BAB9",   
      "eP" = "#94E2EA"    
    ),
    name = "Parameter"
  ) +
  labs(y = "Parameter") +
  theme_minimal()
p1
p2 <- ggplot(historial_df, aes(x = Iteration)) +
  geom_line(aes(y = robustness, color = "Robustness"), linewidth = 1) +
  geom_line(aes(y = historial_df$`min(CoRA)`, color = "min(CoRa)"), linewidth = 1) +
  labs (y = "f(x)") +
  scale_color_manual(values = c("Robustness" = "#F25E8C", "min(CoRa)" = "#FFC7D8"), name = "Metrics") +
  theme_minimal()




# Definir el vector de iteraciones una sola vez
iteraciones_seleccionadas <- c(1, 150, 300, 450, 500)

sub <- historial_df %>% filter(Iteration %in% iteraciones_seleccionadas)

library(ggplot2)
library(dplyr)
library(tidyr)

plot_curves2D <- function(data, titulo, iteraciones) {
  # Procesamiento de datos
  data_processed <- data %>%
    mutate(curve_id = row_number()) %>%
    pivot_longer(
      cols = 10:38,
      names_to = "curve",
      values_to = "CoRa"
    ) %>%
    mutate(curve = as.numeric(curve)) %>%
    group_by(Iteration) %>%
    mutate(unique_par = first(robustness)) %>%
    ungroup() %>%
    # Asegurar que las iteraciones estén en el orden correcto
    mutate(Iteration_factor = factor(Iteration, levels = iteraciones))
  
  # Generar colores automáticamente - gradiente de gris a azul
  n_iter <- length(iteraciones)
  colores_gris_azul <- colorRampPalette(c("#215EA1", "#22BAB9", "#94E2EA", "#FFC7D8", "#F15E8C"))(n_iter)
  names(colores_gris_azul) <- as.character(iteraciones)
  
  # Crear el gráfico con ggplot2
  ggplot(data_processed, aes(x = log10(curve), y = CoRa, group = Iteration_factor, 
                             color = Iteration_factor)) +
    geom_line(alpha = 0.8, linewidth = 1) +
    scale_color_manual(
      values = colores_gris_azul,
      name = "Iteration"
    ) +
    labs(
      x = "mY",
      y = "CoRa",
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_line(color = "#E5ECF6"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = unit(0.5, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    ylim(0, 1)
}

fig1 <- plot_curves2D(sub, "F1", iteraciones_seleccionadas)


p1/p2/fig1






# Leer los tres archivos
opt_p01 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p01_Fig2_mY_mY.txt", col_names = TRUE)
opt_p02 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p02_Fig2_mY_mY.txt", col_names = TRUE)
opt_p03 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p03_Fig2_mY_mY.txt", col_names = TRUE)
opt_p04 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p04_Fig1_mY_mY.txt", col_names = TRUE)
opt_p05 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p05_Fig1_mY_mY.txt", col_names = TRUE)
opt_p06 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p06_Fig1_mY_mY.txt", col_names = TRUE)
opt_p07 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p07_Fig1_mY_mY.txt", col_names = TRUE)

# Procesar cada dataset
opt_p01 <- opt_p01 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p02 <- opt_p02 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p03 <- opt_p03 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p04 <- opt_p04 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p05 <- opt_p05 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p06 <- opt_p06 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p07 <- opt_p07 %>% distinct(mU, mW, .keep_all = TRUE)

# Aplicar log10 a las columnas relevantes
opt_p01[2:4] <- log10(opt_p01[2:4])
opt_p02[2:4] <- log10(opt_p02[2:4])
opt_p03[2:4] <- log10(opt_p03[2:4])
opt_p04[2:4] <- log10(opt_p04[2:4])
opt_p05[2:4] <- log10(opt_p05[2:4])
opt_p06[2:4] <- log10(opt_p06[2:4])
opt_p07[2:4] <- log10(opt_p07[2:4])

# Agregar identificador de grupo a cada dataset
opt_p01$group <- "p01"
opt_p02$group <- "p02"
opt_p03$group <- "p03"
opt_p04$group <- "p04"
opt_p05$group <- "p05"
opt_p06$group <- "p06"
opt_p07$group <- "p07"

# Combinar todos los datos
combined_data <- bind_rows(opt_p01, opt_p02, opt_p03, opt_p04, opt_p05, opt_p06, opt_p07)

# Crear el gráfico de densidad superpuesto
p_density <- ggplot(combined_data, aes(x = mU, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 0.8) +
  xlim(-3, 3) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("#215EA1", "#E31A1C", "#33A02C", "black", "yellow", "purple", "orange")) +
  scale_color_manual(values = c("#215EA1", "#E31A1C", "#33A02C", "black", "yellow", "purple", "orange")) +
  labs(title = "Distribuciones de mU superpuestas",
       x = "mU (log10)",
       y = "Densidad")

# Crear nombre del archivo para la densidad
file_name_density <- "/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_robustness_density_combined.png"

# Guardar el gráfico
ggsave(filename = file_name_density, 
       plot = p_density, 
       width = 8, 
       height = 5, 
       dpi = 300)

# Mensaje de confirmación
cat("Guardado:", file_name_density, "\n")

# Mostrar el gráfico
print(p_density)













# Leer los tres archivos
opt_p01 <- read_tsv("./Output/OUT_OptCoRa_FADv1_p01_Fig1_mY_mY.txt", col_names = TRUE)
opt_p02 <- read_tsv("./Output/OUT_OptCoRa_FADv1_p02_Fig1_mY_mY.txt", col_names = TRUE)
opt_p03 <- read_tsv("./Output/OUT_OptCoRa_FADv1_p03_Fig1_mY_mY.txt", col_names = TRUE)

# Procesar cada dataset
opt_p01 <- opt_p01 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p02 <- opt_p02 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p03 <- opt_p03 %>% distinct(mU, mW, .keep_all = TRUE)

# Aplicar log10 a las columnas relevantes
opt_p01[2:4] <- log10(opt_p01[2:4])
opt_p02[2:4] <- log10(opt_p02[2:4])
opt_p03[2:4] <- log10(opt_p03[2:4])




# Agregar identificador de grupo a cada dataset
opt_p01$group <- "p01"
opt_p02$group <- "p02"
opt_p03$group <- "p03"

# Combinar todos los datos
combined_data <- bind_rows(opt_p01, opt_p02, opt_p03)

# Crear el gráfico de densidad superpuesto
p_density <- ggplot(combined_data, aes(x = mU, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 0.8) +
  xlim(-3, 3) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("#215EA1", "#E31A1C", "#33A02C")) +
  scale_color_manual(values = c("#215EA1", "#E31A1C", "#33A02C")) +
  labs(title = "Distribuciones de mU superpuestas",
       x = "mU (log10)",
       y = "Densidad")

# Crear nombre del archivo para la densidad
file_name_density <- "/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_robustness_density_combined.png"

# Guardar el gráfico
ggsave(filename = file_name_density, 
       plot = p_density, 
       width = 8, 
       height = 5, 
       dpi = 300)

# Mensaje de confirmación
cat("Guardado:", file_name_density, "\n")

# Mostrar el gráfico
print(p_density)








#### 
opt_p01 <- read_tsv("./Output/OUT_OptCoRa_ATFv1_p01_Fig3_mY_mY.txt", col_names = TRUE)
opt_p02 <- read_tsv("./Output/OUT_OptCoRa_FADv1_p01_Fig2_mY_mY.txt", col_names = TRUE)

opt_p01 <- opt_p01 %>% distinct(mU, mW, .keep_all = TRUE)
opt_p02 <- opt_p02 %>% distinct(mU, mW, .keep_all = TRUE)

opt_p01[2:8] <- log10(opt_p01[2:8])
opt_p02[2:8] <- log10(opt_p02[2:8])

opt_p01$group <- "p01"
opt_p02$group <- "p02"

# Combinar todos los datos
combined_data <- bind_rows(opt_p01, opt_p02)

# Crear el gráfico de densidad superpuesto
p_density <- ggplot(combined_data, aes(x = mU, fill = group, color = group)) +
  geom_density(alpha = 0.4, size = 0.8) +
  xlim(-3, 3) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "#E5ECF6"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("#215EA1", "#E31A1C", "#33A02C")) +
  scale_color_manual(values = c("#215EA1", "#E31A1C", "#33A02C")) +
  labs(title = "Distribuciones de mU superpuestas",
       x = "mU (log10)",
       y = "Densidad")

# Crear nombre del archivo para la densidad
file_name_density <- "/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/plot_robustness_density_combined.png"

# Guardar el gráfico
ggsave(filename = file_name_density, 
       plot = p_density, 
       width = 8, 
       height = 5, 
       dpi = 300)

# Mensaje de confirmación
cat("Guardado:", file_name_density, "\n")

# Mostrar el gráfico
print(p_density)



