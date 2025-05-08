library(dbscan)
library(plotly)
library(dplyr)
library(readr)
library(scales)
library(corrplot)
library(geometry)

# Crear la paleta de colores viridis para eps
num_CoRas <- 50  # Número de colores en la paleta
viridis_colors <- viridis_pal()(num_CoRas)

### ATF v1
curves <- read_csv("/media/rodrigo/Nuevo vol/OptCoRa/2025_03_30/p1/bifurcation_results_ATFv1_f.csv", col_names = FALSE)
colnames(curves) <- c("g", "mY", "gY", "mU", "gU", "mW", "gW", "e0", "eP", "eM", "mUs", "eps","os", "ss")

# Asignar colores: rojo si os == 1, de la paleta viridis si no
curves$color <- ifelse(curves$os > 0, "purple", viridis_colors[cut(curves$eps, breaks = num_CoRas, labels = FALSE)])

plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
        marker = list(size = 8),
        type = "scatter3d", mode = "markers",
        color = I(curves$color)) %>%  # I() para indicar que ya es un color definido
  layout(scene = list(
    xaxis = list(title = "mU",  type = "log"),
    yaxis = list(title = "mW",  type = "log"),
    zaxis = list(title = "eP",  type = "log")
  ))

curves <- curves %>% filter(eps == num_CoRas)

p <- plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(size = 5),
             type = "scatter3d", mode = "markers",
             color = I(curves$color)) %>%  # I() para indicar que ya es un color definido
  layout(scene = list(
    xaxis = list(title = "mU", type = "log", range = c(-3, 3)),
    yaxis = list(title = "mW", type = "log", range = c(-3, 3)),
    zaxis = list(title = "eP", type = "log", range = c(-3, 3)),
    aspectmode = "manual",   # Mantener proporciones
    aspectratio = list(x = 1, y = 1, z = 1)  # Proporciones iguales
  ))
