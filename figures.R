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
curves <- read_tsv("./output/OUT_OptCoRa_ATFv1_Fig1_mY_mY.txt", col_names = TRUE)

# Asignar colores: rojo si os == 1, de la paleta viridis si no
curves$color <- ifelse(curves$oscilations > 0, "purple", viridis_colors[cut(curves$`|CoRa<=0.1|`, breaks = num_CoRas, labels = FALSE)])

plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
        marker = list(size = 8),
        type = "scatter3d", mode = "markers",
        color = I(curves$color)) %>%  # I() para indicar que ya es un color definido
  layout(scene = list(
    xaxis = list(title = "mU",  type = "log"),
    yaxis = list(title = "mW",  type = "log"),
    zaxis = list(title = "eP",  type = "log")
  ))

curves <- curves %>% filter(curves$`|CoRa<=0.1|` == num_CoRas)

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
p
