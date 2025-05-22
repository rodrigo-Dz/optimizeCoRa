library(dbscan)
library(plotly)
library(dplyr)
library(readr)
library(scales)
library(corrplot)
library(geometry)
library(viridisLite)  # para escala de colores

### EXPLORATION ###
curves <- read_tsv("./output/OUT_ExplCoRa_ATFv1_p02_Fig1_mY_mY.txt", col_names = TRUE)

cont_colors <- viridis(100)[as.numeric(cut(curves$`proportion<=eps`, breaks = 100))]
final_colors <- ifelse(curves$oscilations > 0, "purple", cont_colors)
curves$color <- final_colors

plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
        marker = list(size = 8, color = curves$color),
        type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = "mU", type = "log"),
    yaxis = list(title = "mW", type = "log"),
    zaxis = list(title = "eP", type = "log")
  ))


#### OPTIMIZATION ###
curves <- read_tsv("./output/OUT_OptCoRa_ATFv1_Fig1_mY_mY.txt", col_names = TRUE)

p <- plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(size = 5),
             type = "scatter3d", mode = "markers",
             color = curves$`proportion<=eps`) %>%  # I() para indicar que ya es un color definido
  layout(scene = list(
    xaxis = list(title = "mU", type = "log", range = c(-3, 3)),
    yaxis = list(title = "mW", type = "log", range = c(-3, 3)),
    zaxis = list(title = "eP", type = "log", range = c(-3, 3)),
    aspectmode = "manual",   # Mantener proporciones
    aspectratio = list(x = 1, y = 1, z = 1)  # Proporciones iguales
  ))
p
