curves <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_OptCoRa_ATFv1_p02_Fig1_mY_mY.txt", col_names = TRUE)
# curves$`proportion<=eps`[1] = 1
cont_colors <- viridis(100)[as.numeric(cut(curves$robustness, breaks = 100))]
final_colors <- ifelse(curves$oscilations > 1, "purple", cont_colors)
curves$color <- final_colors

curves <- curves %>% filter(other_errors == 0)
curves <- curves[-1,]
p <- plot_ly(curves, x = ~mU, y = ~mW, z = ~eP, 
             marker = list(size = 5, color = curves$robustness),
             type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = "mU", type = "log"),
    yaxis = list(title = "mW", type = "log"),
    zaxis = list(title = "eP", type = "log")
  ))
p



library(GGally)
curves <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_OptCoRa_ATFv1_p01_Fig1_mY_mY.txt", col_names = TRUE)
df = curves[1:9]
# Excluir la columna iter (y opcionalmente escalar)
df[2:8] <- log10(df[2:8])

df_scaled <- as.data.frame(df[ , -1])

# Asegúrate de poner output al final o al principio
GGally::ggparcoord(df_scaled[2000:5000,],
                   columns = 1:(ncol(df_scaled)-1),
                   groupColumn = ncol(df_scaled)) +
  scale_color_viridis_c() +
  labs(title = "Parámetros vs output (parallel coordinates)") +
  theme_minimal()



library(plotly)
library(readr)

curves <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_OptCoRa_ATFv1_p01_Fig1_mY_mY.txt")

df <- curves[1:9]
df[2:8] <- log10(df[2:8])

# Parallel coordinates animado
plot_ly(
  type = 'parcoords',
  line = list(color = df$robustness, colorscale = 'Viridis'),
  dimensions = list(
    list(label = 'k1', values = df$mU),
    list(label = 'k2', values = df$mW),
    list(label = 'k3', values = df$eP),
    list(label = 'k4', values = df$gU),
    list(label = 'k5', values = df$gW),
    list(label = 'output', values = df$robustness)
  )
)





library(plotly)

df <- df[2000:2100, ]  # para no saturar

# Normalizar parámetros
params <- names(df)[2:8]
df[params] <- scale(df[params])
# --- Crear figura ---
fig <- plot_ly()

cols <- viridisLite::viridis(nrow(df))
cols <- cols[rank(df$robustness)]  # colores ordenados por output

fig <- plot_ly()
for (i in seq_len(nrow(df))) {
  fig <- fig |> add_trace(
    x = 1:length(params),
    y = as.numeric(df[i, params]),
    z = rep(df$iter[i], length(params)),
    type = "scatter3d",
    mode = "lines",
    line = list(color = cols[i], width = 3),
    showlegend = FALSE
  )
}


# --- Layout del gráfico ---
fig <- fig |>
  layout(
    scene = list(
      zaxis = list(title = "Parámetros", tickvals = 1:length(params), ticktext = params),
      yaxis = list(title = "Valor (normalizado)"),
      xaxis = list(title = "Iteración")
    )
  )

fig




# Cargar librerías
library(umap)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(patchwork)
library(plotly)
library(tsne)

data <- read_tsv("/home/rodrigo/Desktop/optimizeCoRa_29sep/optimizeCoRa/Output/OUT_OptCoRa_ATFv1_p01_Fig1_mY_mY.txt")
features <- data[2:8]
labels <- data[9]
umaps = umap(features, n_components=2, random_state = 15)
layout <- umaps[["layout"]]
layout <- data.frame(layout)
final <- cbind(layout, labels)



fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~robustness, type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='species')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 


fig 



library(corrplot)
library(ggplot2)
library(GGally)
library(viridis)






df <- data[c("mU", "mW", "eP", "gU", "gW", "e0", "eM", "robustness")]

df$mU <- log10(df$mU)
df$mW <- log10(df$mW)
df$eP <- log10(df$eP)
df$gU <- log10(df$gU)
df$gW <- log10(df$gW)
df$e0 <- log10(df$e0)
df$eM <- log10(df$eM)

df <- df[df$robustness > 0.8, ]

# Matriz de correlación
cor_matrix <- cor(df)
corrplot(cor_matrix, method = "color", type = "upper", 
         order = "hclust", tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(100))

# Pair plot con correlaciones
ggpairs(df, columns = 1:(ncol(df)),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5, color = df$robustness)))




library(GGally)
library(ggplot2)

# Función personalizada para puntos con límites fijos
custom_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.3, size = 0.5, ...) +
    xlim(-3, 3) +
    ylim(-3, 3)
}

# Pair plot con límites fijos
ggpairs(df, columns = 1:(ncol(df)-1),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = custom_points)) +
  theme_minimal()
