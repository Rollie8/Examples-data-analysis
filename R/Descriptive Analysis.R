devtools::load_all("/Users/alessio/Documents/Data analytics/MyPackageR/R")
library(MyPackageR)
library(plotly)
library(ggplot2)
library(readxl)
library(dbscan)
library(cluster)
library(reshape2)

D5_colors <- read_excel("D5_colors.xlsx")
View(D5_colors)
D5_colors <- D5_colors[-c(1:2),]
D5_colors$R <- as.integer(D5_colors$R)
D5_colors$G <- as.integer(D5_colors$G)
D5_colors$B <- as.integer(D5_colors$B)


# 1  3D CLUSTERING ####
## 1.1 3D HIERARCHICAL ####
dati3D <- D5_colors
str(dati3D)
distances <- dist(dati3D[, c(3:5)], method = "euclidean")


melted_distance <- melt(as.matrix(distances))

# Plot the distances matrix
ggplot(data = melted_distance, aes(Var1, Var2, fill = value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = c("red", "yellow", "green", "cyan", "blue"),
                       c(min(melted_distance$value), max(melted_distance$value)),  # Assicura che la scala vada dal min al max
                       name = "Distance") + # Scala arcobaleno
  theme_minimal() +  # Tema pulito e moderno
  labs(x = "Samples", y = "Samples", fill = "Distance", title = "Distances Matrix") +  # Etichette
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) +  # Mostra tutti i valori da 1 a 30 sull'asse x
  scale_y_continuous(breaks = c(1,5,10,15,20,25,30)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# Ruota i nomi delle colonne

dendrogram_complete <- hclust(distances, method = "complete")
dendrogram_ward <- hclust(distances, method = "ward.D2")

par(mfrow = c(1,2))
plot(dendrogram_complete, main = "Complete linkage")
plot(dendrogram_ward, main = "Ward's method")
par(mfrow = c(1,1))

plot(dendrogram_ward, main = "Ward's method")
grid()
rect.hclust(dendrogram_ward, k = 7, border = c(1:7))
# plot(dendrogram_ward, main = "Ward's method")
# grid()
# rect.hclust(dendrogram_ward, k = 10, border = c(1:7))

# We chosen 7
clustersH <- cutree(tree = dendrogram_ward, k = 7)

# Let's see the real colors of the points
plot_colors3D <- plot_ly(x = dati3D$R, y = dati3D$G, z = dati3D$B, type = "scatter3d", 
                         mode = "markers",  
                         marker = list(
                         size = 12,        # Dimension of markers
                         color = dati3D$`Hex Code`,  # Hex code colors of markers
                         opacity = 0.8),
                         text = dati3D$COLOR) %>%
                         layout(title = list(text = "Plot of colors", y =0.95))

result <- hexToPoints(dataset = dati3D, clusters = clustersH, clust_name = "clusterH")
dati3D <- result$dataset

# Now let's see them clusterized-3D Plot
plot_clusters3D_H <- plot_ly()
for (cluster_id in unique(dati3D$clusterH)) {
  plot_clusters3D_H <- plot_clusters3D_H %>%
    add_trace(x = dati3D$R[dati3D$clusterH == cluster_id], 
              y = dati3D$G[dati3D$clusterH == cluster_id], 
              z = dati3D$B[dati3D$clusterH == cluster_id], 
              type = "scatter3d", 
              mode = "markers", 
              marker = list(
                size = 12,
                opacity = 0.8,
                color = dati3D$hex_centroid_clusterH[dati3D$clusterH == cluster_id],
                showlegend = TRUE
              ),
              text = paste("Cluster: ", dati3D$clusterH[dati3D$clusterH == cluster_id], 
                           "<br>", "Color:", dati3D$COLOR[dati3D$clusterH == cluster_id]),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
              layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors3D, plot_clusters3D_H, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

# PARALLEL COORDINATES PLOT
(plot_parallel3D_H <- plot_ly(
  type = "parcoords",  # Tipo di grafico parallelo
  line = list(
    color = dati3D$clusterH,  # Colore della linea basato sul colore del centroide
    colorscale = "Rainbow",  # Usa una scala di colori, se preferisci puoi personalizzare o usare "RdBu", "Jet", ecc.
    showscale = TRUE  # Mostra la barra dei colori
  ),
  dimensions = list(  # Ogni dimensione è una variabile che vogliamo mostrare
    list(
      label = "R",
      values = dati3D$R
    ),
    list(
      label = "G",
      values = dati3D$G
    ),
    list(
      label = "B",
      values = dati3D$B
    ),
    list(
      label = "Cluster",
      values = dati3D$clusterH
    )
  )
) %>%
  layout(
    title = "Parallel Coordinates Plot of Clusters",
    font = list(size = 12)
  ))

# Plot the distances matrix
dati3D = dati3D[order(dati3D$clusterH),]

counter = 0.5
rect_data = data.frame()

for (cluster in unique(dati3D$clusterH)) {
  cluster_points <- length(dati3D$clusterH[dati3D$clusterH == cluster])
  
  # Crea una riga per il rettangolo corrente
  rect_data <- rbind(rect_data, data.frame(
    xmin = counter,
    xmax = counter + cluster_points,
    ymin = counter,
    ymax = counter + cluster_points
  ))
  counter = counter + cluster_points
}

similarity_matrix = ggplot() +
  geom_tile(data = melted_distance, aes(Var1, Var2, fill = value)) +
  scale_fill_gradientn(colors = c("red", "yellow", "green", "cyan", "blue"),
                       name = "Distance") +
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", size = 0.5) +
  theme_minimal()

similarity_matrix

## 1.2 3D PARTITIONAL ####


# VALIDATION MEASURES -> SSW & Silhouette
# Cluster_validity3D <- data.frame(Clusters = 2:10, 
#                                 SSW = NA, 
#                                 Silhouette_Avg_width = NA)

# for(i in 2:10){
#   model <- kmeans(dati3D[, c(3:5)], centers = i, nstart = 50, iter.max = 20)
#   Cluster_validity3D$SSW[i-1] <- model$tot.withinss
#   sil <- silhouette(x = model$cluster, dist = distances)
#   width_cluster <- aggregate(sil, by = list(sil[, 1]), FUN = mean)
#   Cluster_validity3D$Silhouette_Avg_width[i-1] <- sum(width_cluster[,4] * table(model$cluster))/sum(table(model$cluster))
# }

Cluster_validity3D <- summaryKmeans(dati3D[ ,c(3:5)], 10)

# After the N°cluster = 9 the silhouette starts to level off, we should use 9 clusters
plot(Cluster_validity3D[, 1:2], type = 'b', pch = 20)
plot(Cluster_validity3D[, c(1,3)], type = 'b', pch = 20)

set.seed(15)
kmeans <- kmeans(dati3D[, c(3:5)], centers = 9, nstart = 50, iter.max = 20)
View(kmeans$centers)
dati3D <- hexToPoints(dataset = dati3D, clusters = kmeans$cluster, clust_name = "clusterKmeans")$dataset
dati3D = dati3D[order(dati3D$clusterKmeans),]

# Now let's see them clusterized
plot_clusters3D_Kmeans <- plot_ly()
for (cluster_id in unique(dati3D$clusterKmeans)){
    plot_clusters3D_Kmeans <- plot_clusters3D_Kmeans %>%
    add_trace(x = dati3D$R[dati3D$clusterKmeans == cluster_id], 
              y = dati3D$G[dati3D$clusterKmeans == cluster_id], 
              z = dati3D$B[dati3D$clusterKmeans == cluster_id], 
              type = "scatter3d", 
              mode = "markers", 
              marker = list(
                size = 12,
                opacity = 0.8,
                color = dati3D$hex_centroid_clusterKmeans[dati3D$clusterKmeans == cluster_id],
                showlegend = TRUE
              ),
              text = paste("Cluster: ", dati3D$clusterKmeans[dati3D$clusterKmeans == cluster_id], 
                           "<br>", "Color:", dati3D$COLOR[dati3D$clusterKmeans == cluster_id]),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors3D, plot_clusters3D_Kmeans, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

# Plot the distances matrix
distances <- dist(dati3D[, c(3:5)], method = "euclidean")
melted_distance <- melt(as.matrix(distances))

ggplot(data = melted_distance, aes(Var1, Var2, fill = value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = c("red", "yellow", "green", "cyan", "blue"),
                       c(min(melted_distance$value), max(melted_distance$value)),  # Assicura che la scala vada dal min al max
                       name = "Distance") + # Scala arcobaleno
  theme_minimal() +  # Tema pulito e moderno
  labs(x = "Samples", y = "Samples", fill = "Distance", title = "Distances Matrix") +  # Etichette
  scale_x_continuous(breaks = c(1,5,10,15,20,25,30)) +  # Mostra tutti i valori da 1 a 30 sull'asse x
  scale_y_continuous(breaks = c(1,5,10,15,20,25,30)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))# Ruota i nomi delle colonne

## 1.3 3D DENSITY-BASED ####


# In orange dbscan widget:
# eps = "neighbordhood distance" and minPts = "core point neighbors"

# DBscan has problems in case of high-dimensional data and varying densities
(dbscan_3D<- dbscan(dati3D[, c(3:5)], eps = 80, minPts = 2))

kNNdistplot(distances, k = 2)
grid()

dati3D <- hexToPoints(dataset = dati3D, clusters = dbscan_3D$cluster, clust_name = "clusterDB")$dataset
dati3D = dati3D[order(dati3D$clusterDB),]
# Now let's see them clusterized
plot_clusters3D_DB <- plot_ly()
for (cluster_id in unique(dati3D$clusterDB)){
  plot_clusters3D_DB <- plot_clusters3D_DB %>%
    add_trace(x = dati3D$R[dati3D$clusterDB == cluster_id], 
              y = dati3D$G[dati3D$clusterDB == cluster_id], 
              z = dati3D$B[dati3D$clusterDB == cluster_id], 
              type = "scatter3d", 
              mode = "markers", 
              marker = list(
                size = 12,
                opacity = 0.8,
                color = dati3D$hex_centroid_clusterDB[dati3D$clusterDB == cluster_id],
                showlegend = TRUE
              ),
              text = paste("Cluster: ", dati3D$clusterDB[dati3D$clusterDB == cluster_id], 
                           "<br>", "Color:", dati3D$COLOR[dati3D$clusterDB == cluster_id]),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

subplot(plot_colors3D, plot_clusters3D_DB, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

# 2 2D CLUSTERING ####
plot(D5_colors[, c(3:5)], pch = 20)
## 2.1  R-B CLUSTERING ####
### 2.1.1  R-B HIERARCHICAL ####
dati_rb <- D5_colors

dist_rb <- dist(dati_rb[, c(3,5)], method = "euclidean")
mod_rb <- hclust(dist_rb, method = "ward.D2")

plot(mod_rb)
grid()
# 1 cluster less 
rect.hclust(mod_rb, k = 6, border = c(1:6))

clusters_rbH <- cutree(mod_rb, k = 6)
dati_rb$clusterH <-clusters_rbH

# Let's see the real colors of the points
plot_colors_rb <- plot_ly(x = dati_rb$R, y = dati_rb$B, type = "scatter", mode = "markers",  
                       marker = list(
                         size = ~dati_rb$G/10 + 8,        # Dimension of markers
                         color = dati_rb$`Hex Code`,  # Hex code colors of markers
                         opacity = 0.5
                       ), text = dati_rb$COLOR)

plot_clusters_rb_H <- plot_ly()
for (cluster_id in unique(dati_rb$clusterH)){
  cluster_data <- dati_rb[dati_rb$clusterH == cluster_id,]
  plot_clusters_rb_H <- plot_clusters_rb_H %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$G/10 + 8,
                opacity = 0.8,
                showlegend = TRUE),
              color = I(cluster_data$clusterH),
              text = paste("Cluster: ", cluster_data$clusterH,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rb, plot_clusters_rb_H, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.1.2  R-B PARTITIONAL ####
Cluster_validity_rb <- summaryKmeans(dati_rb[,c(3,5)], cmax = 10)

# After the N°cluster = 9 the silhouette starts to level off, we should use 9 clusters
plot(Cluster_validity_rb[, 1:2], type = 'b', pch = 20)
plot(Cluster_validity_rb[, c(1,3)], type = 'b', pch = 20)

set.seed(15)
kmeans_rb <- kmeans(dati_rb[, c(3,5)], centers = 9, nstart = 50, iter.max = 20)
dati_rb$clusterKmeans <- kmeans_rb$cluster
dati_rb = dati_rb[order(dati_rb$clusterKmeans),]

plot_clusters_rb_Kmeans <- plot_ly()
for (cluster_id in unique(dati_rb$clusterKmeans)){
  cluster_data <- dati_rb[dati_rb$clusterKmeans == cluster_id,]
  plot_clusters_rb_Kmeans <- plot_clusters_rb_Kmeans %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$G/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),
              color = I(cluster_data$clusterKmeans),
              text = paste("Cluster: ", cluster_data$clusterKmeans,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rb, plot_clusters_rb_Kmeans, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.1.3 R-B DENSITY-BASED #### 
# I prefer 60 as eps
# (dbscan_rb <- dbscan(dati_rb[, c(3,5)], eps = 60, minPts = 2))
(dbscan_rb <- dbscan(dati_rb[, c(3,5)], eps = 38, minPts = 2))

# k here = core point neighbors + 1 on orange
# minpts here = core point neighbors + 2 on orange
kNNdistplot(dist_rb, k = 2)
grid()

dati_rb$clusterDB <- dbscan_rb$cluster

cluster_colors = c("black","red", "yellow",
                    "lightgrey", "cyan", "blue",
                    "pink", "purple", "brown")

plot_clusters_rb_DB <- plot_ly()
for (cluster_id in sort(unique(dati_rb$clusterDB))){
  cluster_data <- dati_rb[dati_rb$clusterDB == cluster_id,]
  plot_clusters_rb_DB <- plot_clusters_rb_DB %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                color = cluster_colors[cluster_id+1],
                size = cluster_data$G/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),# the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterDB,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rb, plot_clusters_rb_DB, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

## 2.2  R-G CLUSTERING ####
### 2.2.1  R-G HIERARCHICAL ####
dati_rg <- D5_colors

dist_rg <- dist(dati_rg[, c(3,4)], method = "euclidean")
mod_rg <- hclust(dist_rg, method = "ward.D2")

plot(mod_rg)
rect.hclust(mod_rg, k = 7, border = c(1:7))

clusters_rgH <- cutree(mod_rg, k = 7)
dati_rg$clusterH <- clusters_rgH

# Let's see the real colors of the points
plot_colors_rg <- plot_ly(x = dati_rg$R, y = dati_rg$G, type = "scatter", mode = "markers",  
                          marker = list(
                            size = ~dati_rg$B/10 + 8,        # Dimension of markers
                            color = dati_rg$`Hex Code`,  # Hex code colors of markers
                            opacity = 0.8
                          ), text = dati_rg$COLOR)

plot_clusters_rg_H <- plot_ly()
for (cluster_id in unique(dati_rg$clusterH)){
  cluster_data <- dati_rg[dati_rg$clusterH == cluster_id,]
  plot_clusters_rg_H <- plot_clusters_rg_H %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$G,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$B/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),
              color = I(cluster_data$clusterH),# the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterH,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rg, plot_clusters_rg_H, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.2.2 R-G PARTITIONAL####
Cluster_validity_rg <- summaryKmeans(dati_rg[, c(3,4)], cmax = 10)

# After the N°cluster = 8 the silhouette starts to level off, we should use 8 clusters
plot(Cluster_validity_rg[, 1:2], type = 'b', pch = 20)
plot(Cluster_validity_rg[, c(1,3)], type = 'b', pch = 20)

set.seed(15)
kmeans_rg <- kmeans(dati_rg[, c(3,4)], centers = 8, nstart = 50, iter.max = 20)
dati_rg$clusterKmeans <- kmeans_rg$cluster
dati_rg = dati_rg[order(dati_rg$clusterKmeans),]

plot_clusters_rg_Kmeans <- plot_ly()
for (cluster_id in unique(dati_rg$clusterKmeans)){
  cluster_data <- dati_rg[dati_rg$clusterKmeans == cluster_id,]
  plot_clusters_rg_Kmeans <- plot_clusters_rg_Kmeans %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$G,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$B/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),
              color = I(cluster_data$clusterKmeans),# the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterKmeans,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rg, plot_clusters_rg_Kmeans, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.2.3 R-G DENSITY-BASED ####
(dbscan_rg <- dbscan(dati_rg[, c(3,4)], eps = 38, minPts = 2))

# k here = core point neighbors + 1 on orange
# minpts here = core point neighbors + 2 on orange
kNNdistplot(dist_rg, k = 2)
grid()

dati_rg$clusterDB <- dbscan_rg$cluster
dati_rg = dati_rg[order(dati_rg$clusterDB),]

plot_clusters_rg_DB <- plot_ly()
for (cluster_id in unique(dati_rg$clusterDB)){
  cluster_data <- dati_rg[dati_rg$clusterDB == cluster_id,]
  plot_clusters_rg_DB <- plot_clusters_rg_DB %>%
    add_trace(x = cluster_data$R, 
              y = cluster_data$G,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$B/10 + 8,
                opacity = 0.5,
                showlegend = TRUE,
                color = cluster_colors[cluster_id+ 1]), # the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterDB,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_rg, plot_clusters_rg_DB, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

## 2.3  G-B CLUSTERING ####
### 2.3.1 G-B HIERARCHICAL ####
dati_gb <- D5_colors

dist_gb <- dist(dati_gb[, c(4,5)], method = "euclidean")
mod_gb <- hclust(dist_gb, method = "ward.D2")

plot(mod_gb)
# We should decide among 7 and 8
rect.hclust(mod_gb, k = 8, border = c(1:8))

clusters_gbH <- cutree(mod_gb, k = 8)
dati_gb$clusterH <- clusters_gbH

# Let's see the real colors of the points
plot_colors_gb <- plot_ly(x = dati_gb$G, y = dati_gb$B, type = "scatter", mode = "markers",  
                          marker = list(
                            size = ~dati_gb$R/15 + 10,        # Dimension of markers
                            color = dati_gb$`Hex Code`,  # Hex code colors of markers
                            opacity = 0.8
                          ), text = dati_gb$COLOR)

plot_clusters_gb_H <- plot_ly()
for (cluster_id in unique(dati_gb$clusterH)){
  cluster_data <- dati_gb[dati_gb$clusterH == cluster_id,]
  plot_clusters_gb_H <- plot_clusters_gb_H %>%
    add_trace(x = cluster_data$G, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$R/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),
              color = I(cluster_data$clusterH),# the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterH,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}
# Show both the plots in the same window
subplot(plot_colors_gb, plot_clusters_gb_H, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.3.2 R-G PARTITIONAL ####
Cluster_validity_gb <- summaryKmeans(dati_gb[, c(4,5)], cmax = 10)

plot(Cluster_validity_gb[, 1:2], type = 'b', pch = 20)
plot(Cluster_validity_gb[, c(1,3)], type = 'b', pch = 20)

# Since on hierarchical we have chosen 8, let's choose 7 here
# However it seems that we have the same clusterization of the hierarchical
set.seed(15)
kmeans_gb<- kmeans(dati_gb[, c(4,5)], centers = 9, nstart = 50, iter.max = 20)
dati_gb$clusterKmeans <-kmeans_gb$cluster
dati_gb = dati_gb[order(dati_gb$clusterKmeans),]

plot_clusters_gb_Kmeans <- plot_ly()
for (cluster_id in unique(dati_gb$clusterKmeans)){
  cluster_data <- dati_gb[dati_gb$clusterKmeans == cluster_id,]
  plot_clusters_gb_Kmeans <- plot_clusters_gb_Kmeans %>%
    add_trace(x = cluster_data$G, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$R/10 + 8,
                opacity = 0.5,
                showlegend = TRUE),
              color = I(cluster_data$clusterKmeans),# the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterKmeans,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}
# Show both the plots in the same window
subplot(plot_colors_gb, plot_clusters_gb_Kmeans, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

### 2.3.3  G-B DENSITY-BASED ####
(dbscan_gb <- dbscan(dati_gb[, c(4,5)], eps = 45, minPts = 2))

# k here = core point neighbors + 1 on orange
# minpts here = core point neighbors + 2 on orange
kNNdistplot(dist_gb, k = 2)
grid()

dati_gb$clusterDB <- dbscan_gb$cluster
dati_gb = dati_gb[order(dati_gb$clusterDB),]

plot_clusters_gb_DB <- plot_ly()
for (cluster_id in unique(dati_gb$clusterDB)){
  cluster_data <- dati_gb[dati_gb$clusterDB == cluster_id,]
  plot_clusters_gb_DB <- plot_clusters_gb_DB %>%
    add_trace(x = cluster_data$G, 
              y = cluster_data$B,
              type = "scatter", 
              mode = "markers", 
              marker = list(
                size = cluster_data$R/10 + 8,
                opacity = 0.5,
                showlegend = TRUE,
                color = cluster_colors[cluster_id+1]), # the color of the cl.0 is equal to cl.8
              text = paste("Cluster: ", cluster_data$clusterDB,
                           "<br>", "Color:", cluster_data$COLOR),
              name = paste("Cluster", cluster_id)) %>% # Nome del tracciato basato sul cluster
    layout(title = list(text = "Plot of clusters", y =0.95))  
}

# Show both the plots in the same window
subplot(plot_colors_gb, plot_clusters_gb_DB, margin = 0.05, 
        nrows = 1, shareX = FALSE, shareY = FALSE,
        titleX = TRUE, titleY = TRUE)

