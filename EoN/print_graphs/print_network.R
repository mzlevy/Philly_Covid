library(network)

# Local --------------------------------------------------------------------
test <- read.csv("~/Philly_Covid/EoN/print_graphs/local.edgelist.csv")
test <- as.character(test[,1])
local <- matrix(0, ncol=999, nrow=999)
for (i in 1:length(test)) {
  edge <- test[i]
  node1 <- as.numeric(strsplit(test[i], " ")[[1]][1])
  node2 <- as.numeric(strsplit(test[i], " ")[[1]][2])
  local[node1, node2] <- 1
  local[node2, node1] <- 1
}
NET<-network(local, directed=FALSE)
plot(NET,vertex.col="black", main="local")

test <- read.csv("~/Desktop/e_local.edgelist.csv")
test <- as.character(test[,1])
local <- matrix(0, ncol=999, nrow=999)
for (i in 1:length(test)) {
  edge <- test[i]
  node1 <- as.numeric(strsplit(test[i], " ")[[1]][1])
  node2 <- as.numeric(strsplit(test[i], " ")[[1]][2])
  local[node1, node2] <- 1
  local[node2, node1] <- 1
}
NET<-network(local, directed=FALSE)
plot(NET,vertex.col="black", main="e_local")

# Full with LD ------------------------------------------------------------------
test <- read.csv("~/Philly_Covid/EoN/print_graphs/O.edgelist.csv")
test <- as.character(test[,1])
local <- matrix(0, ncol=999, nrow=999)
for (i in 1:length(test)) {
  edge <- test[i]
  node1 <- as.numeric(strsplit(test[i], " ")[[1]][1])
  node2 <- as.numeric(strsplit(test[i], " ")[[1]][2])
  local[node1, node2] <- 1
  local[node2, node1] <- 1
}
NET<-network(local, directed=FALSE)
plot(NET,vertex.col="black", main="O")

test <- read.csv("~/Philly_Covid/EoN/print_graphs/SD.edgelist.csv")
test <- as.character(test[,1])
local <- matrix(0, ncol=999, nrow=999)
for (i in 1:length(test)) {
  edge <- test[i]
  node1 <- as.numeric(strsplit(test[i], " ")[[1]][1])
  node2 <- as.numeric(strsplit(test[i], " ")[[1]][2])
  local[node1, node2] <- 1
  local[node2, node1] <- 1
}
NET<-network(local, directed=FALSE)
plot(NET,vertex.col="black", main="SD")

test <- read.csv("~/Philly_Covid/EoN/print_graphs/e_SD.edgelist.csv")
test <- as.character(test[,1])
local <- matrix(0, ncol=999, nrow=999)
for (i in 1:length(test)) {
  edge <- test[i]
  node1 <- as.numeric(strsplit(test[i], " ")[[1]][1])
  node2 <- as.numeric(strsplit(test[i], " ")[[1]][2])
  local[node1, node2] <- 1
  local[node2, node1] <- 1
}
NET<-network(local, directed=FALSE)
plot(NET,vertex.col="black", main="e_SD")

