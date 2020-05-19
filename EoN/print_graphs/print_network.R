library(network)
print_graph <- function(path, name) {
  test <- read.csv(path)
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
  plot(NET,vertex.col="black", main=name)
}
pdf(file=paste0("~/Philly_Covid/EoN/print_graphs/res.pdf"), width=8, height=8)
print_graph("~/Philly_Covid/EoN/print_graphs/local.csv", "local network")
print_graph("~/Philly_Covid/EoN/print_graphs/expanded_local.csv", "expanded local network")
print_graph("~/Philly_Covid/EoN/print_graphs/O.csv", "original network")
print_graph("~/Philly_Covid/EoN/print_graphs/SD.csv", "original network with social distancing")
print_graph("~/Philly_Covid/EoN/print_graphs/expanded_SD.csv", "original network with social distancing, local expansion")
print_graph("~/Philly_Covid/EoN/print_graphs/add_back_long_SD.csv", "original network with social distancing, 10% long distance edges added back")
print_graph("~/Philly_Covid/EoN/print_graphs/add_back_long_expanded_SD.csv", "original network with social distancing, \n10% long distance edges added back, local expansion")
dev.off()

print_local_ld <- function(path_local, path_ld, name_1, name_2) {
  temp <- read.csv(path_local)
  temp <- as.character(temp[,1])
  local <- matrix(0, ncol=999, nrow=999)
  for (i in 1:length(local)) {
    edge <- temp[i]
    node1 <- as.numeric(strsplit(temp[i], " ")[[1]][1])
    node2 <- as.numeric(strsplit(temp[i], " ")[[1]][2])
    local[node1, node2] <- 1
    local[node2, node1] <- 1
  }
  temp <- read.csv(path_ld)
  temp <- as.character(temp[,1])
  ld <- matrix(0, ncol=999, nrow=999)
  for (i in 1:length(ld)) {
    edge <- temp[i]
    node1 <- as.numeric(strsplit(temp[i], " ")[[1]][1])
    node2 <- as.numeric(strsplit(temp[i], " ")[[1]][2])
    ld[node1, node2] <- 1
    ld[node2, node1] <- 1
  }
  local_net <- network(local, directed=FALSE)
  ld_net <- network(ld, directed=FALSE)
  COORDS <- plot(local_net, vertex.col="black")  #this will plot SD nicely and save coordinates
  plot(local_net,vertex.col='blue', jitter=FALSE, coord=COORDS, main=name_1)
  plot(ld_net, vertex.col='blue', jitter=FALSE, coord=COORDS, main=name_2)
}
pdf(file=paste0("~/Philly_Covid/EoN/print_graphs/networks.pdf"), width=8, height=8)
print_local_ld(path_local="~/Philly_Covid/EoN/print_graphs/local.csv", 
               path_ld="~/Philly_Covid/EoN/print_graphs/ld.csv", 
               name_1="local circles",
               name_2="original: local circles and 100% of long distance connections")
print_local_ld(path_local="~/Philly_Covid/EoN/print_graphs/local.csv", 
               path_ld="~/Philly_Covid/EoN/print_graphs/SD_ld.csv", 
               name_1="local circles",
               name_2="SD: local circles and 10% of long distance connections")
print_local_ld(path_local="~/Philly_Covid/EoN/print_graphs/expanded_local.csv", 
               path_ld="~/Philly_Covid/EoN/print_graphs/SD_ld.csv", 
               name_1="fused local circles",
               name_2="fused circles: fused local circles and 10% of long distance connections")
print_local_ld(path_local="~/Philly_Covid/EoN/print_graphs/local.csv", 
               path_ld="~/Philly_Covid/EoN/print_graphs/ld_20.csv", 
               name_1="local circles",
               name_2="LD added back: local circles and 20% of long distance connections")
print_local_ld(path_local="~/Philly_Covid/EoN/print_graphs/expanded_local.csv", 
               path_ld="~/Philly_Covid/EoN/print_graphs/ld_20.csv", 
               name_1="fused local circles", 
               name_2="fused circles + LD added back:\nfused circles and 19% of long distance connections")
dev.off()
