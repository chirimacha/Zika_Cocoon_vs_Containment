"
This code generates networks
"

globalParams <- list(
  graphSavePath="output"
  ,workingDir="~/Zika_Cocoon_vs_Containment"
)

setwd(globalParams$workingDir)

library(igraph)
library(MASS)

makeGrid <- function(numNodes=200) {
  #width = largest non-trivial divisor of numNodes
  width = seq(sqrt(numNodes))
  width = max(width[numNodes %% width== 0])
  height = numNodes / width
  G <- make_lattice(dimvector=c(height,width), circular=F)
  Gmat <- get.adjacency(G)
  rnd <- as.integer(runif(1, max=100000))
  fpath <- paste0(globalParams$graphSavePath, "/", "grid", numNodes, "_rnd", rnd, ".csv")
  print(fpath)
  write.matrix(Gmat, file=fpath, sep=",")
  return(G)
}

##############################
#               MAIN
##############################
print("Making grids...")
for(i in seq(30)) {
  G<-makeGrid()
}

print(G)
plot(G, vertex.size=2, vertex.color="red")
