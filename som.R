############################################################################################################
#
#
# Script to load the shared library of SOM_TRAIN. SOM_TRAIN is a program written in F2003 and by using the 
# iso_binding_c module it is possible to compile functions that can be called from C. These functions
# can be called in R by using the ".C" directive and integrated into a full fledged R function.
# Another alternative would be to use Rcpp to see if the Fortran code can be called from C++ and 
# then integrated into R
#
############################################################################################################
train_som <- function(x,nx,ny,nepoch,alpha,grid_type,distance_type,neigh_type,toroidal) {
  nvar <- ncol(x)
  npat <- nrow(x)
  # Define listst with options
  grid_types <- list(rectangular=0,hexagonal=1)
  distance_types <- list(euclidean=0)
  neigh_types <- list(bubble=0,gaussian=1)
  # 
  grid_option <- grid_types[grid_type][[1]]
  distance_option <- distance_types[distance_type][[1]]
  neigh_option <- neigh_types[neigh_type][[1]]
  
  # Load dynamic library
  dyn.load("self_organized_map_utilities.so")
  # Call to C (Long Live C)
  pp <- matrix(1.0,nrow=(nx*ny),ncol=nvar)
  distortion <- vector("numeric", length = nepoch) #, ncol = 1)
  u_matrix <- matrix(1.0, nrow = (2*nx-1), ncol = (2*ny-1))
  coords <- matrix(0.0, nrow = (nx*ny), ncol = 3)
  number_patterns <- matrix(0L, nrow = nx, ncol = ny)
  node_index <- matrix(0L, nrow = npat, ncol = 3)
  retvals <- .C("train_som_", x = as.matrix(x), nvar = as.integer(nvar), npat = as.integer(npat),  nx = as.integer(nx), 
  ny = as.integer(ny), nepoch = as.integer(nepoch), alpha = as.numeric(alpha),   grid_type = as.integer(grid_option), 
  distance_type = as.integer(distance_option), neigh_type = as.integer(neigh_option), 
  toroidal = as.integer(toroidal), prot = pp, distortion = distortion, u_matrix = u_matrix, coords = coords,
  number_patterns = number_patterns, node_index = node_index)  
  return(retvals)
}

predict_som <- function(prot,nx,ny,new_pat){
  npat <- nrow(new_pat)
  nvar <- ncol(new_pat)
  node_index <- matrix(0L,nrow = npat, ncol = 3)
  #
  dyn.load("self_organized_map_utilities.so")
  #
  retvals <- .C("predict_som_", prot = as.matrix(prot), nx = as.integer(nx), ny = as.integer(ny), 
  new_pat = as.matrix(new_pat), npat = as.integer(npat), nvar = as.integer(nvar), node_index = node_index)
  return(retvals)
}

data(iris)
x_mn <- apply(iris[,1:4],2,min)
x_mx <- apply(iris[,1:4],2,max)
x_rn <- x_mx-x_mn
print(x_rn)
x <- as.matrix(iris[,1:4])
x1 <- scale(x,center = x_mn, scale = x_rn)


#print(x1)
set.seed(12345)
v<-train_som(x1,nx=8,ny=8,nepoch=100,alpha=0.2,grid_type="hexagonal",
             distance_type= "euclidean", neigh_type = "gaussian", toroidal = 1)

# print(v$prot)
# print(v$distortion)
# print(v$u_matrix)
# print(v$coords)
# print(v$node_index)
# print(v$number_patterns)

v1<-predict_som(v$prot,nx = 8,ny = 8,x1)
#print(v1$node_index)

# nx<-8
# ny<-8
# nz<-1
# nvar <- 4
# dist <- vector("numeric",(nx*ny))
# for(ipat in 1:1){
#   dist_min <- 1e7
#   for(inode in 1:(nx*ny)) {
#     dist[inode] <- sum((x1[ipat,]-v1$prot[inode])**2)/nvar
#     if(dist[inode] < dist_min){
#       dist_min <- dist[inode]
#       i_hit <- inode
#     }
#   }  
# }
# print(i_hit)
# index_ <- i_hit
#   cz<-min(1+round((index_-1)/(nx*ny)),nz);
#   cy<-min(1+round((index_-1-(cz-1)*nx*ny)/nx-1),ny);
#   cx<-min(index_-(cz-1)*nx*ny-(cy-1)*nx,nx);
# print(c(cx,cy,cz))
# print(dist)