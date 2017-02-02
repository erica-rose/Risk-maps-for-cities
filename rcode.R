data <- read.csv("sim_data.csv")
# X, Y: coordinates
# block.code: block indicator
# y: infected or not

############
## Plot
plot(data$X, data$Y, pch=16, cex=0.5)
points(data$X[data$y==1], data$Y[data$y==1], col="red", pch=16, cex=1)

#################################
##create origin for every block##
#################################
N <- dim(data)[1]
N.block <- length(unique(data$block.code)) #93 blocks

data$originX <- NA
data$originY <- NA

for(i in unique(data$block.code)){
  
  data$originX[which(data$block.code==i)] <- median(data$X[which(data$block.code==i)])
  data$originY[which(data$block.code==i)] <- median(data$Y[which(data$block.code==i)]) 
}
###################################################
##create x,y matrix for each house on block level##
###################################################
data$Xdiff <- NA
data$Ydiff <- NA

data$Xdiff <- data$X - data$originX
data$Ydiff <- data$Y - data$originY

###############################
##scale x,y of origin matrix ##
###############################
scale.dim.sim <- function(S, data){
  Xscale <- data$originX*S+data$Xdiff
  Yscale <- data$originY*S+data$Ydiff
  return(list(Xscale, Yscale))
}



library(INLA)
S.seq<- seq(1,5,by=0.1)
sim.result.matrix<- matrix(NA, nrow=length(S.seq), ncol=8)

y<- data$y

#################################
#########recover values##########
#################################

for(j in 1:length(S.seq)){
  
  new.coords <- scale.dim.sim(S.seq[j], data=data)

  #define coordinate matrix
  coords <- cbind(new.coords[[1]], new.coords[[2]])
  
  #create mesh 
  mesh <- inla.mesh.2d(coords, max.edge=c(0.1*S.seq[j],0.1*S.seq[j]),cutoff=0.01)
  
  plot(mesh)
  points(new.coords[[1]], new.coords[[2]],pch=16,cex=0.2,col="blue")
  
  A <- inla.spde.make.A(mesh=mesh,loc=coords)
  
  spde <- inla.spde2.matern(mesh, alpha=2)
  
  mesh.index <- inla.spde.make.index(name='spatial', n.spde=spde$n.spde)
  
  stack =
    inla.stack(data=list(y=y.sim[,i]),
               A=list(A),
               effects=
                 list(c(mesh.index,list(Intercept=1))),
               tag="est")
  
  
  formula <- y ~ -1 + Intercept + f(spatial, model=spde)
  
  
  result.sim <- inla(formula,
                     data = inla.stack.data(stack, spde=spde),
                     family = c("binomial"),verbose=TRUE,
                     control.predictor = list(A=inla.stack.A(stack), compute=TRUE,link=1),
                     control.inla = list(reordering = "metis"))
  
  theta.mean <- summary(result.sim)$hyperpar$mean
  theta.sd <- summary(result.sim)$hyperpar$sd
  intercept.mean <- summary(result.sim)$fixed[1]
  intercept.sd <- summary(result.sim)$fixed[2]
  lik <- summary(result.sim)$mlik[2]
  sim.result.matrix[j,] <- cbind(theta.mean[1], theta.mean[2],theta.sd[1], theta.sd[2],intercept.mean, intercept.sd, lik, S.seq[j])
  print(j)
}

sim.result.matrix
S.seq[which.max(sim.result.matrix[,7])]
### Likelihood plot
plot(sim.result.matrix[,8], sim.result.matrix[,7])



