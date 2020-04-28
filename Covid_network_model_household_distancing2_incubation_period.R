#===================================================================================
#  Lots of R packages for graphs
#	these are described in a CRAN task view here:
#	http://cran.r-project.org/web/views/gR.html
#	-we will use the 'network' package
#===================================================================================

#install.packages("network")
library(network)
#help("network-package")


#===================================================================================
#	-Simulated population parameters
#===================================================================================

L<-num_obs<-1000  #number of households and individuals
S<-rep(1,L)
E<-rep(0,L)
I<-rep(0,L)
R<-rep(0,L)


#===================================================================================
#	import neighborhood connections
#	-NEIGH: is a matrix of connections
#	-NEIGH_NET is the same saved as a network
#================================================================================

reset_network <- function() {
  NEIGH<-as.matrix(read.csv('~/Philly_Covid/example_network.csv', header=F))
  
  # force to be symmetric ----------------
  #NEIGH <- Matrix::forceSymmetric(as.matrix(NEIGH),uplo="L")
  NEIGH <- as.matrix(NEIGH)
  NEIGH_NET<-network(NEIGH, directed=FALSE)
  NEIGH_NET <<- NEIGH_NET
  NEIGH <<- NEIGH
  rowSums(NEIGH)
  hist(rowSums(NEIGH)) # look at degree distribution -----------------------------
  
  #===================================================================================
  #	create a network of long distance connections
  #===================================================================================
  
  connectivity<-.01
  LD<-matrix(rbinom(L^2,1,connectivity),L)
  for( i in 1:L)
  {
    LD[i,]<-LD[,i]
  }
  
  diag(LD)<-0
  LD <<- LD
  NET<-network(LD, directed=FALSE)
  rowSums(LD)
  hist(rowSums(LD)) # look at degree distribution --------------------------------
  mean(rowSums(LD))
  
  #===================================================================================
  #	create a single connection network out of NEIGH and LD
  #===================================================================================
  ALL_NET<-(NEIGH+LD)>0  #the unity of connections from long and neighborhood connections
  ALL_NET<-ALL_NET*1       #so that NET is displayed in 0 and 1s
  rowSums(ALL_NET)
  hist(rowSums(ALL_NET)) # look at degree distribution
  
  #ALL_NET<-ALL_NET*0+1   #for homogeneous mixing
  #LD<-ALL_NET
  
  NET<<-network(ALL_NET, directed=FALSE)
  ALL_NET<<-ALL_NET
}
reset_network()


#===================================================================================
#	Plotting functions
#	-COORDS lets it draw network once and saves the coordinates of each node
#	-COLS a vector to define colors in the plot 
#	-visualize: plots network, coloring the I nodes red
#===================================================================================


COORDS<-plot(NEIGH_NET,vertex.col="black")
COLS<-rep(1,L)

visualize_net<-function()
{
  COLS[I==1]<-"red"
  COLS[E==1]<-"orange"
  COLS[S==1]<-"black"
  COLS[R==1]<-"white"
  plot(NET,vertex.col=COLS, jitter=FALSE, coord=COORDS)
}


visualize_nodes<-function()
{
  COLS[I==1]<-"red"
  COLS[E==1]<-"orange"
  COLS[S==1]<-"black"
  COLS[R==1]<-"white"
  lines(COORDS,col=COLS,pch=19,type="p") # much faster than points()
}

#===================================================================================
#	-Set up Simulation 
#	-S,E,I,R vectors
#===================================================================================
setup<-function()
{
  S<<-rep(1,L)
  E<<-rep(0,L)
  I<<-rep(0,L)
  R<<-rep(0,L)
  recoverday<<-rep(NA,L)
}

#===================================================================================
#	Key functions of stochastic simulation
#	-infect changes I from 0 to 1, also changes E and S to 0
#===================================================================================

infect<-function(node,day=i)
{
  I[node]<<-1
  E[node]<<-0
  S[node]<<-0
  cases<<-which(I==1)
  recoverday[node]<<-day+duration
}

expose<-function()
{
  if(length(cases)==0){
    E<<-E*0
    return()
  }
  
  # function exposes all neighbors according to the ALL_NET object
  # maybe check if people are being reexposed after SD. They were infected once.
  
  if(length(cases)==1)
  {
    toexpose<-which(ALL_NET[,cases]==1 & S)
  }
  
  if(length(cases)>1)
  {
    toexpose<-which(rowSums(ALL_NET[,cases])>=1 & S)
  }
  
  E[toexpose]<<-1
  E[cases]<<-0
  S[toexpose]<<-0
  
  #make sure no recovereds are reexpose
  E[R]<<-0
  #R[toexpose]<<-0
}

recover<-function(node)
{
  I[node]<<-0
  E[node]<<-0
  S[node]<<-0
  R[node]<<-1
  cases<<-which(I==1)
  
  #unexpose all connected to recovered nodes
  if(length(node)==1) # typo -----
  {
    tounexpose<-which(ALL_NET[,node]==1)
  }
  
  if(length(node)>1)
  {
    tounexpose<-which(rowSums(ALL_NET[,node])>=1)
  }
  
  E[tounexpose]<<-0
  
  #reexpose those who had more than one infectious contact
  expose()
}



#===================================================================================
#	Simulation Setup
#	-reps: define how many time steps
#	-setup()  set S to 1, I, E, R to 0
#	-PREV: an empty vector to keep track of prevalence at each timestep
#	-par sets up a plot window. par(ask=TRUE) requires a <Enter> between plots (slows sown our movie)
#	-b : probaility infection given exposure / time step (ie rate) 
#	-duration : number of days infectious
#   -b is calculated from R0 and mean duration of infection
#	-Assign index cases
#	-infect it with infect()
#	-expose its neighbors with expose()
#===================================================================================

#define the length in days of the simulation
reps<-100

#set time to day 1
i<-1

#an empty vector to store the prevalence over the course of the simulation
PREV<-rep(0,reps)


#duration of infectiousness
duration<-11

#probability of infection in exposed node per time step
#R0 = duration * b * average number of contacts
R0<-3

b <- R0 / (duration * mean(rowSums(ALL_NET)))

#b<-.05
#sets the whole populations to susceptible
setup()

#set how many index cases and draw them randomly
index<-index_case<-sample(L,20)

#infect the index cases
infect(index_case)

#expose those nodes connected to the index cases
expose()

#plot the starting conditions
#par(ask=FALSE)
#visualize_net()

#===================================================================================
#	Simulation loop
#	-reps: define how many time steps
#   -social distancing is implemented at dist_day. Set to above reps for no distancing
#   -each individual decreases samples which contacts to maintain based on their distance_factor
#   - distance factor is drawn from a beta distribution
#   -The LD long distance contact matrix is reduced by binomial sampling if each contact is to be maintained. 
#   -If one individual eliminates a contact, it is also eliminated for the other
#
#===================================================================================

# incubation period
incubation_period <- 5

sim_loop <- function(dist_day, input_beta) {
  reset_network()
  print(paste0("mean num. neighbors before SD: ", mean(rowSums(ALL_NET))))
  print(paste0("num. edges: ", sum(rowSums(ALL_NET)) / 2))
  infectday <- rep(0, L)
  reps <- 200
  i <- 1
  PREV<-rep(0,reps)
  setup()
  index<-index_case<-sample(L,20)
  infect(index_case)
  expose()
  
  
  for(i in 2:reps)
  {
    #social distancing at 50 cases
    #if(length(cases)>=50) {
    
    if(i==dist_day) {
      print("SD enacted")
      LD2<-LD
      distance_factor<-rbeta(L,input_beta,5)
      for(j in 1:L)
      {
        LD2[j,]<-rbinom(n=L,size=1,prob=distance_factor[j]*LD[j,])
      }
      #set to no contact if either node reduces contact
      for(j in 1:L){
        for(k in 1:L)
        {
          LD2[j,k]<-LD2[j,k]*LD2[k,j]
        }
      }

      ALL_NET<-NEIGH+LD2>0
      
      for (i in nrow(test2)) {
        for (j in ncol(test2)) {
          if (test2[i,j] < 0) {
            print("error")
          }
        }
      }
      #*
      # LD2<-ALL_NET
      # distance_factor<-rbeta(L,input_beta,5)
      # for(j in 1:L)
      # {
      #   LD2[j,]<-rbinom(n=L,size=1,prob=distance_factor[j]*LD[j,])
      # }
      # #set to no contact if either node reduces contact
      # for(j in 1:L){
      #   for(k in 1:L)
      #   {
      #     LD2[j,k]<-LD2[j,k]*LD2[k,j]
      #   }
      # }
      # ALL_NET <- LD2>0
      #*
      
      ALL_NET<<-ALL_NET*1       #so that NET is displayed in 0 and 1s
      NET<<-network(ALL_NET, directed=FALSE)
      #visualize_net()
      print(paste0("Dist. day: ", i))
      print(paste0("mean num. neighbors after SD: ", mean(rowSums(ALL_NET))))
      print(paste0("num. edges: ", sum(rowSums(ALL_NET)) / 2))
    }
    
    torecover<-which(recoverday==i)
    if (length(torecover)>=1) {
      recover(torecover)
    }
    
    if(length(cases)>=1){
      random<-runif(L)
      risk<-E*b		#multiplying the risk by whether or not expsoed
      toinfect<-which(random<risk*E)
      
      infectday[toinfect] = i + incubation_period
      
      if (sum(infectday==i)>=1) {
        infect(which(infectday==i), i)
        #infect(toinfect, i)
        expose()
      }
    }
    
    #visualize_nodes()
    PREV[i]<-(sum(I))/L
    CUMULPREV <- sum(I+R)/L
  }
  return(PREV)
}

#===================================================================================
#	- plot infections over time
#	- vertical line at day of social distancing
#===================================================================================
intervention_times <- sort(rep(c(2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 105), times=5), decreasing=F)
prevs <- list()
prev_dex <- 1
f_size <- list()
f_size_dex <- 1
for (i_time in intervention_times) {
  prevs[[prev_dex]] <- sim_loop(dist_day = i_time)
  prev_dex <- prev_dex + 1
  f_size[[f_size_dex]] <- length(which(R == 1))
  f_size_dex <- f_size_dex + 1
  print(paste0("Final size: ", length(which(R == 1))))
}
prevs_df <- do.call(rbind, prevs)
f_size_df <- do.call(rbind, f_size)

pdf(file=paste0("~/Philly_Covid/plots/interventions_ip_", incubation_period, ".pdf"),width=9, height=12)
par(mfrow=c(4,3))
plot_prev <- function(PREV, first, sd_day) {
  if (first) {
    plot(PREV,typ="l", ylab="Prevalence",xlab="Time step",ylim=c(0,1), main=paste0("SD enacted: ", sd_day))
  } else {
    lines(PREV)
  }
}
for (row_num in 1:nrow(prevs_df)) {
  if (row_num %% 5 == 1) {
    plot_prev(PREV=prevs_df[row_num,], first=T, sd_day=intervention_times[row_num])
  } else {
    plot_prev(PREV=prevs_df[row_num,], first=F, sd_day=intervention_times[row_num])
  }
  abline(v=intervention_times[row_num], col="blue")
}
dev.off()

#===================================================================================
#	- plot the effect of social distancing
#===================================================================================
intervention_times <- sort(rep(c(seq(from=3, to=15, by=3), 105), times=1), decreasing=F)
betas <- rep(5, times=10)
#intervention_times <- sort(rep(c(10, 105), times=5))
#betas <- c(5)
prevs <- list()
prev_dex <- 1
f_size <- list()
f_size_dex <- 1
for (i_time in intervention_times) {
  for (beta in betas) {
    prevs[[prev_dex]] <- sim_loop(dist_day = i_time, input_beta=beta)
    prev_dex <- prev_dex + 1
    f_size[[f_size_dex]] <- length(which(R == 1))
    f_size_dex <- f_size_dex + 1
    print(paste0("Final size: ", length(which(R == 1))))
  }
}
prevs_df <- do.call(rbind, prevs)
f_size_df <- do.call(rbind, f_size)

pdf(file=paste0("~/Desktop/network", ".pdf"),width=12, height=8)
par(mfrow=c(2,3))
index <- 1
for (i_time in intervention_times) {
  for (beta in betas) {
    if (index %% 10 == 1) {
      plot(prevs_df[index,],typ="l", ylab="I/N",xlab="Time step",ylim=c(0,0.3),xlim=c(0, 150), main=paste0("Network model, i_time=", i_time))
      abline(v=i_time, col="green")
    } else {
      lines(prevs_df[index,])
    }
    index <- index + 1
  }
}
dev.off()


# SIR
library(deSolve)
sir <- function(time, state, parameters) 
{
  if (round(time) >= 10) {
    with(
      as.list(c(state, parameters)), 
      {
        dS <- -beta2 * S * I
        dI <- beta2 * S * I - gamma * I
        dR <- gamma * I
        return(list(c(dS, dI, dR)))
      }
    )
  } else {
    with(
      as.list(c(state, parameters)), 
      {
        dS <- -beta * S * I
        dI <- beta * S * I - gamma * I
        dR <- gamma * I
        return(list(c(dS, dI, dR)))
      }
    )
  }
}

input_beta <- 0.0002727273
x <- sort(seq(input_beta/10, input_beta/1, by=input_beta/10), decreasing=T)
first <- T
for (in_beta in x[1]) {
  init <- c(S=1000, I=1, R=0)
  # beta solved as 0.0002727273
  parameters <- c(beta = 0.0002727273, beta2 = in_beta, gamma = 1/11)
  times <- seq(0, 300, by = 1)
  out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
  S<-out$S
  I<-out$I
  R<-out$R
  RESULTS<-data.frame(out$S,out$I,out$R)
  if (first) {
    plot(times,I/1000, ylim=c(0,0.3),xlim=c(0,120),  xlab = "t (days)", ylab = "I/N", main = "Comparison of simple mass-action model to simple network model", typ="l",col="red")
    #abline(v=10, col="green")
    #lines(times,I / 1000, col="red", lty="longdash")
    first <- F
  } else {
    #lines(times,I / 1000, col="red", lty="longdash")
    #lines(times,I/1000, col="black")
  }
}
for (index in 6:10) {
  lines(prevs_df[index,], col="black")
}
legend(80, 0.3, legend=c("Mass-action model", "Network model"),
       col=c("red", "black"), lty=1, cex=0.8)

