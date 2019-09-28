
require(deSolve)
require(GillespieSSA)
require(ggplot2)

one_hit_ODE <- function(times, parms, InitialCellCount){
  
  one_hitmodel<- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dfree <- -free * k1 + conjugate*k2
      dconjugate <- free*k1 - conjugate * (k2 + k3)
      ddead <-conjugate * k3
      return(list(c(dfree, dconjugate, ddead)))
    })
  }
  pars_one_hit <- parms
  yini_one_hit <- c(free = InitialCellCount, conjugate = 0, dead=0)
  
  out <- ode(yini_one_hit, times, one_hitmodel, pars_one_hit)
  
  output_df<-data.frame(out[,"time"],out[,"dead"])
  names(output_df)<-c("time","deaths")
  
  return(output_df)
}

two_hit_ODE <- function(times, parms, InitialCellCount){
  
  two_hitmodel<- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      
      dfree <- -free * k1 + conjugate * k2
      dconjugate <- free * k1-conjugate * (k2 + k3)
      dconjugate_damaged <- conjugate*k3 + free_damaged * k5 - conjugate_damaged * (k4+k6)
      dfree_damaged <- k4 * conjugate_damaged - k5 * free_damaged
      ddead <- k6 * conjugate_damaged
      return(list(c(dfree, dconjugate, dfree_damaged,dconjugate_damaged,ddead)))
    })
  }
  pars_two_hit <- parms
  yini_two_hit <- c(free = InitialCellCount, conjugate = 0, free_damaged = 0, conjugate_damaged = 0, dead=0)
  
  out <- ode(yini_two_hit, times, two_hitmodel, pars_two_hit)
  
  output_df<-data.frame(out[,"time"],
                        out[,"dead"])
  names(output_df)<-c("time","deaths")
  
  return(output_df)
  
}

one_hit_gillespie <- function(parms,InitialCellCount){
  
  ## initial state vector
  x0 <- c(free = InitialCellCount, conjugate = 0, dead = 0)
  
  ## Governing Equations
  a <- c("k1 * free", "k2 * conjugate", "k3 * conjugate")
  
  ## matrix specifying which interactions take place, and what direction.
  ## i.e. dY1 = (+1)*c1*Y1 (-1)*C2*Y1*Y2
  nu <- matrix(c(-1,+1,0,+1,-1,-1,0,0,+1), 
               nrow = 3, byrow = TRUE)
  
  out <- ssa(x0, a, nu, parms, tf = 12, method = "D",simName = "3 State Gillespie")
  
  output_df<-data.frame(out$data[,1], out$data[,"dead"])
  
  names(output_df)<-c("time", "deaths")
  
  return(output_df)
}

two_hit_gillespie <- function(parms,InitialCellCount){
  
  ## rate constants (should look like this)
  ## parms <- c(k1=0.2, k2=1, k3=0.2, k4=0.2,k5=0.2,k6=0.2)
  
  ## initial state vector
  x0 <- c(free = InitialCellCount, 
          conjugate = 0, 
          free_damaged = 0, 
          conjugate_damaged = 0,
          dead = 0)
  
  ## Governing Equations
  a <- c("k1 * free",
         "k2 * conjugate",
         "k3 * conjugate",
         "k4 * conjugate_damaged",
         "k5 * free_damaged",
         "k6 * conjugate_damaged")
  
  ## matrix specifying which interactions take place, and what direction.
  ## i.e. dY1 = (+1)*c1*Y1 (-1)*C2*Y1*Y2
  nu <- matrix(c(-1,+1,0,0,0,0,+1,-1,-1,0,0,0,0,0,0,+1,-1,0,0,0,+1,-1,+1,-1,0,0,0,0,0,+1),
               nrow=5,
               byrow=TRUE)
  
  
  out <- ssa(x0, a, nu, parms, tf = 12,  method ="D", simName ="5 State Gillespie")
  
  output_df <- data.frame(out$data[,1],
                          out$data[,"dead"])
  
  names(output_df)<-c("time","deaths")
  
  return(output_df)
  
}

eventFinder <- function(ODEdata){
  
  total_deaths<-floor(tail(ODEdata$deaths,1))
  
  if(total_deaths==0) return(0)
  
  indices<-rep(0,total_deaths)
  
  for (i in 1:total_deaths){
    index<- min(which(ODEdata$deaths >= i))
    indices[i]<-index
  }
  
  ODEdata <- ODEdata[indices,]
  
  ODEdata$deaths <- floor(ODEdata$deaths)
  
  ODEdata <- unique(ODEdata)
  
  return(ODEdata)
  
}

runSimulation <-function(simtype, parms, InitialCellCount, times=NULL, SD = NULL) {
  
  if (simtype != "one_hit_gillespie" && 
      simtype != "two_hit_gillespie" && 
      simtype != "one_hit_ODE" && 
      simtype != "two_hit_ODE") {
    stop("enter a valid simtype")
  }
  
  if (length(parms)!=3) stop ("check parms")
  
  if(!is.null(SD)) {
    
    randomcellcount<-0
    
    while (randomcellcount<=0){
      randomcellcount<-round(rnorm(1,InitialCellCount,SD))
    }
    
    InitialCellCount <- randomcellcount
  }
  
  if (simtype == "one_hit_gillespie") {
    names(parms) <- c("k1", "k2", "k3")
    out <- one_hit_gillespie(parms,InitialCellCount)
  }
  
  if (simtype == "two_hit_gillespie") {
    parms<-c(parms,parms[2],parms[1],parms[3])
    names(parms) <- c("k1", "k2", "k3", "k4", "k5", "k6")
    out <- two_hit_gillespie(parms,InitialCellCount)
  }
  
  if (simtype == "one_hit_ODE") {
    if(is.null(times)) stop("set 'times' argument")
    names(parms) <- c("k1","k2","k3")
    out <- one_hit_ODE(times,parms,InitialCellCount)
  }
  
  if (simtype == "two_hit_ODE") {
    if(is.null(times)) stop("set 'times' argument")
    parms<-c(parms, parms[2], parms[1], parms[3])
    names(parms) <- c("k1", "k2", "k3", "k4", "k5", "k6")
    out<- two_hit_ODE(times,parms,InitialCellCount)
  }
  out<-eventFinder(out)
  
  return(out)
  
}

multiWell<-function(n,simtype,parms,InitialCellCount,times=NULL, SD = NULL){
  
  
  mylist <- replicate(n, 
                      runSimulation(simtype,parms,InitialCellCount,times,SD), 
                      simplify=FALSE)
  return(mylist)
}

wellPlot<-function(out){
  
  counter1 <-0
  counter2 <-0
  
  filtered_output <- list()
  list_len <- 0
  
  while (list_len<10) {
    counter1 <- counter1 + 1
    if (counter1 > length(out)) break
    
    if (typeof(out[[counter1]]) == "list") {
      counter2 <- counter2 + 1
      filtered_output[[counter2]]<-out[[counter1]]
    }
    
    list_len <- length(filtered_output)
  }
  
  out <- filtered_output
  
  ggplot_df<-data.frame(time=numeric(0),well=numeric(0))
  
  k=0;
  
  for (i in 1:length(out)) {
    for (j in 1:length(out[[i]][,1])) {
      k=k+1
      ggplot_df[k,]<-c(out[[i]]$time[j],i)
      a <- 1
    }
  }
  
  p <- ggplot(ggplot_df, aes(time,well))
  p <- p+geom_point()
  p <- p+scale_y_continuous(breaks=seq(1,length(out),1))
  p <- p+scale_x_continuous(breaks=seq(0,12,2),limits=c(0,12))
  
  print(p)
}

barchartplotter <- function(a){
  
  list_length <- length(a)
  b <- rep(0, list_length)
  
  for (i in 1:list_length) {
    if (typeof(a[[i]]) == "list") {
      b[i]<-tail(a[[i]]$deaths,1)
    }
    else {
      b[i]<-0
    } 
  } 
  
  p <- ggplot(data.frame(b), aes(x=factor(b))) + 
    geom_bar(fill="lightgreen", color="grey50") +
    xlab("deaths")
  
  print(p)

}

