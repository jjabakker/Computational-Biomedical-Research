# ODE model
csc = function(t, state, parameters){
    with(as.list(c(state,parameters)),{
        p3   = 1 - p1 - p2
        q3   = 1 - q1 - q2
        dcsc = (p1 - p3) * rcsc * CSC - dcsc * CSC
        dpc  = (2 * p3 + p2) * rcsc * CSC + (q1 - q3) * rpc * PC - dpc * PC
        ddc  = (q2 + 2 * q3) * rpc * PC - ddc * DC
        list(c(dcsc, dpc, ddc))
    })
}

# Run growth model and plot results
run_model = function(t, state, parameters){
    library(deSolve)
	# run model
 
    out = ode(func=csc, y=state, parms=parameters, times=times)
	# extract results
    time = out[,1]
    CSC  = out[,2]
    PC   = out[,3]
    DC   = out[,4]
    n    = CSC+PC+DC
		
	# create plots
  # dev.new(width=12, height=6)
		
    par(mfrow=c(2,4))
    plot(time, CSC,main="CSC",type="l",xlab="time",ylab="cells")
    plot(time, PC,main="PC",type="l",xlab="time",ylab="cells")    
    plot(time, DC,main="DC",type="l",xlab="time",ylab="cells")
    plot(time, n,main="population",type="l",xlab="time",ylab="cells")
    plot(time, 100*CSC/n,main="percentage CSC",type="l",xlab="time",ylab="%",ylim=c(0,100))
    plot(time, 100*PC/n,main="percentage PC",type="l",xlab="time",ylab="%",ylim=c(0,100))
    plot(time, 100*DC/n,main="percentage DC",type="l",xlab="time",ylab="%",ylim=c(0,100))

    out
  }
    
    
