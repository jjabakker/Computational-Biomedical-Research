# This function uses trace, det and eigen values to determine the type of equilibrium 



are_same <- function(a,b) {
  isTRUE(all.equal(a,b))
}


test_stability <- function (e_values, trace, det) {
  ev1 <- e_values[[1]]
  ev2 <- e_values[[2]]
  
  if (det <= 0) {
    return('Saddle')
  }
  else if (isTRUE(all.equal(det, 0))) {
    return('Centre')
  }
  else if ( trace < 0 ) {   # Stable
    if (class(ev1) == "complex") {
      return ('Stable Spiral')
    }
    else {
      return ('Stable Node')
    }
  }
  else if (trace > 0) {  # Unstable
    if (class(ev1) == "complex") {
      return ('Unstable Spiral')
    }
    else {
      return ('Unstable Node')
    }
  }
  else { # You should not arrive here really.....
    return('Centre 1')
  }
}

# This function uses just the eigen values to determine the type of equilibrium 

test_stability_jacobian <- function (Jacobian) {
  
  ev <- eigen(Jacobian)
  
  ev1 <- ev$values[[1]]
  ev2 <- ev$values[[2]]
  
  ev1
  ev2
  
  if (are_same(Re(ev1),0)  && are_same(Re(ev1),0)) {
    return("Centre")   
  }
  else if (Re(ev1) < 0 && Re(ev2) < 0) { #stable
    if (class(ev1) == "complex") {
      return ('Stable Spiral')
    }
    else {
      return ('Stable Node')
    }
  }
  else if (Re(ev1) > 0 && Re(ev2) > 0) { #unstable
    if (class(ev1) == "complex") {
      return ('Unstable Spiral')
    }
    else {
      return ('Unstable Node')
    }
  }
  else {
    return("Saddle")
  }
}


# Simple function to determine the grace of a 2 by 2 Jacobian
tr <- function (jacobian) {
  return(jacobian[1,1] + jacobian[2,2])
}