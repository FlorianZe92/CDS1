
implizitesEulerVerfahren = function(yAbleitungsGleichung, y0, t0, tEnd, h) {
  
  # y(t+h)-y(t)
  # ----------- = y'
  #      h
  
  yn = y0
  times = seq(t0, tEnd, h)
  values = rep(0, length(times))
  maxError = 1e-3 * h
  
  for (index in 1:length(times)) {
    tn = times[index]
    
    errorFunction = function(yi) {
      return(yAbleitungsGleichung(yi) - (yi - yn) / h)
    }
    
    yn_p1 = findeNullstelle(yn,
                            errorFunction,
                            createDerivative(errorFunction, h),
                            maxError)
    yn = yn_p1
    
    values[index] = yn
    
  }
  
  plot(times, values)
  
  
}

# sucht die Nullstelle von f per Newton-Verfahren
findeNullstelle = function(x0, fVonX, fAbleitungVonX, maxError) {
  error = maxError * 2
  maxSteps = 5000
  stepCount = 0
  xn = x0
  while (error > maxError && stepCount < maxSteps) {
    xn_p1 = newtonSchritt(xn, fVonX, fAbleitungVonX)
    error = abs(xn_p1 - xn)
    xn = xn_p1
    stepCount = stepCount + 1
  }
  return(xn)
}

# x1
newtonSchritt = function(x0, fVonX, fAbleitungVonX) {
  return(x0 - fVonX(x0) / fAbleitungVonX(x0))
}

createDerivative = function(fVonX, h) {
  return(function(x)
    (fVonX(x + h) - fVonX(x)) / h)
}

# test: Pendel, weil wir wissen, was rauskommt?
# y' = -c * x
implizitesEulerVerfahren(function(x)
  sqrt(max(0, 1-x*x)), 0.5, 0, 1, 1e-2)
