
implizitesEulerVerfahren2 = function(yAbleitungsGleichung, y0, t0, tEnd, h) {
  
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
      return(yAbleitungsGleichung(yi,tn) - (yi - yn) / h)
    }
    
    yn_p1 = findeNullstelle(yn,
                            errorFunction,
                            createDerivative(errorFunction, h),
                            maxError)
    yn = yn_p1
    
    values[index] = yn
    
  }
  
  plot(times, values)
  print(values[length(values)])
  
}

# test: Pendel, weil wir wissen, was rauskommt?
# y' = -c * x
# implizitesEulerVerfahren2(function(y,t) sqrt(max(0, 1-y*y)), 0.5, 0, 1, 1e-2)

# funktioniert <3
implizitesEulerVerfahren2(function(y,t) 2*t*exp(y)/(1+t*t), 0, -2, 0, 1e-3)
