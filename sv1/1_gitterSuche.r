data = read.csv("R/sv1/soilrespiration.csv")
temp = unlist(data["temp"], use.names = 0)
resp = unlist(data["resp"], use.names = 0)
log.resp = log(resp)

# log.resp ~ beta1 + beta2 * temp

# erste Aufgabe: erstelle ein Gitter und suche die besten Werte

gitterSuche = function(beta1s, beta2s, errorFunc){
  bestError = Inf
  bestBeta1 = 0
  bestBeta2 = 0
  
  for (beta1 in beta1s) {
    for (beta2 in beta2s) {
      error = sum(errorFunc(beta1, beta2))
      if(error < bestError){
        bestError = error
        bestBeta1 = beta1
        bestBeta2 = beta2
      }
    }
  }
  
  return(c(bestBeta1, bestBeta2, bestError))
  
}

result = gitterSuche(
  seq(-1, 0, length.out = 1001),
  seq(0, 0.1, length.out = 1001),
  function(beta1, beta2) ((beta1 + beta2 * temp) - log.resp)^2
)

print(c('best Error:', result[1], 'beta1:', result[1], 'beta2:', result[2]))