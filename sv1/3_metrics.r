
samples = 101

# Original, Quadrate
print(gitterSuche(seq(-1, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2)
                    ((beta1 + beta2 * temp) - log.resp) ^ 2))

# absolut
print(gitterSuche(seq(-1, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2)
                    abs((beta1 + beta2 * temp) - log.resp)))

# 4. Potenz
print(gitterSuche(seq(-1, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2)
                    ((beta1 + beta2 * temp) - log.resp) ^ 4))

# orthogonale Abstände zur Ausgleichsgeraden
# src: https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line
dist2d = function(p, l0, l1) {
  v1 = l0 - l1
  v2 = p - l0
  m = cbind(v1, v2)
  return(abs(det(m)) / sqrt(sum(v1 * v1)))
}

print(gitterSuche(seq(-1, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2) {
                    sum = 0
                    l0 = c(0, beta1)
                    l1 = c(1, beta1 + beta2)
                    for (i in 1:length(temp)) {
                      p = c(temp[i], log.resp[i])
                      sum = sum + dist2d(p, l0, l1)
                    }
                    return(sum)
                  }))

# Custom, 6. Potenz
print(gitterSuche(seq(-1, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2)
                    ((beta1 + beta2 * temp) - log.resp) ^ 6))

dist2dManhatten = function(p, l0, l1) {
  d = l1-l0
  # die Manhattendistanz von einem Punkt zu einer Linie
  # sei definiert als die Entfernung zwischen dem Punkt und den
  # Punkten, die die gleiche x- bzw y-Koordinate wie p haben.
  # die x-Koordinate vom Punkt mit gleichem y ist:
  pointXWithSameY = l0[1] + d[1] * (p[2]-l0[2])/d[2]
  # die y-Koordinate vom Punkt mit gleichem x ist:
  pointYWithSameX = l0[2] + d[2] * (p[1]-l0[1])/d[1]
  return(abs(p[1]-pointXWithSameY) + abs(p[2]-pointYWithSameX))
}

# Custom, Manhatten-Distanz
print(gitterSuche(seq(-2, 0, length.out = samples),
                  seq(0, 0.1, length.out = samples),
                  function(beta1, beta2) {
                    sum = 0
                    l0 = c(0, beta1)
                    l1 = c(1, beta1 + beta2)
                    for (i in 1:length(temp)) {
                      p = c(temp[i], log.resp[i])
                      sum = sum + dist2dManhatten(p, l0, l1)
                    }
                    return(sum)
                  }))