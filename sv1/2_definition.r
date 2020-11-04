X = matrix(c(rep(1, length(temp)), temp), length(temp), 2)
betaDach = matrix(c(bestBeta1, bestBeta2), 2, 1)

X_betaDach = X %*% betaDach
XtX = t(X) %*% X
invXtX = solve(XtX) # inverse