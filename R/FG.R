FG <- function(t, y, p, delay, vinit) {
  
  vlag <- ifelse((t - delay) <= 0, vinit, lagvalue(t - delay)[2])
  g <- function(x) 1/(1 + exp(-4*x))
  du <- p["alpha"] * (-y["u"] + p["q"] * g(vlag) + p["I"]) 
  dv <- p["C"] * (y["w"] + y["v"] - 1/3 * y["v"]^3) + p["gamma"] * y["u"]
  dw <- (p["a"] - y["v"] - p["b"] * y["w"])/p["C"] 
  
  return(list(c(du, dv, dw)))
}

## define parameters
parms <- c(a=0.9, b=0.9, C=2.0, gamma=1.02,
           I=-2.5, alpha=0.01, q=1)
## define integrations times
times <- seq(from=0, to=600, by = 0.5)
## define initial state
init <- c(u=-2.7, v=-0.8, w=1.9)
