hubercls <- function(r, delta) {
  (1 - r - 0.5 * delta) * (r <= 1 - delta) + 0.5 * (1 - r)^2 / delta * (r <= 1) * (r > 1 - delta)
} 

ercls <- function(r, omega) {
  abs(omega - (r < 0)) * r^2
} 
