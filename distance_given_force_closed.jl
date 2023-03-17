using LinearAlgebra

# The aim of this script is to calculate the maximum distance that can
# be reached by a robot given a certain force.

#
# Hooke's Law
#
# Here we assume that extension occurs in the positive x
hooke(x::Real, k::Real, l::Real) = k * (x - L)
# Elastic constant
K = 0.6
# Rest length (meters)
L = 0.031

#
# Position of the "anchor" Kilobots
#
KBs = [L 0 -L  0;
       0 L  0 -L]

#
# Force magnitude (assuming no friction)
#
function force_polar(x::Real, theta::Real)
  D1 = sqrt(L^2 + x^2 - 2*L*x * cos(theta))
  D2 = sqrt(L^2 + x^2 - 2*L*x * sin(theta))
  D3 = sqrt(L^2 + x^2 + 2*L*x * cos(theta))
  D4 = sqrt(L^2 + x^2 + 2*L*x * sin(theta))
  A = x * ((1-L/D1) + (1-L/D2) + (1-L/D3) + (1-L/D4))
  B = L^2 * (1/D3 - 1/D1)
  C = L^2 * (1/D4 - 1/D2)
  return K * sqrt(A^2 + B^2 + C^2 - 2*A*(B*cos(theta) + C*sin(theta)))
end

function force_cartesian(rho::Real, theta::Real)
  x = [rho*cos(theta); rho*sin(theta)]
  d1 = x - KBs[:,1]
  d2 = x - KBs[:,2]
  d3 = x - KBs[:,3]
  d4 = x - KBs[:,4]
  n1 = sqrt(dot(d1, d1))
  n2 = sqrt(dot(d2, d2))
  n3 = sqrt(dot(d3, d3))
  n4 = sqrt(dot(d4, d4))
  e1 = d1 / n1
  e2 = d2 / n2
  e3 = d3 / n3
  e4 = d4 / n4
  f = (n1-L)*e1 + (n2-L)*e2 + (n3-L)*e3 + (n4-L)*e4
  return K * sqrt(dot(f,f))
end

#
# Parameters
#
# Target forces to consider
FORCE_MIN = 0.0054 # Newton
FORCE_MAX = 0.0167 # Newton
# Gravitational constant
G = 6.67430e−11
# Kilobot mass
M = 0.01 # kg
# Friction coefficient
MU = 1.0
# Friction for any kilobot
FRICTION = MU * M * G # Newton

#
# Look for maximum distance
#
# From the plot we know it's in the direction theta = pi/4 + i*pi/2,
# so it's enough to test theta = pi/4 and the rest is symmetric.
#
function max_dist(target_f::Real)
  for x in range(0.0, L*2, step=0.0001)
    if force_polar(x, pi/4) - FRICTION > target_f
      return x
    end
  end
end

#
# Let's go
#
print_dist(f) = println("Max distance for f = ", f, "N and friction = ", MU, " is ", max_dist(f) * 100, "cm")
print_dist(FORCE_MIN)
print_dist(FORCE_MAX)

#
# Validation
#
validatep(x, theta) = println("Polar force for (", x, ",", theta, ") = ", force_polar(x, theta))
validatec(x, theta) = println("Cartesian force for (", x, ",", theta, ") = ", force_cartesian(x, theta))
validatep(0.028, pi/4)
validatep(0.028, 0)
validatec(0.028, pi/4)
validatec(0.028, 0)
