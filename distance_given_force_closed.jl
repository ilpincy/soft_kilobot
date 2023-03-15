using LinearAlgebra

# The aim of this script is to calculate the maximum distance that can
# be reached by a robot given a certain force.

#
# Hooke's Law
#
# Here we assume that extension occurs in the positive x
hooke(x::Real, k::Real, l::Real) = k * (x - L)
# Elastic constant
K = 0.2
# Rest length (meters)
L = 0.031
L2 = L^2

#
# Position of the "anchor" Kilobots
#
KBs = [L 0 -L  0;
       0 L  0 -L]

#
# Force magnitude (assuming no friction)
#
function force(x::Real, theta::Real)
  L2x2 = L^2 + x^2
  TwoLx = 2*L*x
  costheta = cos(theta)
  sintheta = sin(theta)
  d1 = sqrt(L2x2 - TwoLx * costheta)
  d2 = sqrt(L2x2 + TwoLx * sintheta)
  d3 = sqrt(L2x2 + TwoLx * costheta)
  d4 = sqrt(L2x2 - TwoLx * sintheta)
  return hooke(d1, K, L) + hooke(d2, K, L) + hooke(d3, K, L) + hooke(d4, K, L)
end

#
# Parameters
#
FORCE_MIN = 0.0054 # Newton
FORCE_MAX = 0.0167 # Newton

#
# Look for maximum distance
#
# From the plot we know it's in the direction theta = pi/4 + i*pi/2,
# so it's enough to test theta = pi/2 and the rest is symmetric.
#
function max_dist(target_f::Real)
  for x in range(0.0, L*2, step=0.0001)
    if force(x, pi/4) > target_f
      return x
    end
  end
end

#
# Let's go
#
println("Max distance for f = ", FORCE_MIN, " is ", max_dist(FORCE_MIN))
println("Max distance for f = ", FORCE_MAX, " is ", max_dist(FORCE_MAX))
