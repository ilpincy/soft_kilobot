using LinearAlgebra
using Plots
gr()

# This file produces a heatmap of the force necessary by a kilobot to extend
# the springs at a certain distance off the center
#
# F = -dE / dx
#
# E = energy
# x = variable of interest

#
# Hooke's Law
#
# Here we assume that extension occurs in the positive x
hooke(x::Real, k::Real, l::Real) = k * (L - x)
# Elastic constant
K = 0.6
# Rest length
L = 3.1

#
# Friction
#
# Gravitational constant
G = 6.67430e−11
# Kilobot mass
M = 0.01
# Friction coefficient
MU = 0.5
# Friction for any kilobot
FRICTION = MU * M * G

# Position of the "anchor" Kilobots
KBs = [L 0 -L  0;
       0 L  0 -L]

# Sampling of the space
LENGTH = 200
STEP = L * 2 / LENGTH
f = zeros(LENGTH+1, LENGTH+1)
for dx in range(0, LENGTH)
  for dy in range(0, LENGTH)
    # Sampled position
    x = [STEP * (dx - LENGTH/2); STEP * (dy - LENGTH/2)]
    # Apply Hooke's law to each kilobot and sum the contributions
    fv = [0; 0]
    for i in range(1,4)
      diff = KBs[:,i] - x
      dist = norm(diff)
      if dist > 0.0
        dir = diff / dist
        hk = hooke(dist, K, L)
        springf = hk * dir
        frictionf = [0.0; 0.0]
        if abs(hk) > 0.0
          frictionf = -springf * (FRICTION / hk)
        end
        fv = fv + springf + frictionf
      end
    end
    # Store the force magnitude
    f[dx+1,dy+1] = f[dx+1,dy+1] + norm(fv)
  end  
end

XYS = STEP * (range(0, LENGTH) - ones(LENGTH+1) * (LENGTH/2))
display(heatmap(XYS, XYS, f, aspect_ratio=:equal))
