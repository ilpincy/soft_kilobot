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
hooke(x::Real, k::Real, l::Real) = k * (x - L)
# Elastic constant
K = 0.2
# Rest length
L = 0.031

#
# Position of the "anchor" Kilobots
#
KBs = [L 0 -L  0;
       0 L  0 -L]

X_SAMPLES = 10
X_STEP = L / X_SAMPLES
XS = range(1, X_SAMPLES) * X_STEP
THETA_SAMPLES = 100
THETA_STEP = 2*pi / THETA_SAMPLES
THETAS =
  (range(0, THETA_SAMPLES) - ones(THETA_SAMPLES+1) * (THETA_SAMPLES/2))  *
  THETA_STEP
f = zeros(THETA_SAMPLES+1, X_SAMPLES)
for dx in range(1, X_SAMPLES)
  x = XS[dx]
  for dtheta in range(1, THETA_SAMPLES+1)
    theta = THETAS[dtheta]
    L2x2 = L*L + x*x
    TwoLx = 2*L*x
    costheta = cos(theta)
    sintheta = sin(theta)
    d1 = sqrt(L2x2 - TwoLx*costheta)
    d2 = sqrt(L2x2 + TwoLx*sintheta)
    d3 = sqrt(L2x2 + TwoLx*costheta)
    d4 = sqrt(L2x2 - TwoLx*sintheta)
    f[dtheta,dx] =
      hooke(d1, K, L) +
      hooke(d2, K, L) +
      hooke(d3, K, L) +
      hooke(d4, K, L)
  end
end

XLABELS = string.(transpose(XS))

plt = plot(THETAS, f, labels=XLABELS, xlabel="angle [rad]", ylabel="force [N]")
savefig(plt, "plot.pdf")
display(plt)
