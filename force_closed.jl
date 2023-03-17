using LinearAlgebra
using Plots
using Printf
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
K = 0.6
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
THETA_SAMPLES = 10000
THETA_STEP = 2*pi / THETA_SAMPLES
THETAS =
  (range(0, THETA_SAMPLES) - ones(THETA_SAMPLES+1) * (THETA_SAMPLES/2))  *
  THETA_STEP
f = zeros(THETA_SAMPLES+1, X_SAMPLES)
for dx in range(1, X_SAMPLES)
  x = XS[dx]
  for dtheta in range(1, THETA_SAMPLES+1)
    theta = THETAS[dtheta]
    D1 = sqrt(L^2 + x^2 - 2*L*x * cos(theta))
    D2 = sqrt(L^2 + x^2 - 2*L*x * sin(theta))
    D3 = sqrt(L^2 + x^2 + 2*L*x * cos(theta))
    D4 = sqrt(L^2 + x^2 + 2*L*x * sin(theta))
    A = x * ((1-L/D1) + (1-L/D2) + (1-L/D3) + (1-L/D4))
    B = L^2 * (1/D3 - 1/D1)
    C = L^2 * (1/D4 - 1/D2)
    f[dtheta,dx] =
      K * sqrt(A^2 + B^2 + C^2 - 2*A*(B*cos(theta) + C*sin(theta)))
  end
end

XLABELS = [ "ρ=$(@sprintf "%.2f" x) cm" for x in transpose(XS) .* 100 ]
XTICKS = (
  [-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi],
  ["-π", "-3π/4", "-π/2", "-π/4", 0, "π/4", "π/2", "3π/4", "π"]
)
plt = plot(THETAS, f, labels=XLABELS, xlabel="angle θ [rad]", ylabel="force [N]", xticks=XTICKS)
savefig(plt, "plot.pdf")
display(plt)
