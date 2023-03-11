using LinearAlgebra
using Gnuplot

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
L = 0.031

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
LENGTH = 100
STEP = L * 2 / LENGTH
f = zeros(LENGTH+1, LENGTH+1)
for dx in range(0, LENGTH)
  for dy in range(0, LENGTH)
    # Sampled position
    local x = [STEP * (dx - LENGTH/2); STEP * (dy - LENGTH/2)]
    # Apply Hooke's law to each kilobot and sum the contributions
    fv = [0; 0]
    for i in range(1,4)
      diff = KBs[:,i] - x
      dist = norm(diff)
      if dist > 0.0
        dir = diff / dist
        hk = hooke(dist, K, L)
        springf = hk * dir
        fv = fv + springf
      end
    end
    # Store the force magnitude
    f[dx+1,dy+1] = f[dx+1,dy+1] + norm(fv) - FRICTION
  end  
end

#XYS = STEP * (range(0, LENGTH) - ones(LENGTH+1) * (LENGTH/2))
#plt = heatmap(XYS, XYS, f, aspect_ratio=:equal)

function gradient(m, scale::Real, start::Tuple{Int,Int}, stop::Tuple{Int,Int}, length::Tuple{Int,Int})
  xsize = size(m)[1]
  ysize = size(m)[2]
  u = zeros(xsize, ysize)
  v = zeros(xsize, ysize)
  for r in range(1, xsize)
    v[r,2:ysize] = (m[r,2:ysize] - m[r,1:(ysize-1)]) .* scale
  end
  for c in range(1, ysize)
    u[2:xsize,c] = (m[2:xsize,c] - m[1:(xsize-1),c]) .* scale
  end
  xsamples = trunc.(Int, collect(range(start[1], stop[1], length=length[1]))) .+ 1
  ysamples = trunc.(Int, collect(range(start[2], stop[2], length=length[2]))) .+ 1
  u = u[xsamples, ysamples]
  v = v[xsamples, ysamples]
  return (u,v)
end

function meshgrid(start::Int, stop::Int, length::Int)
  dx = ones(length) .* range(start, stop, length=length)'
  dy = dx'
  dxy = permutedims(cat(dx, dy, dims=3), [3,1,2])
  return reshape(dxy, 2, :)'
end

START = 4
STOP = size(f)[1] - 1 - 4
SAMPLES = 24
xy = meshgrid(START, STOP, SAMPLES)
u,v = gradient(f, 5000, (START,START), (STOP,STOP), (SAMPLES,SAMPLES))
u = reshape(u', 1, :)' .* -1
v = reshape(v', 1, :)' .* -1
df = cat(xy, u, v, dims=2)
@gp    "set size square" "set xrange [0:100]" "set yrange [0:100]" f  "with image notitle"
@gp :- df[:,1] df[:,2] df[:,3] df[:,4] "using 1:2:3:4 with vectors notitle ls -1"
save(term="pngcairo fontscale 0.8", output="heatmap.png")
