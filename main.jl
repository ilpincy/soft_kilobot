using DifferentialEquations
using Distributions
using Plots

# See https://discourse.julialang.org/t/plots-gr-png-results-in-too-many-open-files/48998/11
ENV["GKSwstype"] = "100"

include("SoftKilobot.jl")
import .SoftKilobot as SK

# Duration in seconds
#DURATION = SK.DELTA_T;
#DURATION = 800
DURATION = 120
#DURATION = 1

# Whether the simulation is done
function Done(
  u::Vector{Float64}, # Motion state of the kilobots
  ss::SK.State,       # The SoftKilobot
  t::Integer)         # Time step
  return ((t > DURATION/SK.DELTA_T) || (SK.CoM(u, ss)[1] > SK.ENV_END))
end;

# Noise
#NOISE = nothing
NOISE = Some(Normal(SK.KILOBOT_WMAX * 0.001, SK.KILOBOT_WMAX * 0.001))

# Initialization
ss, u0 = SK.Init(
  [0.4,1],           # Center of the robot distribution
  135.0,             # Orientation of robot distribution (degrees)
  3,                 # Robots per side
  SK.KILOBOT_MASS,   # Robot mass
  SK.KILOBOT_RADIUS, # Robot radius [m]
  SK.KILOBOT_AXLE,   # Robot wheel axle [m]
  SK.KILOBOT_WMAX,   # Max wheel speed
  SK.KILOBOT_FMAX,   # Max linear force
  SK.KILOBOT_TMAX,   # Max torque
  0.99,              # Speed damping
  NOISE,             # Noise on wheel speed
  0.031,             # Spring rest length
  0.1, #0.6,               # Spring stiffness
  0.1);              # Spring damping
u = vec(u0);
println("Initialized");

# Execution
i = 0
steps = []
while (!Done(u, ss, i))
  global i
  tspan = (i * SK.DELTA_T, (i+1) * SK.DELTA_T);
  prob = ODEProblem(SK.Step!, u, tspan, ss);
  sol = solve(prob);
  push!(steps, sol);
  global u = sol.u[end]
  i += 1
  println("Step #", i);
end
println("Done with execution");

anim = Animation()
for i in 1:SK.STEPS_PER_SEC:length(steps)
  step = steps[i]
  fr = SK.Plot(
    reshape(step.u[end], SK.FIELDS, ss.n),
    ss,
    step.t[end],
    vel=true,
    nbr_springs=false,
    nbr_alg_ds=true,
    nbr_alg_as=false);
  frame(anim, fr);
end
gif(anim, "no_noise.gif", fps=25, loop=0);
println("Done with gif");
