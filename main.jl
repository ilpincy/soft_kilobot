include("SoftKilobot.jl")
import .SoftKilobot as SK

##
# Duration in seconds
DURATION = 0.25;

function Done(
  u::Vector{Float64}, # Motion state of the kilobots
  ss::SK.State,       # The SoftKilobot
  t::Integer)         # Time step
  return ((t > DURATION/SK.DELTA_T) || (SK.CoM(u, ss)[1] > SK.ENV_END))
end;

# Initialization
ss, u0 = SK.Init(
  [0.4,1],           # Center of the robot distribution
  135.0,             # Orientation of robot distribution (degrees)
  3,                 # Robots per side
  SK.KILOBOT_MASS,   # Robot mass
  SK.KILOBOT_RADIUS, # Robot radius [m]
  SK.KILOBOT_AXLE,   # Robot wheel axle [m]
  SK.KILOBOT_VMAX,   # Max speed
  0.01,              # Speed damping
  nothing,           # Noise on wheel speed
  0.031,             # Spring rest length
  0.6,               # Spring stiffness
  0.1);              # Spring damping
u = vec(u0);

# Execution
anim = Animation()
i = 0
while (!Done(u, ss, i))
  global i
  tspan = (i * SK.DELTA_T, (i+1) * SK.DELTA_T);
  prob = ODEProblem(SK.Step!, u, tspan, ss);
  sol = solve(prob);
  for j in 1:length(sol.t)
    fr = SK.Plot(
      reshape(sol.u[j], SK.FIELDS, ss.n),
      ss,
      sol.t[j],
      vel=true,
      nbr_springs=false,
      nbr_alg=false);
    frame(anim, fr);
  end
  global u = sol.u[end]
  i += 1
end
gif(anim, "no_noise.gif", fps=31);

# # Initialization
# ss, u0 = SK.Init(
#     [0.4,1],                                                 # Center of the robot distribution
#     135.0,                                                   # Orientation of robot distribution (degrees)
#     1,                                                       # Robots per side
#     SK.KILOBOT_MASS,                                         # Robot mass
#     SK.KILOBOT_RADIUS,                                       # Robot radius [m]
#     SK.KILOBOT_AXLE,                                         # Robot wheel axle [m]
#     SK.KILOBOT_VMAX,                                         # Max speed
#     0.01,                                                    # Speed damping
#     Some(Normal(SK.KILOBOT_VMAX*0.1, SK.KILOBOT_VMAX*0.01)), # Noise on wheel speed
#     0.031,                                                   # Spring rest length
#     0.6,                                                     # Spring stiffness
#     0.1)                                                     # Spring damping
# u = vec(u0)

# # Execution
# i = 0
# @gif while (!Done(u, ss, i))
#     tspan = (i * SK.DELTA_T, (i+1) * SK.DELTA_T);
#     prob = ODEProblem(SK.Step!, u, tspan, ss);
#     sol = solve(prob)
#     for j in 1:length(sol.t)
#         SK.Plot(reshape(sol.u[j], FIELDS, ss.n), sol.t[j], vel=true, nbr_springs=false, nbr_alg=false)
#     end
#     u = sol.u[end]
#     i += 1
# end
