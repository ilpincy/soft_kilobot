module SoftKilobot

########################################

using LinearAlgebra
using Maybe
using Plots

########################################

export Alg_Idx
export Alg_P
export Alg_NbrsDs
export Alg_NbrsAs
export Init
export Alg_Dist
export Alg_Robot
export Alg
export CalcForces
export DiffSteering
export UpdateState
export Step!
export CoM
export Plot
export PlotTrajectories
export ENV_SIDE
export ENV_START
export ENV_END
export DELTA_T
export VEL_X
export VEL_Y
export VEL_W
export POS_X
export POS_Y
export ORNT
export FIELDS
export LWHEEL
export RWHEEL
export KILOBOT_MASS
export KILOBOT_RADIUS
export KILOBOT_AXLE
export KILOBOT_VMAX

########################################

# Environment params
ENV_SIDE  = 2.0;
ENV_START = 0.4;
ENV_END   = 1.6;

# Algorithm parameters
ALG_DIST_159_1      = 0.01; # 1st deformation test
ALG_DIST_159_2      = 0.04; #0.07; # 2nd deformation test
ALG_ANGL_1          = 85.0 * pi / 180.0
ALG_ANGL_9          = 95.0 * pi / 180.0
ALG_ANGL_5_1        = 95.0 * pi / 180.0
ALG_ANGL_5_2        = 85.0 * pi / 180.0
ALG_DIST_234678_1_1 = 0.11; # 1st deformation test
ALG_DIST_234678_1_2 = 0.10; # 1st deformation test
ALG_DIST_234678_2   = 0.07; # 2nd deformation test
ALG_ANGL_234678_1   = 85.0 * pi / 180.0
ALG_ANGL_234678_2   = 95.0 * pi / 180.0

# Control frequency
STEPS_PER_SEC = 31
DELTA_T = 1.0 / STEPS_PER_SEC;

# Motion state field
VEL_X   = 0;
VEL_Y   = 1;
VEL_W   = 2;
POS_X   = 3;
POS_Y   = 4;
ORNT    = 5;
FIELDS  = 6;

# Wheels
LWHEEL = 1;
RWHEEL = 2;

# Robot parameters
KILOBOT_MASS   = 0.024;
KILOBOT_RADIUS = 0.02;
KILOBOT_AXLE   = 0.025;
KILOBOT_WMAX   = 0.015; # Max wheel speed
KILOBOT_FMAX   = 0.5;   # Max linear force
KILOBOT_TMAX   = 0.5;   # Max torque

########################################

"""
The SoftKilobot state information.
"""
mutable struct State
  ### Constant attributes
  N_SIDE::Int                    # Robots per side
  MASS::Float64                  # Robot mass [kg]
  RADIUS::Float64                # Robot radius [m]
  AXLE::Float64                  # Robot wheel axle [m]
  W_MAX::Float64                 # Max wheel speed
  F_MAX::Float64                 # Max linear force
  T_MAX::Float64                 # Max torque
  V_DAMP::Float64                # Speed damping
  S_REST::Float64                # Spring rest length [m]
  S_STIFF::Float64               # Spring stiffness
  S_DAMP::Float64                # Spring damping
  ### Calculated attributes
  n::Int                         # Number of robots
  I::Float64                     # Moment of inertia
  nbr_springs::Array{Array{Int}} # Neighbors due to spring linkage
  p_alg::Array{Int}              # P index according to algorithm
  nbr_alg_ds::Array{Array{Int}}  # Neighbors to calculate the d's, see Fig. 8
  nbr_alg_as::Array{Array{Int}}  # Neighbors to calculate the a's, see Fig. 8

  vs::Matrix{Float64}            # Linear and rotational velocity in the world frame
  w_bias::Matrix{Float64}        # Bias on wheel speed
  w_speeds::Matrix{Float64}      # Wheel speeds
  step_t::Array{Float64}         # Time for the next control step
  State(n_side,mass,radius,axle,w_max,f_max,t_max,v_damp,s_rest,s_stiff,s_damp) = (
    y = new(n_side,mass,radius,axle,w_max,f_max,t_max,v_damp,s_rest,s_stiff,s_damp);
    y.n = n_side^2;
    y.I = 0.5 * y.MASS * y.RADIUS^2;
    y.nbr_springs = [ Int[] for k in 1:y.n ];
    y.p_alg = zeros(Int64, y.n);
    y.nbr_alg_ds = [ Int[] for k in 1:y.n ];
    y.nbr_alg_as = [ Int[] for k in 1:y.n ];
    y.vs = zeros(3, y.n);
    y.w_bias = zeros(2, y.n);
    y.w_speeds = fill(w_max, (2, y.n));
    y.step_t = zeros(y.n);
    return y)
end;

########################################

"""
Get a robot's index given its coordinates in the lattice.
"""
function Alg_Idx(
  i::Int, j::Int, # The coordinates of the current robot
  n_side::Int)    # The number of robots on the side
  (j-1) * n_side + i;
end;

"""
The P indices of the robots, as defined in Fig. 2 of the paper.
"""
function Alg_P(
  i::Int, j::Int, # The coordinates of the current robot
  n_side::Int)    # The number of robots on the side
  # Front robot
  if i == 1 && j == 1; return 1; end
  # Tail robot
  if i == n_side && j == n_side; return 9; end
  # Inner robots
  if i > 1 && i < n_side && j > 1 && j < n_side; return 5; end
  # Corner robots
  if i == n_side && j == 1;      return 3; end
  if i == 1      && j == n_side; return 7; end
  # Edge robots
  if j == 1;      return 2; end
  if j == n_side; return 8; end
  if i == 1;      return 4; end
  if i == n_side; return 6; end
end;

########################################

"""
The indices of the neighbors to calculate the d's, as defined in Fig. 8 of the paper.
"""
function Alg_NbrsDs(
  i::Int, j::Int, # The coordinates of the current robot
  n_side::Int)    # The number of robots on the side
  p = Alg_P(i,j,n_side)
  # Front robot
  if p == 1
    return [Alg_Idx(i+1,j,n_side),Alg_Idx(i,j+1,n_side)]
  end
  # Tail robot
  if p == 9
    return [Alg_Idx(i,j-1,n_side),Alg_Idx(i-1,j,n_side)]
  end
  # Inner robots
  if p == 5
    return [Alg_Idx(i+1,j-1,n_side),Alg_Idx(i-1,j+1,n_side)]
  end
  # Corner robots
  if p == 3
    return [Alg_Idx(i-1,j+1,n_side),Alg_Idx(i-1,j,n_side)]
  end
  if p == 7
    return [Alg_Idx(i+1,j-1,n_side),Alg_Idx(i,j-1,n_side)]
  end
  # Edge robots
  if p == 2
    return [Alg_Idx(i-1,j+1,n_side),Alg_Idx(i-1,j,n_side)]
  end
  if p == 8
    return [Alg_Idx(i+1,j-1,n_side),Alg_Idx(i,j-1,n_side),Alg_Idx(i-1,j,n_side)]
  end
  if p == 4
    return [Alg_Idx(i+1,j-1,n_side),Alg_Idx(i,j-1,n_side)]
  end
  if p == 6
    return [Alg_Idx(i-1,j+1,n_side),Alg_Idx(i,j-1,n_side),Alg_Idx(i-1,j,n_side)]
  end
end;

"""
The indices of the neighbors to calculate the a's, as defined in Fig. 8 of the paper.
"""
function Alg_NbrsAs(
  i::Int, j::Int, # The coordinates of the current robot
  n_side::Int)    # The number of robots on the side
  p = Alg_P(i,j,n_side)
  # Front robot
  if p == 1
    return [Alg_Idx(i+1, j, n_side), Alg_Idx(i, j+1, n_side)]
  end
  # Tail robot
  if p == 9
    return [Alg_Idx(i-1, j, n_side), Alg_Idx(i, j-1, n_side)]
  end
  # Inner robots
  if p == 5
    return [Alg_Idx(i+1, j, n_side), Alg_Idx(i, j+1, n_side), Alg_Idx(i-1, j, n_side), Alg_Idx(i, j-1, n_side)]
  end
  # North east robots
  if p == 2 || p == 3 || p == 6
    return [Alg_Idx(i-1, j, n_side), Alg_Idx(i, j+1, n_side)]
  end
  # South west robots
  if p == 4 || p == 7 || p == 8
    return [Alg_Idx(i+1, j, n_side), Alg_Idx(i, j-1, n_side)]
  end
end;

########################################

"""
SoftKilobot initialization function.
1. Places the robots
2. Creates the spring connections
3. Calculates the neighbors according to the algorithm
4. Initializes the wheel speed biases
Returns a tuple (State, Matrix{Float64}).
"""
function Init(
  cntr::Vector{Float64}, # Center of the robot distribution [m,m]
  ornt::Float64,         # Orientation of robot distribution [deg]
  n_side::Int,           # Robots per side
  mass::Float64,         # Robot mass [kg]
  radius::Float64,       # Robot radius [m]
  axle::Float64,         # Robot wheel axle [m]
  w_max::Float64,        # Max wheel speed [m/s]
  f_max::Float64,        # Max force [N]
  t_max::Float64,        # Max torque [Nm]
  v_damp::Float64,       # Speed damping
  w_noise,               # Noise on wheel speed
  s_rest::Float64,       # Spring rest length
  s_stiff::Float64,      # Spring stiffness
  s_damp::Float64)       # Spring damping
  # Create empty spring system
  ss = State(n_side,mass,radius,axle,w_max,f_max,t_max,v_damp,s_rest,s_stiff,s_damp)
  # Add spring links
  for k in 1:ss.n
    # Horizontal
    kk = (k - 1) % n_side
    if(kk > 0)        push!(ss.nbr_springs[k],k-1) end
    if(kk < n_side-1) push!(ss.nbr_springs[k],k+1) end
    # Vertical
    kk = floor(Int, (k - 1) / n_side)
    if(kk > 0)        push!(ss.nbr_springs[k],k-n_side) end
    if(kk < n_side-1) push!(ss.nbr_springs[k],k+n_side) end
  end
  # Add algorithm neighbors
  for i in 1:n_side
    for j in 1:n_side
      ss.p_alg[Alg_Idx(i, j, n_side)] = Alg_P(i, j, n_side);
      ss.nbr_alg_ds[Alg_Idx(i, j, n_side)] = Alg_NbrsDs(i, j, n_side);
      ss.nbr_alg_as[Alg_Idx(i, j, n_side)] = Alg_NbrsAs(i, j, n_side);
    end
  end
  # Calculate rotation matrix
  ornt *= pi / 180.0
  rm = [cos(ornt) -sin(ornt); sin(ornt) cos(ornt)]
  # Add robots, the ids grow with positive x and negative y
  #    +x ->
  # -y 1 2 3
  #  | 4 5 6
  #  v 7 8 9
  x = zeros(FIELDS, ss.n)
  for i in 0:(n_side-1)
    for j in 0:(n_side-1)
      # Robot id in the lattice
      k = i + n_side * j + 1
      # World position in the lattice
      px = rm * (s_rest * [i, j]) + cntr
      # Robot state
      x[:, k] = [0.0 0.0 0.0 px[1] px[2] 0.0]
      # Wheel noise
      if(w_noise != nothing)
        @? ss.w_bias[LWHEEL, k] = rand(w_noise)
        @? ss.w_bias[RWHEEL, k] = rand(w_noise)
      end
    end
  end
  return ss, x
end;

########################################

function Alg_GetPos(
  k::Int,             # The index of the current robot
  u::Vector{Float64}, # The motion state of the kilobots
  ss::State)          # The SoftKilobot
  i = (k-1) * FIELDS + 1
  return u[POS_X + i:POS_Y + i]
end

"""
Calculates the distances as defined in Fig. 8 of the paper.
"""
function Alg_Dist(
  k::Int,             # The index of the current robot
  u::Vector{Float64}, # The motion state of the kilobots
  ss::State)          # The SoftKilobot
  # Position vector of this robot
  x = Alg_GetPos(k, u, ss)
  # Make list of distances
  ds = []
  for nbr in ss.nbr_alg_ds[k]
    push!(ds, norm(x - Alg_GetPos(nbr, u, ss)))
  end
  return ds
end;

"""
Calculates the angle given three points.
The angle is referred wrt p.
"""
function Alg_CalcAngle(p, n1, n2)
  v1 = n1 - p
  v2 = n2 - p
  return acos(dot(v1, v2) / (norm(v1) * norm(v2)))
end

"""
Calculates the angles as defined in Fig. 8 of the paper.
"""
function Alg_Angles(
  k::Int,             # The index of the current robot
  u::Vector{Float64}, # The motion state of the kilobots
  ss::State)          # The SoftKilobot
  # Position vector of this robot
  x = Alg_GetPos(k, u, ss)
  # Make list of angles
  as = []
  if length(ss.nbr_alg_as[k]) == 2
    # Position vectors of the neighbors
    nx1 = Alg_GetPos(ss.nbr_alg_as[k][1], u, ss)
    nx2 = Alg_GetPos(ss.nbr_alg_as[k][2], u, ss)
    # Calculate angle
    push!(as, Alg_CalcAngle(x, nx1, nx2))
  else
    # Position vectors of the neighbors
    nx1 = Alg_GetPos(ss.nbr_alg_as[k][1], u, ss)
    nx2 = Alg_GetPos(ss.nbr_alg_as[k][2], u, ss)
    # Calculate angle
    push!(as, Alg_CalcAngle(x, nx1, nx2))
    # Position vectors of the neighbors
    nx1 = Alg_GetPos(ss.nbr_alg_as[k][3], u, ss)
    nx2 = Alg_GetPos(ss.nbr_alg_as[k][4], u, ss)
    # Calculate angle
    push!(as, Alg_CalcAngle(x, nx1, nx2))
  end
  return as
end;

########################################

"""
The SoftKilobot motion control algorithm.
"""
function Alg(
  u::Vector{Float64}, # The motion state of the kilobots
  ss::State,          # The SoftKilobot
  t::Float64)         # The current time
  # By default, everybody keeps their old wheel speeds
  # For each robot
  #   If it's time to update the speeds, do it
  for k in 1:ss.n
    if t > ss.step_t[k]
      # Move forward by default
      ss.w_speeds[:,k] = [KILOBOT_WMAX; KILOBOT_WMAX];
      ss.step_t[k] = t + 0.5;
      dsts = Alg_Dist(k, u, ss);
      angls = Alg_Angles(k, u, ss);
      # Position-based corrections
      if ss.p_alg[k] in [1,5,9]
        # Head/tail/interior robot
        d12 = dsts[1] - dsts[2];
        if(d12 > ALG_DIST_159_1)
          # turn right
          ss.w_speeds[:,k] = [KILOBOT_WMAX; -KILOBOT_WMAX];
        elseif(d12 < -ALG_DIST_159_1)
          # turn left
          ss.w_speeds[:,k] = [-KILOBOT_WMAX; KILOBOT_WMAX];
        elseif(
          ((ss.p_alg[k] == 1) && (angls[1] < ALG_ANGL_1) && (dsts[1] > ALG_DIST_159_2) && (dsts[2] > ALG_DIST_159_2)) ||
          ((ss.p_alg[k] == 9) && (angls[1] > ALG_ANGL_9) && (dsts[1] < ALG_DIST_159_2) && (dsts[2] < ALG_DIST_159_2)) ||
          ((ss.p_alg[k] == 5) && ((angls[1] > ALG_ANGL_5_1) || (angls[1] < ALG_ANGL_5_2)))
          )
          # Stop movement
          ss.w_speeds[:,k] = [0.0; 0.0];
        end
      else
        # Edge robots
        if(ss.p_alg[k] == 2)
          println("  angls[1] = ", angls[1] * 180.0 / pi);
          println("  dsts[1]  = ", dsts[1]);
          println("  dsts[2]  = ", dsts[2]);
        end
        # Assume left-hand side robots by default
        turn = [KILOBOT_WMAX; -KILOBOT_WMAX];
        if ss.p_alg[k] in [2,3,6]
          # Right-hand side robots
          turn = -turn;
        end
        # Deformation test
        if ((angls[1] < ALG_ANGL_234678_1) || (dsts[1] > ALG_DIST_234678_1_1))
          # Turn towards center
          ss.w_speeds[:,k] = turn;
          if(ss.p_alg[k] == 2)
            println("  TURN TOWARDS CENTER: ", turn);
          end
        elseif ((angls[1] > ALG_ANGL_234678_2) || (dsts[1] < ALG_DIST_234678_1_2))
          # Turn away from center
          ss.w_speeds[:,k] = -turn;
          if(ss.p_alg[k] == 2)
            println("  TURN AWAY FROM CENTER: ", -turn);
          end
        elseif((dsts[2] < ALG_DIST_234678_2) || (dsts[3] < ALG_DIST_234678_2))
          # Stop movement
          ss.w_speeds[:,k] = [0.0; 0.0];
          if(ss.p_alg[k] == 2)
            println("  STOP");
          end
        end
      end
    end
  end
end;

########################################

"""
Calculates the forces between robots in the SoftKilobot.
Returns a Matrix{Float64}.
"""
function CalcForces(
  u::Vector{Float64}, # The motion state of the kilobots
  ss::State)          # The SoftKilobot
  fs = zeros(2,ss.n)
  for k in 1:ss.n
    # Robot index in state vector to identify block
    idx1 = (k-1) * FIELDS + 1
    # Go through neighbors
    for n in ss.nbr_springs[k]
      # Avoid useless recalculations
      if k < n
        # Robot index in state vector to identify block
        idx2 = (n-1) * FIELDS + 1
        # Get positions
        p1 = u[POS_X + idx1:POS_Y + idx1,:]
        p2 = u[POS_X + idx2:POS_Y + idx2,:]
        # Calculate relative vector
        p12 = p2 - p1
        # Calculate distance
        dist = norm(p12)
        # Calculate force contribution and add it
        f = ss.S_STIFF * (p12 / dist) * (dist - ss.S_REST)
        fs[:,k] += f
        fs[:,n] -= f
      end
    end
  end
  return fs
end;

"""
The differential steering model.
Returns a tuple (Float64,Float64,Float64).
"""
function DiffSteering(
  lspeed::Float64, # Left wheel speed [m/s]
  rspeed::Float64, # Right wheel speed [m/s]
  ornt::Float64,   # Current robot orientation [rad]
  axle::Float64)   # Axle length [m]
  vx = (lspeed + rspeed) * cos(ornt) / 2.0
  vy = (lspeed + rspeed) * sin(ornt) / 2.0
  w  = (rspeed - lspeed) / axle
  return vx, vy, w
end;

"""
Updates the state of a SoftKilobot.
"""
function UpdateState(
  du::Vector{Float64},     # Difference in motion state of the kilobots
  u::Vector{Float64},      # Motion state of the kilobots
  ss::State,               # The SoftKilobot
  fs::Matrix{Float64})     # Forces [N]
  for k in 1:ss.n
    # Robot index in state vector to identify block
    idx = (k-1) * FIELDS + 1
    # Linear/rotational velocities set in previous step
    vx_old = u[VEL_X + idx]
    vy_old = u[VEL_Y + idx]
    vw_old = u[VEL_W + idx]
    # Calculate linear speeds
    vx, vy, vw = DiffSteering(
      ss.w_speeds[LWHEEL,k] + ss.w_bias[1,k],
      ss.w_speeds[RWHEEL,k] + ss.w_bias[2,k],
      u[ORNT + idx],
      ss.AXLE)
    # Calculate velocity update
    ss.vs[1,k] = -ss.V_DAMP * vx_old + (vx - vx_old + fs[1,k]) / ss.MASS
    ss.vs[2,k] = -ss.V_DAMP * vy_old + (vy - vy_old + fs[2,k]) / ss.MASS
    ss.vs[1,k] = clamp(ss.vs[1,k], -ss.F_MAX, ss.F_MAX)
    ss.vs[2,k] = clamp(ss.vs[2,k], -ss.F_MAX, ss.F_MAX)
    # State update
    du[VEL_X + idx] = ss.vs[1,k] - vx_old
    du[VEL_Y + idx] = ss.vs[2,k] - vy_old
    du[VEL_W + idx] = (vw - vw_old) / ss.I
    du[POS_X + idx] = vx_old
    du[POS_Y + idx] = vy_old
    du[ORNT  + idx] = vw_old
  end
  return du
end;

########################################

"""
Step function to calculate the new state of the SoftKilobot.
"""
function Step!(
  du, # Difference in motion state of the kilobots
  u,  # Motion state of the kilobots
  ss, # The SoftKilobot
  t)  # Time [s]
  # Calculate target speeds using algorithm
  Alg(u, ss, t)
  # Calculate forces due to springs
  fs = CalcForces(u, ss);
  # Update states
  du = UpdateState(du, u, ss, fs)
  # All done
  return du
end;

########################################

"""
Calculates the center of mass of a SoftKilobot
"""
function CoM(
  u::Vector{Float64}, # Motion state of the kilobots
  ss::State)          # The SoftKilobot
  com = zeros(2)
  for k in 1:ss.n
    # Robot index in state vector to identify block
    idx = (k-1) * FIELDS + 1
    # Calculate CoM
    com[1] += u[POS_X + idx]
    com[2] += u[POS_Y + idx]
  end
  return com / ss.n
end;

########################################

"""
Plots the SoftKilobot
"""
function Plot(
  x::Matrix{Float64}, # The poses of the SoftKilobot robots
  ss::State,          # The SoftKilobot
  t::Float64;         # The current time
  vel=false,          # Whether to plot the velocity vectors
  nbr_springs=false,  # Whether to plot the springs
  nbr_alg_ds=false,   # Whether to plot the algorithm-related distance links
  nbr_alg_as=false)   # Whether to plot the algorithm-related angle links
  # General settings
  p = plot(title=string("t = ", t), legend=:none)
  plot!(p, xlims=(0,2), ylims=(0,2))
  # Start and end lines
  plot!(p, [ENV_START ENV_END; ENV_START ENV_END], [0 0; ENV_SIDE ENV_SIDE])
  # Robots as circles
  scatter!(p, x[POS_X+1,:], x[POS_Y+1,:])
  # Plot velocities
  if(vel)
    i = 2
    plot!(p, [x[POS_X+1,i]; (x[POS_X+1,i]+x[VEL_X+1,i]*100)], [x[POS_Y+1,i]; (x[POS_Y+1,i]+x[VEL_Y+1,i]*100)], label="new")
    # for i in 1:ss.n
    #   # New velocity
    #   plot!(p, [x[POS_X+1,i]; (x[POS_X+1,i]+x[VEL_X+1,i]*100)], [x[POS_Y+1,i]; (x[POS_Y+1,i]+x[VEL_Y+1,i]*100)], label="new")
    #   # Old velocity
    #   plot!(p, [x[POS_X+1,i]; (x[POS_X+1,i]+ss.vs[1,i])], [x[POS_Y+1,i]; (x[POS_Y+1,i]+ss.vs[2,i])], label="old")
    # end
  end
  # Plot springs as segments
  if(nbr_springs)
    for i in 1:ss.n
      for j in ss.nbr_springs[i]
        if(i < j)
          plot!(p, [x[POS_X+1,i]; x[POS_X+1,j]], [x[POS_Y+1,i]; x[POS_Y+1,j]])
        end
      end
    end
  end
  # Plot algorithm neighbors as segments
  if(nbr_alg_ds)
    for i in 1:ss.n
      for j in ss.nbr_alg_ds[i]
        plot!(p, [x[POS_X+1,i]; x[POS_X+1,j]], [x[POS_Y+1,i]; x[POS_Y+1,j]])
      end
    end
  end
  # Plot algorithm neighbors as segments
  if(nbr_alg_as)
    for i in 1:ss.n
      for j in ss.nbr_alg_as[i]
        plot!(p, [x[POS_X+1,i]; x[POS_X+1,j]], [x[POS_Y+1,i]; x[POS_Y+1,j]])
      end
    end
  end
  return p
end;

########################################

function PlotTrajectories(
  x::Matrix{Float64}, # The poses of the SoftKilobot robots
  ss::State)          # The SoftKilobot
    # General settings
    p = plot(title=string("t = ", t))
    plot!(p, xlims=(0,2), ylims=(0,2), legend=:none)
    # Start and end lines
    plot!(p, [ENV_START ENV_END; ENV_START ENV_END], [0 0; ENV_SIDE ENV_SIDE])  
end;

########################################

end # module SoftKilobot
