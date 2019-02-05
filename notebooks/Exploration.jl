#!/usr/bin/env python
# coding: utf-8

# In[1]:


using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


# In[ ]:


using LinearAlgebra
using QPControl
using RigidBodyDynamics
using RigidBodyDynamics.PDControl
using RigidBodyDynamics.Contact
using StaticArrays
using AtlasRobot
using Test
using RigidBodySim
using MathOptInterface
using OSQP.MathOptInterfaceOSQP: OSQPSettings
const MOI = MathOptInterface
using OSQP
using PushRecovery
BLAS.set_num_threads(1)


# In[3]:


@time AtlasRobot.mechanism(add_flat_ground=true);


# In[4]:


@time mechanism = AtlasRobot.mechanism(add_flat_ground=true);


# In[5]:


# create low level controller
lowlevel = let
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-5)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-5)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 5000)
    MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    lowlevel = MomentumBasedController{4}(mechanism, optimizer,
        floatingjoint = findjoint(mechanism, "pelvis_to_world"));
    for body in bodies(mechanism)
        for point in RigidBodyDynamics.contact_points(body)
            position = location(point)
            normal = FreeVector3D(default_frame(body), 0.0, 0.0, 1.0)
            μ = point.model.friction.μ
            contact = addcontact!(lowlevel, body, position, normal, μ)
            contact.maxnormalforce[] = 1e6 # TODO
            contact.weight[] = 1e-3
        end
    end
    lowlevel
end;


# In[6]:


nominalstate = MechanismState(mechanism)
AtlasRobot.setnominal!(nominalstate)


# In[7]:


# ICP stuff
icptraj = let
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, OSQPSettings.Verbose(), false)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-6)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-8)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 10000)
    MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    foot_polygon_sides = 4
    num_segments = 1
    ICPTrajectoryGenerator{Float64, foot_polygon_sides}(optimizer, num_segments)
end
linear_momentum_controller = ICPController(mechanism, icptraj, center_of_mass(nominalstate).v[3] - 0.05);


# In[8]:


# create high level controller
feet = findbody.(Ref(mechanism), ["l_foot", "r_foot"])
pelvis = findbody(mechanism, "pelvis")
controller = PushRecoveryController(lowlevel, feet, pelvis, nominalstate, linear_momentum_controller);


# In[9]:


# create visualizer
using MeshCat
using MeshCatMechanisms

if !(@isdefined gui) || !(@isdefined pushapplier) || !any(isopen, gui.visualizer.visualizer.core.scope.pool.connections)
    pushapplier = PushApplier(mechanism, Point3D(default_frame(pelvis), 0.0, 0.0, 0.0));    
    visuals = URDFVisuals(AtlasRobot.urdfpath(); package_path = [AtlasRobot.packagepath()])
    mvis = MechanismVisualizer(mechanism, visuals)
    gui = GUI(mvis, usernode=Widget(pushapplier, max_force=100.0, max_Δt=0.3))
    open(gui)
end

set_configuration!(gui.visualizer, configuration(nominalstate))


# In[10]:


# create ODEProblem
state = MechanismState(mechanism)
AtlasRobot.setnominal!(state)
Δt = 1 / 300
pcontroller = PeriodicController(similar(velocity(state)), Δt, controller)
# TODO: add damping
dynamics = Dynamics(mechanism, SumController(similar(velocity(state)), (pcontroller, pushapplier)))
callback = CallbackSet(RealtimeRateLimiter(poll_interval=pi / 100), CallbackSet(gui))
problem = ODEProblem(dynamics, state, (0., Inf); callback=callback)


# In[ ]:


# simulate
@time sol = solve(problem, Tsit5(), abs_tol = 1e-8, dt = 1e-6);
last(sol.t)


# In[ ]:


setanimation!(mvis, sol)


# In[ ]:


using BenchmarkTools
AtlasRobot.setnominal!(state)
τ = similar(velocity(state));
benchresult = @benchmark $controller($τ, 0.0, $state)
@test benchresult.allocs == 0
benchresult

