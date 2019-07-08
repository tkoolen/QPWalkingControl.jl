module QPWalkingControl

export
    HumanoidQPController,
    ICPController,
    ICPTrajectoryGenerator,
    ICPWalkingStateMachine,
    PDCoMController,
    PushApplier,
    Widget, # from InteractBase
    PushRecoveryVisualizer,
    add_sole_frames!,
    make_foot_polygons,
    icp,
    horizontal_projection,
    transfer_weight!

# consider moving to RBD.jl
export
    in_frame,
    @framechecked,
    unwrap

# consider moving to RigidBodySim.jl
export
    JointDamping

using LinearAlgebra
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using StaticArrays
using Rotations
using QPControl
using QPControl.Trajectories
using Parametron
using PlanarConvexHulls
using StaticUnivariatePolynomials
import MathOptInterface

using Base.Iterators: cycle
using RigidBodyDynamics.Graphs: TreePath, target
using RigidBodyDynamics: frame_definition
using StaticUnivariatePolynomials: derivative
using QPControl.Trajectories: exponential_integral, breaks, subfunctions

const SUP = StaticUnivariatePolynomials
const MOI = MathOptInterface
const RBD = RigidBodyDynamics

const SPoint3D{T} = Point3D{SVector{3, T}}

include("frames.jl")
include("util.jl")

include("setup.jl")
include("planarcontactmode.jl")
include("controller.jl")

include("damping.jl")
include("pushapplier.jl")

include("foot_trajectories/basic_foot_trajectory.jl")
include("foot_trajectories/foot_trajectory_gen.jl")

include("icp/icp.jl")
include("icp/control.jl")
include("icp/trajectory_gen.jl")
include("icp/trajectory_init.jl")
include("icp/statemachine.jl")

include("com_tracking/control.jl")

include("visualization.jl")

using .Visualization

end # module
