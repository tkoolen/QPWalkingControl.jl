module PushRecovery

export
    PushRecoveryController,
    ICPController,
    ICPTrajectoryGenerator,
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
using QPControl.Trajectories: exponential_integral

const SUP = StaticUnivariatePolynomials
const MOI = MathOptInterface
const RBD = RigidBodyDynamics

const Vec2{T} = SVector{2, T}
const SPoint3D{T} = Point3D{SVector{3, T}}

include("frames.jl")
include("util.jl")

include("setup.jl")
include("planarcontactmode.jl")
include("controller.jl")

include("pushapplier.jl")

include("icp/icp.jl")
include("icp/control.jl")
include("icp/trajectory_gen.jl")
include("icp/trajectory_init.jl")

include("foot_trajectories/basic_foot_trajectory.jl")

include("visualization.jl")

using .Visualization

end # module
