module PushRecovery

export
    PushRecoveryController,
    ICPController,
    ICPTrajectoryGenerator,
    PushApplier,
    Widget, # from InteractBase
    add_sole_frames!,
    foot_polygons,
    init_back_and_forth_shift!

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
import MathOptInterface

using Base.Iterators: cycle
using RigidBodyDynamics.Graphs: TreePath, target
using RigidBodyDynamics: frame_definition

const MOI = MathOptInterface
const RBD = RigidBodyDynamics

const Vec2{T} = SVector{2, T}
const SPoint3D{T} = Point3D{SVector{3, T}}

include("util.jl")

include("setup.jl")
include("planarcontactmode.jl")
include("controller.jl")

include("pushapplier.jl")
include("visualization.jl")

include("icp/icp.jl")
include("icp/control.jl")
include("icp/trajectory_gen.jl")
include("icp/trajectory_init.jl")

using .Visualization

end # module
