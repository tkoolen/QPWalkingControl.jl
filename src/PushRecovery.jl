module PushRecovery

export
    PushRecoveryController,
    ICPController,
    ICPTrajectoryGenerator,
    PushApplier,
    Widget, # from InteractBase
    add_sole_frames!,
    foot_polygons

using LinearAlgebra
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using RigidBodyDynamics.PDControl
using StaticArrays
using Rotations
using QPControl
using Parametron
using PlanarConvexHulls
import MathOptInterface

using RigidBodyDynamics.Graphs: TreePath, target
using RigidBodyDynamics: frame_definition

const MOI = MathOptInterface
const RBD = RigidBodyDynamics

const Vec2{T} = SVector{2, T}
const SPoint3D{T} = Point3D{SVector{3, T}}

include("util.jl")

include("setup.jl")

include("icpcontrol.jl")
include("icptrajectory.jl")
include("contactmode.jl")
include("controller.jl")

include("pushapplier.jl")
include("visualization.jl")

using .Visualization

end # module
