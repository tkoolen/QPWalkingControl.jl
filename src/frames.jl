struct FrameAnnotated{T}
    frame::CartesianFrame3D
    val::T
end

function Base.show(io::IO, x::FrameAnnotated)
    print(io, x.val)
    print(io, " (in frame: $(x.frame))")
end

in_frame(frame::CartesianFrame3D, val) = FrameAnnotated(frame, val)

unwrap(x) = x
unwrap(x::FrameAnnotated) = x.val

@inline function same_frame(args::Tuple, frame::Union{Missing, CartesianFrame3D}=missing)
    isempty(args) && return frame
    head = args[1]
    tail = Base.tail(args)
    if head isa FrameAnnotated
        if frame !== missing && head.frame !== frame
            throw(ArgumentError("Frame mismatch.")) # TODO: nicer message
        end
        return same_frame(tail, head.frame)
    else
        return same_frame(tail, frame)
    end
end

function framechecked(f, args...)
    frame = same_frame(args)
    in_frame(frame, f(map(unwrap, args)...))
end

macro framechecked(expr)
    expr.head == :call || return :(error("Expected a function call"))
    f = esc(expr.args[1])
    args = map(esc, expr.args[2 : end])
    return quote
        framechecked($f, $(args...))
    end
end
