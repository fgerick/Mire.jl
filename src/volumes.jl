"""
    Volume{T}

Defines `abstract type` of the volume.
"""
abstract type Volume{T} end


"""
    Ellipsoid{T<:Number} <: Volume{T}

`Volume` type `Ellipsoid`. Create by calling `Ellipsoid(a,b,c)`, with `a,b,c` the semi-axes.
`a,b,c` can be of any `Number` type.

Examples:
`Ellipsoid(1.1,1.0,0.9)` is a `Ellipsoid{Float64}`,
`Ellipsoid(1//1,1//2,1//5)` is a `Ellispoid{Rational{Int64}}`.
"""
struct Ellipsoid{T} <: Volume{T}
    a::T
    b::T
    c::T
end

Ellipsoid(a, b, c) = Ellipsoid(promote(a, b, c)...)
Ellipsoid(a::Int, b::Int, c::Int) = Ellipsoid(a // 1, b // 1, c // 1)


"""
    Sphere{T<:Number} <: Volume{T}

`Volume` type `Sphere`. Create by calling `Sphere{T}()`, with `T` any `Number` type.
Default: `Sphere(T)` gives a `Sphere{Float64}()`. For other types use, e.g.
`Sphere{Rational{BigInt}}()`.
"""
struct Sphere{T} <: Volume{T} end

Sphere() = Sphere{Float64}()
