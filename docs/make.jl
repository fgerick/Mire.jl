using Documenter, Mire

makedocs(sitename = "Mire.jl - Modes In Rotating Ellipsoids Julia package")
deploydocs(
    repo = "github.com/fgerick/Mire.jl.git",
    push_preview=false
)
