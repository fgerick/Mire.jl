using Documenter, Mire

makedocs(sitename = "Mire.jl - Modes In Rotating Ellipsoids Julia package";

        pages = ["Home" => "index.md",
                "Examples" => ["Inertial modes in Ellipsoid" => "man/examples/inertialmodes.md", "Modified Malkus modes in Ellipsoid" => "man/examples/modifiedmalkusmodes.md"], "Functions" => "man/functions.md"])
deploydocs(
    repo = "github.com/fgerick/Mire.jl.git"
)
