# Authoring script: pack examples/data/*.nc into a GitHub-Release-hosted
# artifact and bind it into examples/Artifacts.toml.
#
# Usage: julia --project=<repo root> examples/scripts/pack_data.jl [tag]
#
# After it succeeds, upload the tarball (CONFIRM with the user first):
#   gh release create <tag> examples/example-data.tar.gz --title "<tag>" \
#       --notes "Example data for the UXarray comparison tutorial"
import Pkg
import SHA  # Pkg.SHA is not accessible as a module property on this Julia (1.10)

const TAG = length(ARGS) >= 1 ? ARGS[1] : "example-data-v1"
const URL = "https://github.com/VANvonZHANG/ManifoldFields.jl/releases/download/$(TAG)/example-data.tar.gz"

const EXAMPLES = dirname(@__DIR__)
const DATA = joinpath(EXAMPLES, "data")
const TOML = joinpath(EXAMPLES, "Artifacts.toml")
const TARBALL = joinpath(EXAMPLES, "example-data.tar.gz")

const FILES = [
    "oQU480.grid.nc",
    "oQU480.data.nc",
    "oQU480.ugrid.nc",
    "psi_healpix.nc",
    "latlon_analytic.nc",
    "cubedsphere_analytic.nc",
    "healpix_analytic.nc",
]

for f in FILES
    isfile(joinpath(DATA, f)) ||
        error("missing $(joinpath(DATA, f)) — run prepare_ugrid.py and make_julia_data.jl first")
end

hash = Pkg.create_artifact() do dir
    for f in FILES
        cp(joinpath(DATA, f), joinpath(dir, f))
    end
end

adir = Pkg.artifact_path(hash)
# Flat layout: Julia 1.10's artifact unpacker does NOT strip a top-level
# directory, so files must sit at the archive root for the unpacked tree
# hash to match the manifest's git-tree-sha1.
run(`tar -czf $TARBALL -C $adir $FILES`)
sha = bytes2hex(open(SHA.sha256, TARBALL))
Pkg.Artifacts.bind_artifact!(TOML, "example-data", hash;
                             download_info = [(URL, sha)], lazy = true, force = true)
println("artifact tree hash: ", hash)
println("tarball sha256:     ", sha)
println("tarball:            ", TARBALL, " (", round(filesize(TARBALL) / 1e6; digits = 1), " MB)")
println("manifest:           ", TOML)
println("upload:  gh release create $(TAG) $(TARBALL) --title \"$(TAG)\" --notes \"Example data for the UXarray comparison tutorial\"")
