using Oceananigans
using Oceananigans.Units
using Printf
using AquaSkyLES

Nx = Nz = 128
Lz = 4 * 1024
grid = RectilinearGrid(size=(Nx, Nz), x=(0, 2Lz), z=(0, Lz), topology=(Periodic, Flat, Bounded))

p₀ = 101325 # Pa
θ₀ = 288 # K
reference_state = AquaSkyLES.ReferenceState(base_pressure=p₀, potential_temperature=θ₀)
buoyancy = AquaSkyLES.MoistAirBuoyancy(; reference_state)

ρ₀ = AquaSkyLES.base_density(buoyancy) # air density at z=0
cₚ = buoyancy.thermodynamics.dry_air.heat_capacity
Q₀ = 1000 # heat flux in W / m²
Jθ = Q₀ / (ρ₀ * cₚ) # temperature flux
θ_bcs = FieldBoundaryConditions(bottom=FluxBoundaryCondition(Jθ))

vapor_flux = FluxBoundaryCondition(1e-3)
q_bcs = FieldBoundaryConditions(bottom=vapor_flux)

advection = WENO() #(momentum=WENO(), θ=WENO(), q=WENO(bounds=(0, 1)))
tracers = (:θ, :q)
model = NonhydrostaticModel(; grid, advection, buoyancy,
                            tracers = (:θ, :q),
                            boundary_conditions = (θ=θ_bcs, q=q_bcs))

Lz = grid.Lz
Δθ = 5 # K
Tₛ = reference_state.θ # K
θᵢ(x, z) = Tₛ + Δθ * z / Lz + 1e-2 * Δθ * randn()
qᵢ(x, z) = 1e-2 + 1e-5 * rand()
set!(model, θ=θᵢ, q=qᵢ)

simulation = Simulation(model, Δt=10, stop_time=2hours)
conjure_time_step_wizard!(simulation, cfl=0.7)

T = AquaSkyLES.TemperatureField(model)
qˡ = AquaSkyLES.CondensateField(model, T)
qᵛ★ = AquaSkyLES.SaturationField(model, T)
δ = Field(model.tracers.q - qᵛ★)

function progress(sim)
    compute!(T)
    compute!(qˡ)
    compute!(δ)
    q = sim.model.tracers.q
    θ = sim.model.tracers.θ
    u, v, w = sim.model.velocities

    umax = maximum(u)
    vmax = maximum(v)
    wmax = maximum(w)

    qmin = minimum(q)
    qmax = maximum(q)
    qˡmax = maximum(qˡ)
    δmax = maximum(δ)

    θmin = minimum(θ)
    θmax = maximum(θ)

    msg = @sprintf("Iter: %d, t = %s, max|u|: (%.2e, %.2e, %.2e)",
                    iteration(sim), prettytime(sim), umax, vmax, wmax)

    msg *= @sprintf(", extrema(q): (%.2e, %.2e), max(qˡ): %.2e, min(δ): %.2e, extrema(θ): (%.2e, %.2e)",
                     qmin, qmax, qˡmax, δmax, θmin, θmax)

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

outputs = merge(model.velocities, model.tracers, (; T, qˡ, qᵛ★))

ow = JLD2Writer(model, outputs,
                filename = "free_convection.jld2",
                schedule = IterationInterval(10),
                overwrite_existing = true)

simulation.output_writers[:jld2] = ow

run!(simulation)

θt = FieldTimeSeries("free_convection.jld2", "θ")
Tt = FieldTimeSeries("free_convection.jld2", "T")
qt = FieldTimeSeries("free_convection.jld2", "q")
qˡt = FieldTimeSeries("free_convection.jld2", "qˡ")
times = qt.times
Nt = length(θt)

using GLMakie, Printf

n = Observable(length(θt))

θn = @lift θt[$n]
qn = @lift qt[$n]
Tn = @lift Tt[$n]
qˡn = @lift qˡt[$n]
title = @lift "t = $(prettytime(times[$n]))"

fig = Figure(size=(800, 400), fontsize=12)
axθ = Axis(fig[1, 1], xlabel="x (m)", ylabel="z (m)")
axq = Axis(fig[1, 2], xlabel="x (m)", ylabel="z (m)")
axT = Axis(fig[2, 1], xlabel="x (m)", ylabel="z (m)")
axqˡ = Axis(fig[2, 2], xlabel="x (m)", ylabel="z (m)")

fig[0, :] = Label(fig, title, fontsize=22, tellwidth=false)

Tmin = minimum(Tt)
Tmax = maximum(Tt)

hmθ = heatmap!(axθ, θn, colorrange=(Tₛ, Tₛ+Δθ))
hmq = heatmap!(axq, qn, colorrange=(0, 2e-2), colormap=:magma)
hmT = heatmap!(axT, Tn, colorrange=(Tmin, Tmax))
hmqˡ = heatmap!(axqˡ, qˡn, colorrange=(0, 2e-4), colormap=:magma)

# Label(fig[0, 1], "θ", tellwidth=false)
# Label(fig[0, 2], "q", tellwidth=false)
# Label(fig[0, 1], "θ", tellwidth=false)
# Label(fig[0, 2], "q", tellwidth=false)

Colorbar(fig[1, 0], hmθ, label = "θ [K]", vertical=true)
Colorbar(fig[1, 3], hmq, label = "q", vertical=true)
Colorbar(fig[2, 0], hmT, label = "T [K]", vertical=true)
Colorbar(fig[2, 3], hmqˡ, label = "qˡ", vertical=true)

fig

record(fig, "free_convection.mp4", 1:Nt, framerate=12) do nn
    n[] = nn
end
