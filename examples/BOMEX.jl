using Oceananigans
using Oceananigans.Units
using Printf
using AquaSkyLES

Nx = 64
Nz = 75
Ny = 64
Lz = Nz * 40
Lz = Nz * 100

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, 2Lz),
                       y = (0, 2Lz),
                       z = (0, Lz),
                       topology = (Periodic, Periodic, Bounded))

p₀ = 101325 # Pa
θ₀ = 285 # K
reference_state = AquaSkyLES.ReferenceState(base_pressure=p₀, potential_temperature=θ₀)
buoyancy = AquaSkyLES.MoistAirBuoyancy(; reference_state)
thermodynamics = buoyancy.thermodynamics
condensation = thermodynamics.condensation
