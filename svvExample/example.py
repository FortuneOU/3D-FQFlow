import pyvista as pv
from svv.domain.domain import Domain
from svv.tree.tree import Tree
 from svv.simulation.simulation import Simulation


# Creating the Tissue Domain
cube = Domain(pv.Cube(center=(0.0, 0.0, 0.0), x_length=50.0, y_length=50.0, z_length=50.0,))
cube.create()
cube.solve()
cube.build()

# Creating the Vascular Tree Object
t = Tree()
t.set_domain(cube)
t.set_root()
t.n_add(5)

# Create Simulation Container for Synthetic Tree Object
 sim = Simulation(t)

 # Build All Surface/Volume Meshes for 3D CFD Simulation
 sim.build_meshes()

# Detect walls, outlets, and shared boundaries
sim.extract_faces(crease_angle=60.0)

# Inspect names stored in the GeneralMesh container
print(sim.fluid_domain_meshes[0].faces.keys())

# Build and persist 3-D CFD input files
sim.construct_3d_fluid_simulation()
sim.write_3d_fluid_simulation()

# Result: meshes + XML under sim.file_path (e.g. ./simulations/example/)

from svv.simulation.fluid.rom.zero_d.zerod_tree import export_0d_simulation

export_0d_simulation(t, outdir=sim.file_path + "/rom_0d")

