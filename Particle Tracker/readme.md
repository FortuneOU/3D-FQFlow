
# Particle Tracking Simulation in Velocity Fields

This repository provides a Python script for 3D particle tracking simulation based on velocity or acoustic field data. It reads sequential `.vtu` mesh files, calculates particle trajectories using 2nd-order Runge-Kutta (RK2) integration, handles stagnant particles by reseeding them at the inlets, and supports resuming from saved breakpoints.

## 🛠 Dependencies

Before running the script, please ensure you have the required Python libraries installed:

```bash
pip install numpy scipy pyvista
```

## 🚀 How to Use

To adapt this script to your own data, you need to configure the **folder path** for your `.vtu` files and determine the **inlet coordinates** of your vessel model using ParaView.

### Step 1: Modify the Data Folder Path
By default, the script reads velocity field files named `result_*.vtu` from a specific directory. 

Open the Python script, locate the `main_optimized` function, and change the `folder` variable to the actual path where your `.vtu` files are stored:

```python
def main_optimized(resume_file=None):  
    # CHANGE THIS to the folder path containing your .vtu files
    folder = "path/to/your/VelocityField"  
    file_pattern = os.path.join(folder, "result_*.vtu") 
```

### Step 2: Find Vessel Inlet Coordinates using ParaView
To ensure particles are correctly spawned at the vessel inlets, you need to find their exact spatial coordinates:
1. Open one of your `.vtu` files using **ParaView**.
2. Make sure the 3D vessel structure is visible.
3. Use the **Hover Cells On** tool (hover your mouse over the inlet) or use the **Probe Location** filter to inspect the coordinates of the inlet.
4. Note down the coordinate ranges (e.g., X ranges from -5.0 to -4.5, Y ranges from -5.2 to -4.7, and the Z plane is at -4.8). 
*(Note: If your model has multiple inlets, please record the coordinate ranges for each inlet.)*

### Step 3: Update Inlet Positions in the Code
Once you have obtained the coordinates from ParaView, you must update the random generation ranges in the Python script. There are **TWO** places that require modification:

**1. Initial Particle Seeding**
Find the `Initial Particles` section inside the `main_optimized` function and replace the `np.random.uniform` values with your measured ranges:
```python
# Configuration for Inlet 1
particles1 = np.column_stack([  
    np.random.uniform(X_MIN_1, X_MAX_1, n_particles1), # Replace with Inlet 1 X range
    np.random.uniform(Y_MIN_1, Y_MAX_1, n_particles1), # Replace with Inlet 1 Y range
    np.full(n_particles1, Z_VAL_1)                     # Replace with Inlet 1 Z coordinate
])  

# Configuration for Inlet 2 (If you have multiple inlets)
particles2 = np.column_stack([  
    np.random.uniform(X_MIN_2, X_MAX_2, n_particles2), # Replace with Inlet 2 X range
    np.random.uniform(Y_MIN_2, Y_MAX_2, n_particles2), # Replace with Inlet 2 Y range
    np.full(n_particles2, Z_VAL_2)                     # Replace with Inlet 2 Z coordinate
])  
```

**2. Stagnant Particle Reseeding**
When particles stagnate (speed drops near zero), the script will respawn them at the inlets. Scroll down to the **Stagnant particles judgment and reset** section inside the main loop and update the coordinates to match exactly what you set above:
```python
# Reseed at Inlet 1
current_positions[idx1, 0] = np.random.uniform(X_MIN_1, X_MAX_1, len(idx1))  
current_positions[idx1, 1] = np.random.uniform(Y_MIN_1, Y_MAX_1, len(idx1))  
current_positions[idx1, 2] = Z_VAL_1  

# Reseed at Inlet 2
if (~mask1).any(): 
    idx2 = reseed_indices[~mask1]  
    current_positions[idx2, 0] = np.random.uniform(X_MIN_2, X_MAX_2, len(idx2))  
    current_positions[idx2, 1] = np.random.uniform(Y_MIN_2, Y_MAX_2, len(idx2))  
    current_positions[idx2, 2] = Z_VAL_2  
```

## 💻 Execution & Breakpoints

**Start a New Simulation:**
Run the script directly via terminal or command prompt:
```bash
python script_name.py
```

**Resume from a Breakpoint:**
The script automatically saves particle trajectories into `.mat` files every 5 steps. If the execution is interrupted or crashes, you do not need to start over. You can resume the simulation from the last saved `.mat` file by passing `resume` and the filename as arguments:
```bash
python script_name.py resume particle_trajectories_cycle01_part02.mat
```
