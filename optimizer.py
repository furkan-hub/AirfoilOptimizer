import casadi as ca
import aerosandbox as asb
import matplotlib.pyplot as plt
import numpy as np

# Constants and Parameters
W_FUSELAGE = 68.64  # Fuselage weight [N]
K_WING = 40  # Wing weight coefficient [N/(m²*AR)]
RHO = 1.225  # Air density [kg/m³]
G = 9.81  # Gravity [m/s²]
E_OSWALD = 0.9  # Oswald efficiency factor
REYNOLDS = 1e6  # Reynolds number
MACH = 0.1  # Mach number
SIGMA_ALLOWABLE = 100e6  # Allowable bending stress [Pa]

# NACA 4-digit airfoils to analyze
#NACA_AIRFOILS = ["0012", "2412", "4412", "23012"] 
NACA_AIRFOILS= [str(i).zfill(4) for i in range(0,9999,4)]

def get_airfoil_params(naca_code):
    """Analyze airfoil and extract CL_max and CD0 using AeroSandbox."""
    af = asb.Airfoil(name="naca"+naca_code)
    alphas = np.linspace(-5, 15, 50)
    cl_list = []
    cd_list = []
    
    for alpha in alphas:
        try:
            analysis = af.get_aero_from_neuralfoil(
                alpha=alpha, Re=REYNOLDS, mach=MACH, model_size="xlarge"
            )
            cl_list.append(analysis["CL"])
            cd_list.append(analysis["CD"])
        except:
            cl_list.append(0)
            cd_list.append(0)
    
    cl_array = np.array(cl_list)
    cd_array = np.array(cd_list)
    valid = np.isfinite(cl_array) & np.isfinite(cd_array)
    
    max_cl_idx = np.nanargmax(cl_array[valid]) if valid.any() else 0
    cl_max = cl_array[valid][max_cl_idx] if valid.any() else 0
    cd0 = np.min(cd_array[valid]) if valid.any() else 0.05
    
    return {"CL_max": cl_max, "CD0": cd0}

def optimize_wing(airfoil_params):
    """Run optimization with structural constraints."""
    opti = ca.Opti()
    
    # Design variables
    S = opti.variable()  # Wing area [m²]
    AR = opti.variable() # Aspect ratio
    V = opti.variable()  # Cruise speed [m/s]
    
    # Derived geometry
    b = ca.sqrt(S * AR)  # Wingspan [m]
    c = S / b            # Mean chord [m]
    
    # Total weight calculation (including wing weight)
    W_wing = K_WING * S * AR
    W_total = W_FUSELAGE + W_wing
    
    # Aerodynamic model
    CL = (2 * W_total) / (RHO * S * V**2)
    CD = airfoil_params["CD0"] + (CL**2) / (ca.pi * AR * E_OSWALD)
    L_over_D = CL / CD
    
    # Structural bending constraint
    bending_moment = W_total * b / 8  # Simplified bending moment
    sigma = (bending_moment * 0.12*c) / ((1/12)*c*(0.12*c)**3)  # Stress calculation
    opti.subject_to(sigma <= SIGMA_ALLOWABLE)
    
    # Objective and constraints
    opti.minimize(-L_over_D)
    opti.subject_to([
        V >= ca.sqrt((2 * W_total) / (RHO * S * airfoil_params["CL_max"])),  # Stall speed
        S >= 0.300, S <= 0.320,                   # Practical wing area limits
        AR >= 4, AR <= 20,                   # Realistic aspect ratio bounds
        V >= 20, V <= 30,                   # Speed limits
        CL <= 0.95*airfoil_params["CL_max"]   # CL margin
    ])
    
    # Initial guess
    opti.set_initial(S, 30)
    opti.set_initial(AR, 8)
    opti.set_initial(V, 100)
    
    # Solve
    opti.solver("ipopt")
    try:
        sol = opti.solve()
        return {
            "S": sol.value(S),
            "AR": sol.value(AR),
            "V": sol.value(V),
            "L_over_D": sol.value(L_over_D),
            "W_total": sol.value(W_total)
        }
    except:
        return None
import casadi as ca
import aerosandbox as asb
import matplotlib.pyplot as plt
import numpy as np

# ... [Keep all constants and functions identical until the visualization section] ...

# Analysis and visualization
results = {}
for naca in NACA_AIRFOILS:
    print(f"Processing {naca}...")
    params = get_airfoil_params(naca)
    if params["CL_max"] < 1.0:
        print(f"Skipping {naca} - insufficient CL")
        continue
    opt_result = optimize_wing(params)
    if opt_result:
        results[naca] = opt_result

# Sort airfoils by L/D ratio and select top 5
sorted_airfoils = sorted(results.items(), key=lambda x: -x[1]["L_over_D"])[:5]
airfoil_names = [item[0] for item in sorted_airfoils]
airfoil_data = [item[1] for item in sorted_airfoils]

# Create condensed visualization
plt.figure(figsize=(14, 10))

# L/D Ratio
plt.subplot(2, 2, 1)
plt.bar(airfoil_names, [d["L_over_D"] for d in airfoil_data])
plt.title("Top 5 Airfoils: Lift-to-Drag Ratio")
plt.ylabel("L/D")
plt.grid(True, alpha=0.3)

# Aspect Ratios
plt.subplot(2, 2, 2)
plt.bar(airfoil_names, [d["AR"] for d in airfoil_data])
plt.axhline(12, color='r', linestyle='--', label="Typical Transport AR")
plt.title("Aspect Ratio Comparison")
plt.ylabel("AR")
plt.legend()
plt.grid(True, alpha=0.3)

# Wing Loading
plt.subplot(2, 2, 3)
wing_loadings = [d["W_total"]/d["S"]/G for d in airfoil_data]
plt.bar(airfoil_names, wing_loadings)
plt.title("Wing Loading Analysis")
plt.ylabel("Wing Loading (kg/m²)")
plt.grid(True, alpha=0.3)

# Cruise Speed
plt.subplot(2, 2, 4)
plt.bar(airfoil_names, [d["V"] for d in airfoil_data])
plt.title("Optimal Cruise Speeds")
plt.ylabel("Speed (m/s)")
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Detailed visualization of best airfoil
best_airfoil = airfoil_names[0]
best_data = airfoil_data[0]

# Airfoil geometry plot
af = asb.Airfoil(naca4=best_airfoil if len(best_airfoil) == 4 else "2412")
x = np.linspace(0, 1, 200)
yu = [af.upper_coordinates(xi)[1] for xi in x]
yl = [af.lower_coordinates(xi)[1] for xi in x]

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].plot(x, yu, 'b', label='Upper Surface')
ax[0].plot(x, yl, 'r', label='Lower Surface')
ax[0].set_title(f"{best_airfoil} Airfoil Geometry")
ax[0].set_xlabel("Chord")
ax[0].set_ylabel("Thickness")
ax[0].grid(True)
ax[0].axis("equal")
ax[0].legend()

# Performance parameters table
cell_text = [
    [f"{best_data['S']:.1f} m²"],
    [f"{best_data['AR']:.1f}"],
    [f"{best_data['V']:.1f} m/s"],
    [f"{best_data['L_over_D']:.1f}"],
    [f"{best_data['W_total']/1000:.1f} kN"]
]

ax[1].axis('off')
ax[1].table(
    cellText=cell_text,
    rowLabels=["Wing Area", "Aspect Ratio", "Cruise Speed", "L/D Ratio", "Total Weight"],
    colLabels=["Optimal Value"],
    loc='center', 
    cellLoc='center',
    colWidths=[0.3]*2
)
ax[1].set_title("Best Configuration Parameters")

plt.tight_layout()
plt.show()