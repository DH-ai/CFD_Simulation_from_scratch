import numpy as np
import matplotlib.pyplot as plt
from thin_airfoil_theory import ThinAirfoilAnalyzer

# Create analyzer instance
analyzer = ThinAirfoilAnalyzer()

# Create points along the chord
x_points = np.linspace(0, 1, 100)

# NACA Airfoil parameters (for comparison)
m_naca = 0.04  # maximum camber
p_naca = 0.4   # location of maximum camber

# Generate NACA camber line
y_naca = analyzer.generate_camber_line(x_points, airfoil_type='naca', m=m_naca, p=p_naca)

# Define three custom airfoils
# Custom Airfoil 1: Reflexed camber line
def custom_airfoil1(x):
    # Reflexed camber line with maximum camber at 0.3c and reflex at 0.7c
    max_camber = 0.06
    return max_camber * (3*(x/0.3)**2 - 2*(x/0.3)**3) * (x <= 0.3) + \
           max_camber * (1 - 1.5*((x-0.3)/0.7) + 0.5*((x-0.3)/0.7)**3) * (x > 0.3)

# Custom Airfoil 2: Front-loaded camber line
def custom_airfoil2(x):
    # Front-loaded camber line with maximum camber at 0.2c
    max_camber = 0.05
    return max_camber * np.sin(np.pi*x/0.4) * (x <= 0.4) + \
           max_camber * np.sin(np.pi*0.4/0.4) * (1 - (x - 0.4)/(1 - 0.4)) * (x > 0.4)

# Custom Airfoil 3: S-shaped camber line
def custom_airfoil3(x):
    # S-shaped camber with positive camber in front, negative in rear
    max_camber = 0.035
    return max_camber * np.sin(2 * np.pi * x)

# Generate camber lines for custom airfoils
y_custom1 = custom_airfoil1(x_points)
y_custom2 = custom_airfoil2(x_points)
y_custom3 = custom_airfoil3(x_points)

# 1. Plot camber lines for comparison
plt.figure(figsize=(10, 6))
plt.plot(x_points, y_naca, 'k-', linewidth=2, label='NACA 4412')
plt.plot(x_points, y_custom1, 'b-', linewidth=2, label='Custom 1: Reflexed')
plt.plot(x_points, y_custom2, 'r-', linewidth=2, label='Custom 2: Front-loaded')
plt.plot(x_points, y_custom3, 'g-', linewidth=2, label='Custom 3: S-shaped')
plt.grid(True, alpha=0.3)
plt.xlabel('x/c')
plt.ylabel('y/c')
plt.title('Comparison of Airfoil Camber Lines')
plt.legend()
plt.axis('equal')
plt.savefig("custom_camber_comparison.png", dpi=300, bbox_inches='tight')

# 2. Calculate Cl vs alpha for all airfoils
alpha_range = np.linspace(-3, 12, 16)

cl_naca = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_naca)
cl_custom1 = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_custom1)
cl_custom2 = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_custom2)
cl_custom3 = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_custom3)

# Plot Cl vs alpha comparison
plt.figure(figsize=(10, 6))
plt.plot(alpha_range, cl_naca, 'k-', linewidth=2, label='NACA 4412')
plt.plot(alpha_range, cl_custom1, 'b-', linewidth=2, label='Custom 1: Reflexed')
plt.plot(alpha_range, cl_custom2, 'r-', linewidth=2, label='Custom 2: Front-loaded')
plt.plot(alpha_range, cl_custom3, 'g-', linewidth=2, label='Custom 3: S-shaped')
plt.grid(True, alpha=0.3)
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient (Cl)')
plt.title('Comparison of Cl vs Alpha')
plt.legend()
plt.savefig("custom_cl_comparison.png", dpi=300, bbox_inches='tight')

# 3. Calculate and plot vector fields at alpha = 3 degrees
# Custom Airfoil 1
vf_plot1, _ = analyzer.plot_vector_field(x_points, y_custom1, 3, 
                                       title="Custom 1 (Reflexed): Velocity Field at α=3°")
vf_plot1.savefig("custom1_vector_field.png", dpi=300, bbox_inches='tight')

# Custom Airfoil 2
vf_plot2, _ = analyzer.plot_vector_field(x_points, y_custom2, 3, 
                                       title="Custom 2 (Front-loaded): Velocity Field at α=3°")
vf_plot2.savefig("custom2_vector_field.png", dpi=300, bbox_inches='tight')

# Custom Airfoil 3
vf_plot3, _ = analyzer.plot_vector_field(x_points, y_custom3, 3, 
                                       title="Custom 3 (S-shaped): Velocity Field at α=3°")
vf_plot3.savefig("custom3_vector_field.png", dpi=300, bbox_inches='tight')

print("Custom airfoil analysis complete!")
