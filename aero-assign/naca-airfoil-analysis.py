import numpy as np
import matplotlib.pyplot as plt
from thin_airfoil_theory import ThinAirfoilAnalyzer

# Create analyzer instance
analyzer = ThinAirfoilAnalyzer()

# Create points along the chord
x_points = np.linspace(0, 1, 100)

# NACA Airfoil (assuming NACA 4412 from Assignment 1, modify if needed)
m = 0.04  # maximum camber
p = 0.4   # location of maximum camber

# Generate camber line
y_points = analyzer.generate_camber_line(x_points, airfoil_type='naca', m=m, p=p)

# 1. Plot camber line
camber_plot = analyzer.plot_camber_line(x_points, y_points, title="NACA 4412 Camber Line")
camber_plot.savefig("naca_camber_line.png", dpi=300, bbox_inches='tight')

# 2. Calculate and plot camber line slope
slopes = analyzer.compute_camber_slope(x_points, y_points)
slope_plot = analyzer.plot_camber_slope(x_points, slopes, title="NACA 4412 Camber Line Slope")
slope_plot.savefig("naca_camber_slope.png", dpi=300, bbox_inches='tight')

# 3. Calculate Cl vs alpha and compare with CFD
alpha_range = np.linspace(-3, 12, 16)
cl_values = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_points)

# CFD data from Assignment 1 (placeholder - replace with your actual data)
# Format: alpha_cfd (angles in degrees), cl_cfd (lift coefficients)
alpha_cfd = np.array([-3, 0, 3, 6, 9, 12])
cl_cfd = np.array([0.1, 0.4, 0.7, 1.0, 1.2, 1.3])  # Replace with your actual CFD results

# Plot Cl vs alpha comparison
cl_plot = analyzer.plot_cl_vs_alpha(alpha_range, cl_values, 
                                   title="NACA 4412: Cl vs Alpha Comparison",
                                   cfd_data=(alpha_cfd, cl_cfd),
                                   label="Thin Airfoil Theory")
cl_plot.savefig("naca_cl_vs_alpha.png", dpi=300, bbox_inches='tight')

# 4. Calculate and plot vector field at alpha = 3 degrees
vf_plot, (u, v, X, Y) = analyzer.plot_vector_field(x_points, y_points, 3, 
                                                  title="NACA 4412: Velocity Field at α=3°")
vf_plot.savefig("naca_vector_field.png", dpi=300, bbox_inches='tight')

# 5. Calculate circulation using line integral
circulation_line = analyzer.calculate_circulation_line_integral(u, v, X, Y)
print(f"Circulation (line integral): {circulation_line:.4f}")

# 6. Calculate bound circulation
gamma = analyzer.calculate_circulation_distribution(x_points, y_points, 3)
circulation_bound = analyzer.calculate_bound_circulation(x_points, gamma)
print(f"Circulation (bound vorticity): {circulation_bound:.4f}")

print("Analysis complete!")
