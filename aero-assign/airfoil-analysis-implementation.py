import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

# Implementation of all thin airfoil theory functions from scratch
class ThinAirfoilTheory:
    def __init__(self):
        """Initialize the ThinAirfoilTheory class with default parameters."""
        self.rho = 1.225  # air density in kg/m^3
        self.V_inf = 1.0  # freestream velocity in m/s
        
    def generate_camber_line(self, x_points, airfoil_type='naca', m=0.04, p=0.4):
        """
        Generate camber line for a NACA 4-digit airfoil.
        
        Args:
            x_points: array of x coordinates along chord (0 to 1)
            airfoil_type: 'naca' for NACA 4-digit
            m: maximum camber
            p: location of maximum camber
            
        Returns:
            y_points: array of camber line y coordinates
        """
        y_points = np.zeros_like(x_points)
        
        if airfoil_type.lower() == 'naca':
            for i, x in enumerate(x_points):
                if x < p:
                    # Formula for camber line before p
                    y_points[i] = m * (2 * p * x - x**2) / p**2
                else:
                    # Formula for camber line after p
                    y_points[i] = m * (1 - 2 * p + 2 * p * x - x**2) / (1 - p)**2
                    
        return y_points
    
    def compute_camber_slope(self, x_points, y_points):
        """
        Calculate the slope of the camber line using central differences.
        
        Args:
            x_points: array of x coordinates
            y_points: array of y coordinates
            
        Returns:
            slopes: array of slopes at each point
        """
        # Use numpy gradient function for central differences
        slopes = np.gradient(y_points, x_points)
        return slopes
    
    def plot_camber_line(self, x_points, y_points, title="Camber Line"):
        """
        Plot the camber line.
        
        Args:
            x_points: array of x coordinates
            y_points: array of y coordinates
            title: plot title
            
        Returns:
            fig: matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(x_points, y_points, 'b-', linewidth=2)
        ax.set_title(title)
        ax.set_xlabel('x/c')
        ax.set_ylabel('y/c')
        ax.grid(True)
        ax.set_aspect('equal')
        
        # Set y-axis limits to better visualize camber (which is typically small)
        y_max = max(abs(np.max(y_points)), abs(np.min(y_points)))
        ax.set_ylim(-2*y_max, 2*y_max)
        
        return fig
    
    def plot_camber_slope(self, x_points, slopes, title="Camber Line Slope"):
        """
        Plot the slope of the camber line.
        
        Args:
            x_points: array of x coordinates
            slopes: array of slopes
            title: plot title
            
        Returns:
            fig: matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(x_points, slopes, 'r-', linewidth=2)
        ax.set_title(title)
        ax.set_xlabel('x/c')
        ax.set_ylabel('Slope (dy/dx)')
        ax.grid(True)
        
        return fig
    
    def compute_cl_vs_alpha(self, alpha_range, x_points, y_points):
        """
        Calculate lift coefficient across a range of angles of attack.
        Based on thin airfoil theory.
        
        Args:
            alpha_range: array of angles of attack in degrees
            x_points: array of x coordinates
            y_points: array of y coordinates
            
        Returns:
            cl_values: array of lift coefficients
        """
        # Calculate camber slope
        slopes = self.compute_camber_slope(x_points, y_points)
        
        # Calculate the zero-lift angle of attack (alpha_0)
        # In thin airfoil theory, alpha_0 = -1 * arctangent(slope at quarter-chord)
        quarter_chord_idx = np.abs(x_points - 0.25).argmin()
        alpha_0 = -np.arctan(slopes[quarter_chord_idx]) * 180 / np.pi
        
        # Calculate lift slope (2π per radian in thin airfoil theory)
        cl_slope = 2 * np.pi / 180  # Convert to per degree
        
        # Calculate lift coefficient for each angle of attack
        cl_values = cl_slope * (alpha_range - alpha_0)
        
        return cl_values
    
    def plot_cl_vs_alpha(self, alpha_range, cl_values, title="Cl vs Alpha", 
                         cfd_data=None, label="Theory"):
        """
        Plot lift coefficient vs angle of attack.
        
        Args:
            alpha_range: array of angles of attack
            cl_values: array of lift coefficients
            title: plot title
            cfd_data: tuple (alpha_cfd, cl_cfd) for comparison
            label: label for theory data
            
        Returns:
            fig: matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(alpha_range, cl_values, 'b-', linewidth=2, label=label)
        
        if cfd_data is not None:
            alpha_cfd, cl_cfd = cfd_data
            ax.scatter(alpha_cfd, cl_cfd, color='red', s=50, marker='o', label='CFD Data')
        
        ax.set_title(title)
        ax.set_xlabel('Angle of Attack (degrees)')
        ax.set_ylabel('Lift Coefficient (Cl)')
        ax.grid(True)
        ax.legend()
        
        return fig
    
    def calculate_velocity_field(self, x_points, y_points, alpha_deg, X, Y):
        """
        Calculate velocity field around an airfoil using thin airfoil theory.
        
        Args:
            x_points: array of x coordinates of camber line
            y_points: array of y coordinates of camber line
            alpha_deg: angle of attack in degrees
            X, Y: meshgrid for vector field evaluation
            
        Returns:
            u, v: velocity components at each point in X, Y
        """
        # Convert angle of attack to radians
        alpha = alpha_deg * np.pi / 180
        
        # Calculate camber slope
        slopes = self.compute_camber_slope(x_points, y_points)
        
        # Initialize velocity arrays
        u = np.ones_like(X) * self.V_inf * np.cos(alpha)
        v = np.ones_like(Y) * self.V_inf * np.sin(alpha)
        
        # For each point on the airfoil, add the effect of a vortex panel
        for i in range(len(x_points)-1):
            x1, y1 = x_points[i], y_points[i]
            x2, y2 = x_points[i+1], y_points[i+1]
            
            # Panel length
            panel_length = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            
            # Panel angle
            panel_angle = np.arctan2(y2-y1, x2-x1)
            
            # Average slope for this panel
            avg_slope = (slopes[i] + slopes[i+1]) / 2
            
            # Circulation strength (proportional to local angle of attack)
            gamma = 2 * self.V_inf * (alpha + np.arctan(avg_slope))
            
            # For each point in the grid, add vortex influence
            for j in range(X.shape[0]):
                for k in range(X.shape[1]):
                    # Distance from panel to point
                    dx = X[j,k] - (x1 + x2)/2
                    dy = Y[j,k] - (y1 + y2)/2
                    r_squared = dx**2 + dy**2
                    
                    if r_squared > 0.001:  # Avoid singularity
                        # Induced velocity by a vortex
                        u_ind = -gamma * dy / (2 * np.pi * r_squared)
                        v_ind = gamma * dx / (2 * np.pi * r_squared)
                        
                        # Add to velocity field
                        u[j,k] += u_ind * panel_length
                        v[j,k] += v_ind * panel_length
        
        return u, v
    
    def plot_vector_field(self, x_points, y_points, alpha_deg, title="Velocity Field"):
        """
        Plot velocity vector field around an airfoil.
        
        Args:
            x_points: array of x coordinates of camber line
            y_points: array of y coordinates of camber line
            alpha_deg: angle of attack in degrees
            title: plot title
            
        Returns:
            fig: matplotlib figure object
            (u, v, X, Y): velocity components and meshgrid
        """
        # Create a coarser grid for vector field
        x_grid = np.linspace(-0.5, 1.5, 30)
        y_grid = np.linspace(-0.5, 0.5, 20)
        X, Y = np.meshgrid(x_grid, y_grid)
        
        # Calculate velocity field
        u, v = self.calculate_velocity_field(x_points, y_points, alpha_deg, X, Y)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot camber line
        ax.plot(x_points, y_points, 'k-', linewidth=2)
        
        # Scale vectors for better visualization
        magnitude = np.sqrt(u**2 + v**2)
        scale = 0.8 * min(np.diff(x_grid)[0], np.diff(y_grid)[0]) / np.max(magnitude)
        
        # Plot vector field
        q = ax.quiver(X, Y, u, v, magnitude, scale=1/scale, cmap='viridis', width=0.002)
        
        # Adding a color bar 
        cbar = fig.colorbar(q, ax=ax)
        cbar.set_label('Velocity Magnitude (m/s)')

        ax.set_title(title)
        ax.set_xlabel('x/c')
        ax.set_ylabel('y/c')
        ax.set_aspect('equal')
        ax.grid(True)
        
        return fig, (u, v, X, Y)
    
    def calculate_circulation_line_integral(self, u, v, X, Y):
        """
        Calculate circulation around the airfoil using line integral.
        
        Args:
            u, v: velocity components at each point
            X, Y: meshgrid on which u, v are defined
            
        Returns:
            circulation: calculated circulation
        """
        # Create a circular path around the airfoil
        theta = np.linspace(0, 2*np.pi, 100)
        radius = 0.4
        center_x, center_y = 0.5, 0.0
        
        path_x = center_x + radius * np.cos(theta)
        path_y = center_y + radius * np.sin(theta)
        
        # Interpolate velocity at path points
        u_path = np.zeros_like(path_x)
        v_path = np.zeros_like(path_y)
        
        for i in range(len(path_x)):
            # Find closest grid point
            idx_x = np.abs(X[0,:] - path_x[i]).argmin()
            idx_y = np.abs(Y[:,0] - path_y[i]).argmin()
            
            u_path[i] = u[idx_y, idx_x]
            v_path[i] = v[idx_y, idx_x]
        
        # Calculate tangential velocity
        dx = np.gradient(path_x)
        dy = np.gradient(path_y)
        ds = np.sqrt(dx**2 + dy**2)
        
        tangential_velocity = (u_path * dy - v_path * dx) / ds
        
        # Calculate circulation
        circulation = np.sum(tangential_velocity * ds)
        
        return circulation
    
    def calculate_circulation_distribution(self, x_points, y_points, alpha_deg):
        """
        Calculate the distribution of bound vorticity along the airfoil.
        
        Args:
            x_points: array of x coordinates
            y_points: array of y coordinates
            alpha_deg: angle of attack in degrees
            
        Returns:
            gamma: array of vorticity strengths
        """
        # Convert to radians
        alpha = alpha_deg * np.pi / 180
        
        # Calculate camber slope
        slopes = self.compute_camber_slope(x_points, y_points)
        
        # Calculate local angles of attack
        local_alpha = alpha + np.arctan(slopes)
        
        # Calculate vorticity distribution (proportional to local angle of attack)
        gamma = 2 * self.V_inf * local_alpha
        
        return gamma
    
    def calculate_bound_circulation(self, x_points, gamma):
        """
        Calculate total bound circulation from vorticity distribution.
        
        Args:
            x_points: array of x coordinates
            gamma: array of vorticity strengths
            
        Returns:
            circulation: total circulation
        """
        # Integrate gamma along the chord
        dx = np.diff(x_points)
        gamma_avg = (gamma[:-1] + gamma[1:]) / 2
        
        circulation = np.sum(gamma_avg * dx)
        
        return circulation

# Main code to run the analysis
if __name__ == "__main__":
    # Create analyzer instance
    analyzer = ThinAirfoilTheory()

    # Create points along the chord
    x_points = np.linspace(0, 1, 100)

    # NACA Airfoil (NACA 4412)
    # NACA 0011 for 23b00007
    # m = 0 # maximum camber
    # p = 15  # location of maximum camber

    m = 6.9# maximum camber
    p =  25.5 # location of maximum camber


    # Generate camber line
    y_points = analyzer.generate_camber_line(x_points, airfoil_type='naca', m=m, p=p)

    # 1. Plot camber line
    camber_plot = analyzer.plot_camber_line(x_points, y_points, title="Camber Line")
    camber_plot.savefig("naca_camber_line.png", dpi=300, bbox_inches='tight')

    # 2. Calculate and plot camber line slope
    slopes = analyzer.compute_camber_slope(x_points, y_points)
    slope_plot = analyzer.plot_camber_slope(x_points, slopes, title="Camber Line Slope")
    slope_plot.savefig("naca_camber_slope.png", dpi=300, bbox_inches='tight')

    # 3. Calculate Cl vs alpha and compare with CFD
    alpha_range = np.linspace(-3, 12, 16)
    cl_values = analyzer.compute_cl_vs_alpha(alpha_range, x_points, y_points)

    # CFD data (placeholder - replace with actual data)
    alpha_cfd = np.array([-3, 0, 3, 6, 9, 12])
    cl_cfd = np.array([0.1, 0.4, 0.7, 1.0, 1.2, 1.3])

    # Plot Cl vs alpha comparison
    cl_plot = analyzer.plot_cl_vs_alpha(alpha_range, cl_values, 
                                      title="Cl vs Alpha Comparison",
                                      cfd_data=(alpha_cfd, cl_cfd),
                                      label="Thin Airfoil Theory")
    cl_plot.savefig("naca_cl_vs_alpha.png", dpi=300, bbox_inches='tight')

    # 4. Calculate and plot vector field at alpha = 3 degrees
    vf_plot, (u, v, X, Y) = analyzer.plot_vector_field(x_points, y_points, 3, 
                                                    title="Velocity Field at α=3°")
    vf_plot.savefig("naca_vector_field.png", dpi=300, bbox_inches='tight')

    # 5. Calculate circulation using line integral
    circulation_line = analyzer.calculate_circulation_line_integral(u, v, X, Y)
    print(f"Circulation (line integral): {circulation_line:.4f}")

    # 6. Calculate bound circulation
    gamma = analyzer.calculate_circulation_distribution(x_points, y_points, 3)
    circulation_bound = analyzer.calculate_bound_circulation(x_points, gamma)
    print(f"Circulation (bound vorticity): {circulation_bound:.4f}")

    print("Analysis complete!")
