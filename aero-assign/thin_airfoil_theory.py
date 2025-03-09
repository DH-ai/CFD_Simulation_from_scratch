
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson, quad  # Fixed: changed simps to simpson
from scipy.interpolate import interp1d 

class ThinAirfoilAnalyzer:
    def __init__(self):
        self.density = 1.225  # air density in kg/m^3
        self.chord = 1.0  # normalized chord length

    def generate_camber_line(self, x_points, airfoil_type='custom', **kwargs):
        """
        Generates points along the camber line based on airfoil type
        
        Parameters:
        x_points (array): Array of x coordinates along chord
        airfoil_type (str): 'naca' or 'custom'
        kwargs: Parameters specific to airfoil type
            For NACA: m, p
            For custom: camber_function
            
        Returns:
        array: y coordinates of camber line
        """
        if airfoil_type.lower() == 'naca':
            m = kwargs.get('m', 0.02)  # maximum camber as fraction of chord
            p = kwargs.get('p', 0.4)   # location of maximum camber
            
            y_points = np.zeros_like(x_points)
            for i, x in enumerate(x_points):
                if x < p:
                    y_points[i] = m * (2 * p * x - x**2) / p**2
                else:
                    y_points[i] = m * (1 - 2 * p + 2 * p * x - x**2) / (1 - p)**2
            
            return y_points
        
        elif airfoil_type.lower() == 'custom':
            camber_function = kwargs.get('camber_function', lambda x: np.zeros_like(x))
            return camber_function(x_points)
    
    def plot_camber_line(self, x_points, y_points, title="Airfoil Camber Line"):
        """
        Plots the camber line
        
        Parameters:
        x_points (array): x coordinates
        y_points (array): y coordinates (camber line)
        title (str): Plot title
        """
        plt.figure(figsize=(10, 6))
        plt.plot(x_points, y_points, 'b-', linewidth=2)
        plt.plot(x_points, np.zeros_like(x_points), 'k--', alpha=0.3)  # chord line
        plt.grid(True, alpha=0.3)
        plt.xlabel('x/c')
        plt.ylabel('y/c')
        plt.title(title)
        plt.axis('equal')
        return plt.gcf()
    
    def compute_camber_slope(self, x_points, y_points):
        """
        Computes the slope of the camber line at given points
        
        Parameters:
        x_points (array): x coordinates
        y_points (array): y coordinates of camber line
        
        Returns:
        array: slopes at each x point
        """
        # Use central differences for interior points
        # and forward/backward differences at the endpoints
        slopes = np.zeros_like(x_points)
        
        # For interior points
        for i in range(1, len(x_points) - 1):
            slopes[i] = (y_points[i+1] - y_points[i-1]) / (x_points[i+1] - x_points[i-1])
        
        # For endpoints
        if len(x_points) > 1:
            slopes[0] = (y_points[1] - y_points[0]) / (x_points[1] - x_points[0])
            slopes[-1] = (y_points[-1] - y_points[-2]) / (x_points[-1] - x_points[-2])
        
        return slopes
    
    def plot_camber_slope(self, x_points, slopes, title="Camber Line Slope"):
        """
        Plots the slope of the camber line
        
        Parameters:
        x_points (array): x coordinates
        slopes (array): slopes at each x point
        title (str): Plot title
        """
        plt.figure(figsize=(10, 6))
        plt.plot(x_points, slopes, 'r-', linewidth=2)
        plt.grid(True, alpha=0.3)
        plt.xlabel('x/c')
        plt.ylabel('dy/dx')
        plt.title(title)
        return plt.gcf()
    
    def compute_cl(self, alpha_deg, x_points, y_points):
        """
        Computes the lift coefficient using thin airfoil theory
        
        Parameters:
        alpha_deg (float): Angle of attack in degrees
        x_points (array): x coordinates
        y_points (array): y coordinates of camber line
        
        Returns:
        float: Lift coefficient
        """
        alpha = np.radians(alpha_deg)  # Convert to radians
        
        # Calculate camber line slopes
        slopes = self.compute_camber_slope(x_points, y_points)
        
        # Thin airfoil theory: Cl = 2π(α - α0)
        # where α0 is the zero-lift angle of attack
        
        # Calculate α0 (angle where camber line is aligned with the flow)
        # For thin airfoil theory: α0 = -1/π * ∫(0 to π) [dz/dx * (1 - cos(θ))] dθ
        
        # Convert x from [0,1] to theta [0,π] using x = (1-cos(θ))/2
        theta = np.arccos(1 - 2 * x_points)
        
        # Interpolate slopes at theta points
        slope_interp = interp1d(x_points, slopes, kind='linear', fill_value="extrapolate")
        
        # Calculate integrand for α0
        integrand = slope_interp(x_points) * (1 - np.cos(theta))
        
        # Perform numerical integration
        alpha_0 = -1/np.pi * simpson(integrand, theta)
        
        # Calculate Cl using thin airfoil theory formula
        cl = 2 * np.pi * (alpha - alpha_0)
        
        return cl
    
    def compute_cl_vs_alpha(self, alpha_range, x_points, y_points):
        """
        Computes lift coefficients for a range of angles of attack
        
        Parameters:
        alpha_range (array): Array of angles of attack in degrees
        x_points (array): x coordinates
        y_points (array): y coordinates of camber line
        
        Returns:
        array: Lift coefficients for each angle of attack
        """
        cl_values = np.zeros_like(alpha_range, dtype=float)
        
        for i, alpha in enumerate(alpha_range):
            cl_values[i] = self.compute_cl(alpha, x_points, y_points)
        
        return cl_values
    
    def plot_cl_vs_alpha(self, alpha_range, cl_values, title="Lift Coefficient vs Angle of Attack", 
                         cfd_data=None, label="Thin Airfoil Theory"):
        """
        Plots the lift coefficient vs angle of attack
        
        Parameters:
        alpha_range (array): Array of angles of attack in degrees
        cl_values (array): Lift coefficients for each angle of attack
        title (str): Plot title
        cfd_data (tuple): (alpha_cfd, cl_cfd) from CFD simulations
        label (str): Label for the plot legend
        """
        plt.figure(figsize=(10, 6))
        plt.plot(alpha_range, cl_values, 'b-', linewidth=2, label=label)
        
        if cfd_data is not None:
            alpha_cfd, cl_cfd = cfd_data
            plt.plot(alpha_cfd, cl_cfd, 'ro--', linewidth=2, label='CFD Simulation')
        
        plt.grid(True, alpha=0.3)
        plt.xlabel('Angle of Attack (degrees)')
        plt.ylabel('Lift Coefficient (Cl)')
        plt.title(title)
        plt.legend()
        return plt.gcf()
    
    def calculate_circulation_distribution(self, x_points, y_points, alpha_deg):
        """
        Calculates circulation distribution along the camber line
        
        Parameters:
        x_points (array): x coordinates
        y_points (array): y coordinates of camber line
        alpha_deg (float): Angle of attack in degrees
        
        Returns:
        array: Circulation distribution at each point
        """
        alpha = np.radians(alpha_deg)
        slopes = self.compute_camber_slope(x_points, y_points)
        
        # Calculate circulation distribution
        # Γ(x) = 2πV∞c [ α(1-x/c)^(1/2) + (1/π)∫(x/c to 1) (dz/dx)/(x-ξ)^(1/2) dξ ]
        
        V_inf = 1.0  # Free stream velocity
        c = self.chord
        
        gamma = np.zeros_like(x_points)
        
        for i, x in enumerate(x_points):
            # Term 1: α(1-x/c)^(1/2)
            term1 = alpha * np.sqrt(1 - x)
            
            # Term 2: (1/π)∫(x/c to 1) (dz/dx)/(x-ξ)^(1/2) dξ
            # Simplification for numerical integration
            
            # Only integrate over points >= x
            mask = x_points >= x
            if np.any(mask):
                x_subset = x_points[mask]
                slopes_subset = slopes[mask]
                
                # Simple numerical integration (trapezoidal rule)
                term2 = 0
                for j in range(len(x_subset) - 1):
                    dx = x_subset[j+1] - x_subset[j]
                    integrand_j = slopes_subset[j] / np.sqrt(x_subset[j] - x) if x_subset[j] > x else 0
                    integrand_j1 = slopes_subset[j+1] / np.sqrt(x_subset[j+1] - x) if x_subset[j+1] > x else 0
                    term2 += 0.5 * (integrand_j + integrand_j1) * dx
                
                term2 = term2 / np.pi
            else:
                term2 = 0
            
            gamma[i] = 2 * np.pi * V_inf * c * (term1 + term2)
        
        return gamma
    
    def calculate_velocity_field(self, x_grid, y_grid, x_points, y_points, alpha_deg):
        """
        Calculates velocity field around the airfoil
        
        Parameters:
        x_grid, y_grid (array): Grid points
        x_points, y_points (array): Camber line points
        alpha_deg (float): Angle of attack in degrees
        
        Returns:
        tuple: (u, v) velocity components at each grid point
        """
        alpha = np.radians(alpha_deg)
        V_inf = 1.0  # Free stream velocity
        
        # Freestream components
        u_inf = V_inf * np.cos(alpha) * np.ones_like(x_grid)
        v_inf = V_inf * np.sin(alpha) * np.ones_like(x_grid)
        
        # Calculate circulation distribution
        gamma = self.calculate_circulation_distribution(x_points, y_points, alpha_deg)
        
        # Vortex-induced velocity
        u_vortex = np.zeros_like(x_grid)
        v_vortex = np.zeros_like(x_grid)
        
        # Place vortices along the camber line and calculate induced velocity
        for i in range(len(x_points)):
            dx = x_grid - x_points[i]
            dy = y_grid - y_points[i]
            r_squared = dx**2 + dy**2
            
            # Skip points very close to the vortex to avoid singularity
            mask = r_squared > 1e-6
            
            # Vortex strength proportional to local circulation gradient
            strength = gamma[i] / len(x_points)  # Simplified approach
            
            # Induced velocity components from a vortex
            u_vortex[mask] += -strength * dy[mask] / (2 * np.pi * r_squared[mask])
            v_vortex[mask] += strength * dx[mask] / (2 * np.pi * r_squared[mask])
        
        # Total velocity components
        u = u_inf + u_vortex
        v = v_inf + v_vortex
        
        return u, v
    
    def plot_vector_field(self, x_points, y_points, alpha_deg, domain_size=(4, 3), 
                          grid_density=30, title="Velocity Field"):
        """
        Plots vector field around the airfoil
        
        Parameters:
        x_points, y_points (array): Camber line points
        alpha_deg (float): Angle of attack in degrees
        domain_size (tuple): Size of domain in multiples of chord length (x, y)
        grid_density (int): Number of grid points in each direction
        title (str): Plot title
        """
        # Create grid
        x_size, y_size = domain_size
        x_min, x_max = -x_size/4, x_size*3/4
        y_min, y_max = -y_size/2, y_size/2
        
        x = np.linspace(x_min, x_max, grid_density)
        y = np.linspace(y_min, y_max, grid_density)
        
        X, Y = np.meshgrid(x, y)
        
        # Calculate velocity field
        u, v = self.calculate_velocity_field(X, Y, x_points, y_points, alpha_deg)
        
        # Calculate speed for color mapping
        speed = np.sqrt(u**2 + v**2)
        
        # Plot
        plt.figure(figsize=(12, 8))
        plt.streamplot(X, Y, u, v, density=1.5, color=speed, cmap='viridis', linewidth=1)
        plt.quiver(X, Y, u, v, scale=25, width=0.002, color='gray', alpha=0.7)
        
        # Plot the airfoil
        plt.plot(x_points, y_points, 'k-', linewidth=2)
        
        plt.colorbar(label='Velocity Magnitude')
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.xlabel('x/c')
        plt.ylabel('y/c')
        plt.title(title)
        
        return plt.gcf(), (u, v, X, Y)
    
    def calculate_circulation_line_integral(self, u, v, X, Y):
        """
        Calculates circulation around airfoil using velocity line integral
        
        Parameters:
        u, v (array): Velocity components
        X, Y (array): Grid coordinates
        
        Returns:
        float: Circulation value
        """
        # Define a contour around the airfoil
        R = 1.5  # Radius of the circular contour
        
        theta = np.linspace(0, 2*np.pi, 1000)
        x_contour = R * np.cos(theta)
        y_contour = R * np.sin(theta)
        
        # Fixed: Proper handling of interpolation with correct dimensions
        # Create interpolation functions for u and v over the grid
        # First flatten the arrays to use griddata
        x_flat = X.flatten()
        y_flat = Y.flatten()
        u_flat = u.flatten()
        v_flat = v.flatten()
        
        # Use a simpler approach with nearest neighbor interpolation
        from scipy.interpolate import NearestNDInterpolator
        
        u_interp = NearestNDInterpolator(list(zip(x_flat, y_flat)), u_flat)
        v_interp = NearestNDInterpolator(list(zip(x_flat, y_flat)), v_flat)
        
        # Get interpolated values at contour points
        u_contour = u_interp(x_contour, y_contour)
        v_contour = v_interp(x_contour, y_contour)
        
        # Calculate line integral around the contour
        # Γ = ∮ v⋅dl = ∮ (u dx + v dy)
        
        # Calculate dl components
        dx = np.diff(np.append(x_contour, x_contour[0]))
        dy = np.diff(np.append(y_contour, y_contour[0]))
        
        # Calculate circulation using line integral
        circulation = np.sum(u_contour * dx + v_contour * dy)
        
        return circulation
    
    def calculate_bound_circulation(self, x_points, gamma):
        """
        Calculates bound circulation by integrating circulation distribution along camber line
        
        Parameters:
        x_points (array): x coordinates
        gamma (array): Circulation distribution
        
        Returns:
        float: Bound circulation
        """
        # Total bound circulation is the sum of all vortex strengths
        bound_circulation = simpson(gamma, x_points)
        
        return bound_circulation

    def main(self):
        """
        Main function to run the program
        """
        print("Thin Airfoil Theory Analysis")
        print("---------------------------")
        
        # User input
        airfoil_type = input("Enter airfoil type (NACA or custom): ").lower()
        
        if airfoil_type == "naca":
            # NACA input
            naca_code = input("Enter NACA 4-digit code: ")
            m = int(naca_code[0]) / 100
            p = int(naca_code[1]) / 10
            
            # Generate points
            x_points = np.linspace(0, 1, 100)
            y_points = self.generate_camber_line(x_points, airfoil_type='naca', m=m, p=p)
            
        elif airfoil_type == "custom":
            # Custom input - simple polynomial or predefined function
            print("Define the camber line as y = f(x)")
            print("1. Parabolic camber")
            print("2. Circular arc camber")
            print("3. S-shaped camber")
            print("4. Custom polynomial")
            
            choice = int(input("Enter your choice (1-4): "))
            
            x_points = np.linspace(0, 1, 100)
            
            if choice == 1:
                # Parabolic camber
                max_camber = float(input("Enter maximum camber (fraction of chord): "))
                pos_max_camber = float(input("Enter position of maximum camber (fraction of chord): "))
                
                def camber_function(x):
                    y = np.zeros_like(x)
                    for i, xi in enumerate(x):
                        if xi < pos_max_camber:
                            y[i] = max_camber * (2 * pos_max_camber * xi - xi**2) / pos_max_camber**2
                        else:
                            y[i] = max_camber * (1 - 2 * pos_max_camber + 2 * pos_max_camber * xi - xi**2) / (1 - pos_max_camber)**2
                    return y
                
            elif choice == 2:
                # Circular arc camber
                max_camber = float(input("Enter maximum camber (fraction of chord): "))
                
                def camber_function(x):
                    # Circular arc with maximum camber at x = 0.5
                    return max_camber * np.sin(np.pi * x)
                
            elif choice == 3:
                # S-shaped camber
                max_camber = float(input("Enter maximum camber (fraction of chord): "))
                
                def camber_function(x):
                    # S-shaped camber using sin function
                    return max_camber * np.sin(2 * np.pi * x)
                
            elif choice == 4:
                # Custom polynomial
                coefs_str = input("Enter polynomial coefficients (space-separated, highest degree first): ")
                coefs = [float(c) for c in coefs_str.split()]
                
                def camber_function(x):
                    y = np.zeros_like(x)
                    for i, c in enumerate(coefs):
                        y += c * x**(len(coefs) - i - 1)
                    return y
            
            y_points = camber_function(x_points)
        
        else:
            print("Invalid airfoil type")
            return
        
        # Calculate and plot camber line
        self.plot_camber_line(x_points, y_points, title=f"{airfoil_type.upper()} Airfoil Camber Line")
        
        # Calculate and plot camber line slope
        slopes = self.compute_camber_slope(x_points, y_points)
        self.plot_camber_slope(x_points, slopes, title=f"{airfoil_type.upper()} Airfoil Camber Line Slope")
        
        # Calculate and plot lift coefficient vs angle of attack
        alpha_range = np.linspace(-3, 12, 16)
        cl_values = self.compute_cl_vs_alpha(alpha_range, x_points, y_points)
        self.plot_cl_vs_alpha(alpha_range, cl_values, title=f"{airfoil_type.upper()} Airfoil Cl vs Alpha")
        
        # Calculate and plot vector field
        alpha_test = 3  # Test at 3 degrees
        vf_plot, (u, v, X, Y) = self.plot_vector_field(x_points, y_points, alpha_test, 
                                                      title=f"{airfoil_type.upper()} Airfoil Velocity Field at α={alpha_test}°")
        
        # Calculate circulation using line integral
        circulation_line = self.calculate_circulation_line_integral(u, v, X, Y)
        print(f"Circulation (line integral): {circulation_line:.4f}")
        
        # Calculate circulation using bound vorticity
        gamma = self.calculate_circulation_distribution(x_points, y_points, alpha_test)
        circulation_bound = self.calculate_bound_circulation(x_points, gamma)
        print(f"Circulation (bound vorticity): {circulation_bound:.4f}")
        
        plt.savefig('plot.png')
        

# Example usage of the class
if __name__ == "__main__":
    analyzer = ThinAirfoilAnalyzer()
    analyzer.main()


