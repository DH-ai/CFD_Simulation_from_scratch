# Thin Airfoil Theory Implementation
# This program implements thin airfoil theory to analyze airfoil performance
# Author: Claude Assistant

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson, quad  # Fixed: changed simps to simpson
from scipy.interpolate import interp1d

class ThinAirfoilAnalyzer:
    def __init__(self):
        self.density = 1.225        # Air density in kg/m^3
        self.chord = 1.0            # Normalized chord length
        self.V_inf = 1.0            # Freestream velocity
        # Define chordwise x-points (from 0 to 1)
        self.x_points = np.linspace(0, 1, 100)

    def generate_camber_line(self, airfoil_type, **kwargs):
        """
        Generate camber line coordinates based on airfoil type.
        
        Parameters:
          airfoil_type: 'naca', 'custom', or 'tabulated'
          For 'naca': expects m (max camber fraction) and p (location of max camber)
          For 'custom': expects a callable function (camber_function)
          For 'tabulated': expects arrays tab_x and tab_y (will be interpolated)
          
        Returns:
          y: Array of camber line y-coordinates.
        """
        if airfoil_type.lower() == 'naca':
            m = kwargs.get('m', 0.02)
            p = kwargs.get('p', 0.4)
            y = np.zeros_like(self.x_points)
            mask = self.x_points < p
            # For x < p use first formula
            y[mask] = m * (2 * p * self.x_points[mask] - self.x_points[mask]**2) / p**2
            # For x >= p use second formula
            y[~mask] = m * (1 - 2*p + 2*p*self.x_points[~mask] - self.x_points[~mask]**2) / (1-p)**2
            return y
        elif airfoil_type.lower() == 'custom':
            camber_function = kwargs.get('camber_function', lambda x: np.zeros_like(x))
            return camber_function(self.x_points)
        elif airfoil_type.lower() == 'tabulated':
            tab_x = kwargs['tab_x']
            tab_y = kwargs['tab_y']
            f = interp1d(tab_x, tab_y, bounds_error=False, fill_value="extrapolate")
            return f(self.x_points)
        else:
            return np.zeros_like(self.x_points)

    def compute_camber_slope(self, y_points):
        """
        Compute the slope (dy/dx) of the camber line using finite differences.
        
        Parameters:
          y_points: Array of camber line y-values.
          
        Returns:
          dy: Array of slopes at each x-point.
        """
        x = self.x_points
        dy = np.zeros_like(x)
        # Central differences for interior points:
        dy[1:-1] = (y_points[2:] - y_points[:-2]) / (x[2:] - x[:-2])
        # Forward difference for the first point:
        dy[0] = (y_points[1] - y_points[0]) / (x[1] - x[0])
        # Backward difference for the last point:
        dy[-1] = (y_points[-1] - y_points[-2]) / (x[-1] - x[-2])
        return dy

    def compute_cl(self, alpha_deg, y_points):
        """
        Compute the lift coefficient (Cl) using thin-airfoil theory:
          Cl = 2π (α - α0)
        where α0 is determined via an integral of the camber slope.
        
        Parameters:
          alpha_deg: Angle of attack in degrees.
          y_points: Camber line coordinates.
          
        Returns:
          Cl: Lift coefficient.
        """
        alpha = np.radians(alpha_deg)
        dy_dx = self.compute_camber_slope(y_points)
        # Convert x from [0,1] to θ in [0,π] using: x = (1-cosθ)/2 => θ = arccos(1-2x)
        theta = np.arccos(1 - 2 * self.x_points)
        # Compute the integrand: dy/dx * (1-cosθ)
        integrand = dy_dx * (1 - np.cos(theta))
        alpha_0 = -simpson(integrand, theta) / np.pi
        Cl = 2 * np.pi * (alpha - alpha_0)
        return Cl

    def calculate_cm(self, alpha_deg, y_points):
        """
        Compute the moment coefficient about the aerodynamic center (≈ quarter chord)
        using a Fourier-series expansion of the camber line.
        
        Steps:
          1. Reparameterize x using: x = (1 - cosθ)/2, θ ∈ [0,π]
          2. Compute Fourier sine coefficients:
             A_n = (2/π)*∫[0,π] z((1-cosθ)/2) sin(nθ) dθ.
          3. Then, Cm_ac = -π * sum_{n=1}^{N} [ n * A_n / (n+1) ].
        
        Returns:
          Cm_ac: Moment coefficient about the aerodynamic center.
          
        Note: In thin-airfoil theory, Cm_ac is independent of α.
        """
        N = 10  # number of Fourier terms to use
        # Use self.x_points (if available) as x domain; else generate a linear space.
        if hasattr(self, "x_points"):
            x_pts = self.x_points
        else:
            x_pts = np.linspace(0, 1, len(y_points))
        # Create a fine grid in θ.
        theta = np.linspace(0, np.pi, 200)
        # Reparameterize: x = (1 - cosθ) / 2.
        x_reparam = (1 - np.cos(theta)) / 2
        # Interpolate the camber line onto this new x.
        z = np.interp(x_reparam, x_pts, y_points)
        A = np.zeros(N)
        for n in range(1, N+1):
            integrand = z * np.sin(n * theta)
            A[n-1] = (2/np.pi) * simpson(integrand, theta)
        Cm_ac = -np.pi * np.sum([ (n * A[n-1]) / (n+1) for n in range(1, N+1) ])
        return Cm_ac

    def calculate_velocity_field(self, y_points, alpha_deg=3):
        """
        Compute the velocity field around the airfoil using a simplified vortex panel method.
        Domain: x from -1 to 3, y from -1.5 to 1.5.
        Uses a reduced grid (25×25) for faster computation.
        
        Parameters:
          alpha_deg: Angle of attack (°).
          y_points: Camber line coordinates.
          
        Returns:
          u, v: Velocity components.
          X, Y: Meshgrid arrays.
        """
        alpha = np.radians(alpha_deg)
        # Reduced grid size for quick computation.
        x_grid = np.linspace(-1, 3, 25)
        y_grid = np.linspace(-1.5, 1.5, 25)
        X, Y = np.meshgrid(x_grid, y_grid)
        # Freestream velocity.
        u = np.ones_like(X) * self.V_inf * np.cos(alpha)
        v = np.ones_like(X) * self.V_inf * np.sin(alpha)

        # Compute circulation distribution based on camber slope.
        slopes = self.compute_camber_slope(y_points)
        gamma = 2 * self.V_inf * (alpha + np.arctan(slopes))
        
        # Loop over each panel along the camber line.
        for i in range(len(self.x_points) - 1):
            x1, y1 = self.x_points[i], y_points[i]
            x2, y2 = self.x_points[i+1], y_points[i+1]
            panel_length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            panel_center_x = (x1 + x2) / 2
            panel_center_y = (y1 + y2) / 2
            for j in range(X.shape[0]):
                for k in range(X.shape[1]):
                    dx = X[j, k] - panel_center_x
                    dy = Y[j, k] - panel_center_y
                    r_sq = dx**2 + dy**2
                    if r_sq > 1e-3:
                        u_ind = -gamma[i] * dy / (2 * np.pi * r_sq)
                        v_ind = gamma[i] * dx / (2 * np.pi * r_sq)
                        u[j, k] += u_ind * panel_length
                        v[j, k] += v_ind * panel_length
        return u, v, X, Y

    def plot_vector_field(self, y_points, alpha_deg=3):
        """
        Plot velocity vectors and streamlines.
        Masks a narrow band around the camber line (in the central region) to avoid streamlines crossing it.
        
        Parameters:
          alpha_deg: Angle of attack (°).
          y_points: Camber line coordinates.
          
        Returns:
          u, v, X, Y: Velocity field and grid arrays.
        """
        u, v, X, Y = self.calculate_velocity_field(y_points, alpha_deg)
        speed = np.sqrt(u**2 + v**2)
        
        # Mask a narrow band around the camber line in the central region.
        tol = 0.0002
        f_interp = interp1d(self.x_points, y_points, bounds_error=False, fill_value="extrapolate")
        mask = (np.abs(Y - f_interp(X)) < tol) & (X > 0.1) & (X < 0.9)
        u[mask] = np.nan
        v[mask] = np.nan
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.quiver(X, Y, u, v, scale=50, width=0.002, color='gray')
        ax.streamplot(X, Y, u, v, density=2, color=speed, cmap='jet', linewidth=1)
        sm = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=np.nanmin(speed), vmax=np.nanmax(speed)))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Velocity Magnitude")
        ax.plot(self.x_points, y_points, 'k-', linewidth=2)
        ax.set_title(f"Velocity Field at α = {alpha_deg}°\n(Velocity magnitude shown by colour)")
        ax.set_xlabel('x/c')
        ax.set_ylabel('y/c')
        ax.axis('equal')
        plt.show()
        return u, v, X, Y

# Example usage of the class
if __name__ == "__main__":
    analyzer = ThinAirfoilAnalyzer()
    analyzer.main()