import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson, quad
from scipy.interpolate import interp1d, NearestNDInterpolator

# =============================================================================
# Unified ThinAirfoilAnalyzer Class
# =============================================================================
class ThinAirfoilAnalyzer:
    def __init__(self):
        self.density = 1.225        # Air density (kg/m^3)
        self.chord = 1.0            # Normalized chord length
        self.V_inf = 1.0            # Freestream velocity (magnitude)
        # Chordwise x-points from leading edge (0) to trailing edge (1)
        self.x_points = np.linspace(0, 1, 100)

    def generate_camber_line(self, airfoil_type, **kwargs):
        """
        Generate camber line coordinates based on airfoil type.
        For 'naca': expects m (max camber fraction) and p (location of max camber).
        For 'custom': expects a callable function (camber_function).
        For 'tabulated': expects arrays tab_x and tab_y (will be interpolated).
        """
        if airfoil_type.lower() == 'naca':
            m = kwargs.get('m', 0.02)
            p = kwargs.get('p', 0.4)
            y = np.zeros_like(self.x_points)
            mask = self.x_points < p
            y[mask] = m * (2 * p * self.x_points[mask] - self.x_points[mask]**2) / p**2
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
        Compute the camber line slope (dy/dx) using finite differences.
        """
        x = self.x_points
        dy = np.zeros_like(x)
        dy[1:-1] = (y_points[2:] - y_points[:-2]) / (x[2:] - x[:-2])
        dy[0] = (y_points[1] - y_points[0]) / (x[1] - x[0])
        dy[-1] = (y_points[-1] - y_points[-2]) / (x[-1] - x[-2])
        return dy

    def compute_cl(self, alpha_deg, y_points):
        """
        Compute the lift coefficient (Cl) using thin-airfoil theory:
          Cl = 2π (α - α0)
        where α0 is determined by integrating the camber slope.
        """
        alpha = np.radians(alpha_deg)
        dy_dx = self.compute_camber_slope(y_points)
        theta = np.arccos(1 - 2 * self.x_points)  # using: x = (1-cosθ)/2
        integrand = dy_dx * (1 - np.cos(theta))
        alpha_0 = -simpson(integrand, x=theta) / np.pi
        Cl = 2 * np.pi * (alpha - alpha_0)
        return Cl

    def compute_cl_vs_alpha(self, alpha_range, y_points, cfd_data=None, do_plot=True):
        """
        Compute and optionally plot Cl vs. angle of attack.
        The do_plot flag (default True) determines whether a new plot is created.
        """
        cl_values = np.array([self.compute_cl(alpha, y_points) for alpha in alpha_range])
        if do_plot:
            plt.figure(figsize=(10, 6))
            plt.plot(alpha_range, cl_values, 'b-', linewidth=2, label="Thin Airfoil Theory")
            if cfd_data is not None:
                plt.plot(cfd_data['alpha'], cfd_data['cl'], 'ro--', linewidth=2, label="ANSYS CFD")
            plt.title("Lift Coefficient (Cl) vs. Angle of Attack")
            plt.xlabel("Angle of Attack (°)")
            plt.ylabel("Cl")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.show()
        return cl_values

    def calculate_cm(self, alpha_deg, y_points):
        """
        Compute the moment coefficient about the aerodynamic center (≈ quarter chord)
        using a Fourier-series expansion of the camber line.
        """
        N = 10
        x_pts = self.x_points
        theta = np.linspace(0, np.pi, 200)
        x_reparam = (1 - np.cos(theta)) / 2
        z = np.interp(x_reparam, x_pts, y_points)
        A = np.zeros(N)
        for n in range(1, N+1):
            integrand = z * np.sin(n * theta)
            A[n-1] = (2/np.pi) * simpson(integrand, x=theta)
        Cm_ac = -np.pi * np.sum([(n * A[n-1]) / (n+1) for n in range(1, N+1)])
        return Cm_ac

    def calculate_velocity_field(self, y_points, alpha_deg=3):
        """
        Compute the velocity field around the airfoil using a simplified vortex panel method.
        Domain: x from -1 to 3, y from -1.5 to 1.5.
        Uses a reduced grid (25×25) for faster computation.
        """
        effective_V = -abs(self.V_inf)
        alpha = np.radians(alpha_deg)
        x_grid = np.linspace(-1, 3, 25)
        y_grid = np.linspace(-1.5, 1.5, 25)
        X, Y = np.meshgrid(x_grid, y_grid)
        u = np.ones_like(X) * effective_V * np.cos(alpha)
        v = np.ones_like(X) * effective_V * np.sin(alpha)
        slopes = self.compute_camber_slope(y_points)
        gamma = 2 * effective_V * (alpha + np.arctan(slopes))
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
                        u_ind = gamma[i] * dy / (2 * np.pi * r_sq)
                        v_ind = -gamma[i] * dx / (2 * np.pi * r_sq)
                        u[j, k] += u_ind * panel_length
                        v[j, k] += v_ind * panel_length
        # Reverse the sign for display.
        u_plot = -u
        v_plot = -v
        return u_plot, v_plot, X, Y

    def calculate_circulation_distribution(self, y_points, alpha_deg):
        """
        Calculate the circulation distribution along the camber line.
        Returns:
          gamma: Array of circulation values at each x-point.
        """
        alpha = np.radians(alpha_deg)
        slopes = self.compute_camber_slope(y_points)
        V_inf = self.V_inf
        c = self.chord
        gamma = np.zeros_like(self.x_points)
        for i, x in enumerate(self.x_points):
            term1 = alpha * np.sqrt(1 - x)
            mask = self.x_points >= x
            if np.any(mask):
                x_subset = self.x_points[mask]
                slopes_subset = slopes[mask]
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

    def calculate_circulation_line_integral(self, u, v, X, Y):
        """
        Calculate circulation using a velocity line integral around a closed circular contour.
        Returns:
          circulation: Circulation computed by integrating u dx + v dy.
        """
        R = 1.5  # radius of the circular contour
        theta = np.linspace(0, 2*np.pi, 1000)
        x_contour = R * np.cos(theta)
        y_contour = R * np.sin(theta)
        dx = np.diff(np.append(x_contour, x_contour[0]))
        dy = np.diff(np.append(y_contour, y_contour[0]))
        u_interp = NearestNDInterpolator(np.column_stack((X.flatten(), Y.flatten())), u.flatten())
        v_interp = NearestNDInterpolator(np.column_stack((X.flatten(), Y.flatten())), v.flatten())
        u_contour = u_interp(x_contour, y_contour)
        v_contour = v_interp(x_contour, y_contour)
        circulation = np.sum(u_contour * dx + v_contour * dy)
        return circulation

    def calculate_bound_circulation(self, y_points, alpha_deg):
        """
        Calculate bound circulation by integrating the circulation distribution along the camber line.
        Returns:
          bound_circ: Total bound circulation.
        """
        gamma = self.calculate_circulation_distribution(y_points, alpha_deg)
        bound_circ = simpson(gamma, x=self.x_points)
        return bound_circ

    def plot_vector_field(self, y_points, alpha_deg=3):
        """
        Plot velocity vectors and streamlines.
        Returns:
          u, v, X, Y: Velocity field (for plotting) and grid arrays.
        """
        u, v, X, Y = self.calculate_velocity_field(y_points, alpha_deg)
        speed = np.sqrt(u**2 + v**2)
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

    def main(self):
        """
        Command-line interactive main function.
        """
        print("Thin Airfoil Theory Analysis")
        print("---------------------------")
        mode = input("Enter 1 for single airfoil analysis or 3 to compare three airfoils: ").strip()
        if mode == "3":
            airfoil_params = []
            for i in range(1, 4):
                print(f"\nDefine Airfoil {i}:")
                atype = input("Enter airfoil type (NACA or custom): ").lower().strip()
                if atype == "naca":
                    naca_code = input("Enter NACA 4-digit code: ").strip()
                    if len(naca_code) != 4 or not naca_code.isdigit():
                        print("Invalid NACA code!")
                        return
                    m_val = int(naca_code[0]) / 100
                    p_val = int(naca_code[1]) / 10
                    y_pts = self.generate_camber_line('naca', m=m_val, p=p_val)
                    label = f"NACA {naca_code}"
                elif atype == "custom":
                    print("Define the camber line as y = f(x):")
                    print("1. Parabolic camber")
                    print("2. Circular arc camber")
                    print("3. S-shaped camber")
                    print("4. Custom polynomial")
                    try:
                        choice = int(input("Enter your choice (1-4): "))
                    except:
                        print("Invalid choice!")
                        return
                    if choice == 1:
                        max_camber = float(input("Enter maximum camber (fraction of chord): "))
                        pos_max_camber = float(input("Enter position of maximum camber (fraction of chord): "))
                        def camber_function(x): 
                            y = np.zeros_like(x)
                            for j, xi in enumerate(x):
                                if xi < pos_max_camber:
                                    y[j] = max_camber * (2 * pos_max_camber * xi - xi**2) / pos_max_camber**2
                                else:
                                    y[j] = max_camber * (1 - 2 * pos_max_camber + 2 * pos_max_camber * xi - xi**2) / (1 - pos_max_camber)**2
                            return y
                        y_pts = camber_function(self.x_points)
                        label = f"Custom Parabolic {i}"
                    elif choice == 2:
                        max_camber = float(input("Enter maximum camber (fraction of chord): "))
                        def camber_function(x):
                            return max_camber * np.sin(np.pi * x)
                        y_pts = camber_function(self.x_points)
                        label = f"Custom Circular {i}"
                    elif choice == 3:
                        max_camber = float(input("Enter maximum camber (fraction of chord): "))
                        def camber_function(x):
                            return max_camber * np.sin(2 * np.pi * x)
                        y_pts = camber_function(self.x_points)
                        label = f"Custom S-shaped {i}"
                    elif choice == 4:
                        coefs_str = input("Enter polynomial coefficients (space-separated, highest order first): ")
                        coefs = [float(c) for c in coefs_str.split()]
                        def camber_function(x):
                            y = np.zeros_like(x)
                            for j, c in enumerate(coefs):
                                y += c * x**(len(coefs) - j - 1)
                            return y
                        y_pts = camber_function(self.x_points)
                        label = f"Custom Poly {i}"
                    else:
                        print("Invalid choice!")
                        return
                else:
                    print("Invalid airfoil type.")
                    return
                airfoil_params.append((y_pts, label))
            
            # Overlay camber line plots.
            plt.figure(figsize=(8, 4))
            for (y_pts, label) in airfoil_params:
                plt.plot(self.x_points, y_pts, label=label, linewidth=2)
            plt.title("Comparison of Airfoil Camber Lines")
            plt.xlabel("x/c")
            plt.ylabel("y/c")
            plt.axis("equal")
            plt.grid(True)
            plt.legend()
            plt.show()
            
            # Overlay Cl vs. angle of attack plots on one figure.
            alpha_range = np.linspace(-3, 12, 16)
            plt.figure(figsize=(10, 6))
            for (y_pts, label) in airfoil_params:
                cl_vals = self.compute_cl_vs_alpha(alpha_range, y_pts, do_plot=False)
                plt.plot(alpha_range, cl_vals, label=label, linewidth=2)
            plt.title("Comparison of Cl vs. Angle of Attack")
            plt.xlabel("Angle of Attack (°)")
            plt.ylabel("Cl")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.show()
            
            # Plot vector field for each airfoil separately at α = 3°.
            for (y_pts, label) in airfoil_params:
                print(f"Plotting velocity field for {label} at α = 3°")
                self.plot_vector_field(y_pts, alpha_deg=3)
            
            # Compute and print circulations for each airfoil.
            for (y_pts, label) in airfoil_params:
                u, v, X, Y = self.calculate_velocity_field(y_pts, 3)
                circ_line = self.calculate_circulation_line_integral(u, v, X, Y)
                gamma = self.calculate_circulation_distribution(y_pts, 3)
                bound_circ = simpson(gamma, x=self.x_points)
                print(f"{label}  Circulation (line integral) = {-1*circ_line:.4f}, Bound circulation = {bound_circ:.4f}")
        elif mode == "1":
            # Single airfoil analysis.
            airfoil_type = input("Enter airfoil type (NACA or custom): ").lower().strip()
            if airfoil_type == "naca":
                naca_code = input("Enter NACA 4-digit code: ").strip()
                if len(naca_code) != 4 or not naca_code.isdigit():
                    print("Invalid NACA code!")
                    return
                m = int(naca_code[0]) / 100
                p = int(naca_code[1]) / 10
                y_points = self.generate_camber_line('naca', m=m, p=p)
            elif airfoil_type == "custom":
                print("Define the camber line as y = f(x):")
                print("1. Parabolic camber")
                print("2. Circular arc camber")
                print("3. S-shaped camber")
                print("4. Custom polynomial")
                try:
                    choice = int(input("Enter your choice (1-4): "))
                except:
                    print("Invalid choice!")
                    return
                if choice == 1:
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
                    y_points = camber_function(self.x_points)
                elif choice == 2:
                    max_camber = float(input("Enter maximum camber (fraction of chord): "))
                    def camber_function(x):
                        return max_camber * np.sin(np.pi * x)
                    y_points = camber_function(self.x_points)
                elif choice == 3:
                    max_camber = float(input("Enter maximum camber (fraction of chord): "))
                    def camber_function(x):
                        return max_camber * np.sin(2 * np.pi * x)
                    y_points = camber_function(self.x_points)
                elif choice == 4:
                    coefs_str = input("Enter polynomial coefficients (space-separated, highest order first): ")
                    coefs = [float(c) for c in coefs_str.split()]
                    def camber_function(x):
                        y = np.zeros_like(x)
                        for i, c in enumerate(coefs):
                            y += c * x**(len(coefs) - i - 1)
                        return y
                    y_points = camber_function(self.x_points)
                else:
                    print("Invalid choice!")
                    return
            else:
                print("Invalid airfoil type.")
                return

            # Plot the camber line.
            plt.figure(figsize=(8, 4))
            plt.plot(self.x_points, y_points, 'k-', linewidth=2)
            plt.title("Airfoil Camber Line")
            plt.xlabel("x/c")
            plt.ylabel("y/c")
            plt.axis('equal')
            plt.grid(True)
            plt.show()

            # Plot the camber slope.
            slopes = self.compute_camber_slope(y_points)
            plt.figure(figsize=(8, 4))
            plt.plot(self.x_points, slopes, 'r-', linewidth=2)
            plt.title("Camber Line Slope (dy/dx)")
            plt.xlabel("x/c")
            plt.ylabel("dy/dx")
            plt.axis('equal')
            plt.grid(True)
            plt.show()

            # Plot Cl vs. angle of attack with optional CFD data overlay.
            alpha_range = np.linspace(-3, 12, 16)
            cfd_choice = input("Do you want to overlay CFD data? (y/n): ").strip().lower()
            cfd_data = None
            if cfd_choice == 'y':
                cfd_alpha = input("Enter CFD alpha values (comma-separated): ").split(",")
                cfd_cl = input("Enter corresponding CFD Cl values (comma-separated): ").split(",")
                cfd_data = {'alpha': np.array([float(val) for val in cfd_alpha]),
                            'cl': np.array([float(val) for val in cfd_cl])}
            cl_values = self.compute_cl_vs_alpha(alpha_range, y_points, cfd_data=cfd_data)
            
            # Compute circulations using the velocity field at α = 3°.
            alpha_test = 3
            u, v, X, Y = self.calculate_velocity_field(y_points, alpha_test)
            circ_line = self.calculate_circulation_line_integral(u, v, X, Y)
            print(f"Circulation (line integral) = {-1*circ_line:.4f}")
            gamma = self.calculate_circulation_distribution(y_points, alpha_test)
            bound_circ = simpson(gamma, x=self.x_points)
            print(f"Bound circulation (integrated) = {bound_circ:.4f}")

            # Plot the velocity field at α = 3°.
            self.plot_vector_field(y_points, alpha_deg=3)
        else:
            print("Invalid mode selected!")

# =============================================================================
# Main Program Entry Point (CLI Version)
# =============================================================================
if __name__ == "__main__":
    analyzer = ThinAirfoilAnalyzer()
    analyzer.main()
