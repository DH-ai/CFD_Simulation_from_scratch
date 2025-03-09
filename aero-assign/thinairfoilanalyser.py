import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Use the TkAgg backend for Tkinter
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from thin_airfoil_theory import ThinAirfoilAnalyzer

# ----------------------------------------------------------------
#  GUI Class
#  This class creates a GUI with inputs for airfoil type, velocity, and angle of attack.
#  It uses both text entries and sliders for live updates.
# ----------------------------------------------------------------
class AirfoilGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Thin Airfoil GUI (Bonus)")
        self.geometry("1050x700")

        self.analyzer = ThinAirfoilAnalyzer()

        # Define GUI variables
        self.airfoil_type_var = tk.StringVar(value="naca")
        self.naca_code_var = tk.StringVar(value="4412")
        self.custom_poly_var = tk.StringVar(value="0.05, 0, 0")
        self.velocity_var = tk.DoubleVar(value=1.0)
        self.alpha_var = tk.DoubleVar(value=3.0)
        self.tabulated_points_var = tk.StringVar(value="0.0 0.0\n0.5 0.02\n1.0 0.0")

        self.build_ui()

        # Create a matplotlib figure embedded in the GUI
        self.fig = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Update plot on start
        self.update_plot()

    def build_ui(self):
        """
        Build the control panel on the left side with inputs and sliders.
        """
        control_frame = ttk.Frame(self, padding=5)
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        # Airfoil type radio buttons
        ttk.Label(control_frame, text="Airfoil Type:").pack(pady=2)
        ttk.Radiobutton(control_frame, text="NACA 4-digit", variable=self.airfoil_type_var, value="naca").pack(anchor=tk.W)
        ttk.Radiobutton(control_frame, text="Custom Polynomial", variable=self.airfoil_type_var, value="custom").pack(anchor=tk.W)
        ttk.Radiobutton(control_frame, text="Tabulated Points", variable=self.airfoil_type_var, value="tabulated").pack(anchor=tk.W)

        # NACA code input
        ttk.Label(control_frame, text="NACA code (e.g. 4412):").pack(pady=2)
        naca_entry = ttk.Entry(control_frame, textvariable=self.naca_code_var, width=10)
        naca_entry.pack()

        # Custom polynomial input
        ttk.Label(control_frame, text="Custom Polynomial Coeffs:").pack(pady=2)
        custom_entry = ttk.Entry(control_frame, textvariable=self.custom_poly_var, width=20)
        custom_entry.pack()

        # Tabulated points input
        ttk.Label(control_frame, text="Tabulated Points (x y):").pack(pady=2)
        tabulated_text = ttk.Entry(control_frame, textvariable=self.tabulated_points_var, width=25)
        tabulated_text.pack()

        # Velocity input: text entry and slider
        ttk.Label(control_frame, text="Velocity (m/s):").pack(pady=2)
        vel_entry = ttk.Entry(control_frame, textvariable=self.velocity_var, width=6)
        vel_entry.pack(pady=2)
        ttk.Label(control_frame, text="Velocity Slider:").pack(pady=2)
        vel_slider = ttk.Scale(control_frame, from_=0.5, to=5.0, variable=self.velocity_var, orient=tk.HORIZONTAL, command=self.on_slider_change)
        vel_slider.pack(fill=tk.X, pady=2)

        # Angle of attack input: text entry and slider
        ttk.Label(control_frame, text="Angle of Attack (°):").pack(pady=2)
        alpha_entry = ttk.Entry(control_frame, textvariable=self.alpha_var, width=6)
        alpha_entry.pack(pady=2)
        ttk.Label(control_frame, text="Angle Slider:").pack(pady=2)
        alpha_slider = ttk.Scale(control_frame, from_=-5, to=20, variable=self.alpha_var, orient=tk.HORIZONTAL, command=self.on_slider_change)
        alpha_slider.pack(fill=tk.X, pady=2)

        # Update Plot button
        update_button = ttk.Button(control_frame, text="Update Plot", command=self.update_plot)
        update_button.pack(pady=5)

        # Labels for displaying Cl and Cm
        self.cl_label = ttk.Label(control_frame, text="Cl: 0.00")
        self.cl_label.pack(pady=2)
        self.cm_label = ttk.Label(control_frame, text="Cm: 0.00")
        self.cm_label.pack(pady=2)

        # Frame for embedding the plot
        self.plot_frame = ttk.Frame(self)
        self.plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def on_slider_change(self, event=None):
        """
        Update the plot when a slider changes.
        """
        self.update_plot()

    def parse_custom_polynomial(self):
        """
        Parse the custom polynomial coefficients from text.
        Returns a function f(x) that computes the polynomial.
        """
        coeffs_str = self.custom_poly_var.get().replace(" ", "")
        if not coeffs_str:
            return lambda x: np.zeros_like(x)
        try:
            coeffs = [float(c) for c in coeffs_str.split(",")]
        except:
            messagebox.showerror("Error", "Invalid polynomial coefficients!")
            return lambda x: np.zeros_like(x)
        def poly_func(x):
            y = np.zeros_like(x)
            order = len(coeffs) - 1
            for i, coef in enumerate(coeffs):
                y += coef * x**(order - i)
            return y
        return poly_func

    def parse_tabulated_points(self):
        """
        Parse the tabulated points input (one point per line: "x y").
        Returns arrays for x and y.
        """
        text = self.tabulated_points_var.get().strip()
        lines = text.split("\n")
        tab_x, tab_y = [], []
        for ln in lines:
            parts = ln.strip().split()
            if len(parts) == 2:
                x_val, y_val = float(parts[0]), float(parts[1])
                tab_x.append(x_val)
                tab_y.append(y_val)
        if len(tab_x) < 2:
            messagebox.showerror("Error", "Need at least two points for tabulated input!")
            return None, None
        return np.array(tab_x), np.array(tab_y)

    def update_plot(self):
        """
        Read all user inputs, compute the camber line and flow parameters,
        update Cl and Cm labels, and update the plot.
        """
        # Update freestream velocity and angle-of-attack from the shared variables
        self.analyzer.V_inf = float(self.velocity_var.get())
        alpha_deg = float(self.alpha_var.get())

        # Determine which airfoil type is selected and generate the camber line accordingly
        airfoil_type = self.airfoil_type_var.get()
        if airfoil_type == "naca":
            code = self.naca_code_var.get().strip()
            if len(code) == 4 and code.isdigit():
                m = int(code[0]) / 100
                p = int(code[1]) / 10
                y_points = self.analyzer.generate_camber_line('naca', m=m, p=p)
            else:
                messagebox.showerror("Error", "Invalid NACA code!")
                return
        elif airfoil_type == "custom":
            func = self.parse_custom_polynomial()
            y_points = self.analyzer.generate_camber_line('custom', function=func)
        elif airfoil_type == "tabulated":
            tab_x, tab_y = self.parse_tabulated_points()
            if tab_x is None:
                return
            y_points = self.analyzer.generate_camber_line('tabulated', tab_x=tab_x, tab_y=tab_y)
        else:
            y_points = np.zeros_like(self.analyzer.x_points)

        # Compute lift (Cl) and moment (Cm) coefficients
        Cl = self.analyzer.compute_cl(alpha_deg, y_points)
        Cm = self.analyzer.calculate_cm(alpha_deg, y_points)

        # Update the Cl and Cm labels
        self.cl_label.config(text=f"Cl: {Cl:.3f}")
        self.cm_label.config(text=f"Cm: {Cm:.3f}")

        # Clear the current plot
        self.ax.clear()
        # Plot the camber line (airfoil geometry)
        self.ax.plot(self.analyzer.x_points, y_points, 'k-', linewidth=2, label="Camber Line")
        # Set the axis limits and labels
        self.ax.set_xlim(-1, 3)
        self.ax.set_ylim(-1.5, 1.5)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_title(f"Airfoil Flow at α={alpha_deg}°, V={self.analyzer.V_inf} m/s\nCl={Cl:.3f}, Cm={Cm:.3f}")
        self.ax.set_xlabel("x/c")
        self.ax.set_ylabel("y/c")
        self.ax.legend()

        # Compute and plot the velocity field (vectors and streamlines)
        u, v, X, Y = self.analyzer.calculate_velocity_field(y_points, alpha_deg)
        speed = np.sqrt(u**2 + v**2)
        # Mask a narrow band around the camber line only in the central region (avoid masking near tip/tail)
        tol = 0.0002
        f_interp = interp1d(self.analyzer.x_points, y_points, bounds_error=False, fill_value="extrapolate")
        mask = (np.abs(Y - f_interp(X)) < tol) & (X > 0.1) & (X < 0.9)
        u[mask] = np.nan
        v[mask] = np.nan

        # Plot velocity vectors using quiver
        self.ax.quiver(X, Y, u, v, scale=50, width=0.002, color='gray')
        # Plot streamlines, colored by velocity magnitude
        self.ax.streamplot(X, Y, u, v, density=2, color=speed, cmap='jet', linewidth=1)

        # Redraw the canvas with the new plot
        self.canvas.draw()

# ----------------------------------------------------------------
#  Main Program Entry Point
# ----------------------------------------------------------------
if __name__ == "__main__":
    app = AirfoilGUI()
    app.mainloop()
