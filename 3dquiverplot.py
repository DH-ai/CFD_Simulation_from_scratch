from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
# Define a 3D grid
z = np.linspace(-5, 5, 10)
y= np.linspace(-5, 5, 10)
x = np.linspace(-5, 5, 10)
X, Y = np.meshgrid(x, y)

# Define velocity components
U = -Y
V = X
W =  np.sqrt(U**2 + V**2)  # Some variation in the Z-direction

# Create a 3D quiver plot
plt.figure(figsize=(6,6))
quiver = plt.quiver(X, Y, U, V, W, cmap='viridis')

# Add colorbar
plt.colorbar(quiver, label="Velocity Magnitude")

# Labels and title
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Colored Quiver Plot")
plt.grid()
plt.show(

    
)
