import numpy as np
import matplotlib.pyplot as plt

# Define a grid
x = np.linspace(-5, 5, 10)
y = np.linspace(-5, 5, 10)
X, Y = np.meshgrid(x, y)

# Define vector components (e.g., velocity components)
U = -Y  # Horizontal component
V = X   # Vertical component

# Create quiver plot
plt.figure(figsize=(6,6))
# plt.quiver(X, Y, U, V, scale=20, color='b')

# Labels and title
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Quiver Plot Example")
plt.grid()
plt.show()
