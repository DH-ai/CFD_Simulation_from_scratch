import csv
import regex as re
import numpy as np

import matplotlib.pyplot as plt

CSV_FILE = './user-000.csv'
class Airfoil:
    def __init__(self):
        self.airfoil = np.array([[]]).reshape(-1,2)
        self.name = None
        self.maximum_camber = None # M in percentage
        self.maximum_camber_position = None # P in percentage
        self.thickness = None # T in percentage
        self.chord_length = None # in mm default
        self.x = None


    def  csv_air_foil(self, csv_file):
        with open(csv_file, 'r') as file:
            
            reader = csv.reader(file)
            iter_reader = iter(reader)
            for row in iter_reader:
                # print(row) 
                
                self.airfoil = np.append(self.airfoil,[[float(row[0]),float(row[1])]],axis=0)   
            # self.airfoil.reshape(-1,2)
            np.delete(self.airfoil, 0, axis=0) 





af = Airfoil()
af.csv_air_foil(CSV_FILE)
# print(af.airfoil)

# plt.plot(af.airfoil[:,0],af.airfoil[:,1])
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()


x = np.linspace(-5, 5, 500)  # X values
y = np.linspace(-5, 5, 500)  # Y values
X ,Y= np.meshgrid(x, y)
Z = -((X**2 + Y**2) / 10)
# Create a 3D figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Create 3D contour plot
ax.contour3D(X, Y, Z, levels=300, cmap='viridis')

# Labels and title
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.set_title("3D Contour Plot")

plt.show()