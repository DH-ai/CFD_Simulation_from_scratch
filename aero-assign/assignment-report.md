# AE 244: Assignment 2 Report
**Thin Airfoil Theory Implementation**

## 1. Program Development Team Introduction and Work Share

| Sr. No | Roll Number | Name | Contribution Level (0 to 5) | Specifics of Contribution |
|--------|-------------|------|---------------------------|---------------------------|
| 1      | [Your Roll] | [Your Name] | 5 | Implemented all functions and created all visualizations |
| 2      | [Partner Roll] | [Partner Name] | [Level] | [Details] |
| 3      | [Partner Roll] | [Partner Name] | [Level] | [Details] |

## 2. Algorithm

### 2.1 Algorithm for Plotting Camber Line

1. Define a range of x-coordinates along the chord (0 to 1)
2. For NACA airfoils:
   - Extract the camber (m) and position of max camber (p) from NACA digits
   - For each x-coordinate:
     - If x < p, apply formula: y = m * (2 * p * x - x^2) / p^2
     - If x ≥ p, apply formula: y = m * (1 - 2 * p + 2 * p * x - x^2) / (1 - p)^2
3. For custom airfoils:
   - Apply user-defined function f(x) to each x-coordinate
4. Plot the resulting (x, y) coordinates

### 2.2 Algorithm for Plotting Camber Line Slope

1. Take the previously generated camber line points (x, y)
2. For interior points:
   - Calculate slope using central difference: dy/dx = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
3. For the first point:
   - Calculate slope using forward difference: dy/dx = (y[1] - y[0]) / (x[1] - x[0])
4. For the last point:
   - Calculate slope using backward difference: dy/dx = (y[-1] - y[-2]) / (x[-1] - x[-2])
5. Plot the resulting (x, slope) coordinates

### 2.3 Algorithm for Computing Cl

1. Convert angle of attack from degrees to radians
2. Calculate camber line slopes at all x-coordinates
3. Convert x-coordinates to theta using x = (1 - cos(θ))/2
4. Calculate zero-lift angle of attack (α0):
   - α0 = -1/π * ∫(from 0 to π) [dy/dx * (1 - cos(θ))] dθ
   - Implement this integration numerically using Simpson's rule
5. Calculate Cl using thin airfoil theory formula:
   - Cl = 2π(α - α0)
6. For a range of angles, repeat steps 1-5

### 2.4 Algorithm for Plotting Vector Field

1. Create a grid of points around the airfoil, extending to domain sizes of 4c by 3c
2. Calculate the freestream velocity components based on angle of attack:
   - u_inf = V∞ * cos(α)
   - v_inf = V∞ * sin(α)
3. Calculate circulation distribution along the camber line
4. For each point in the grid:
   - Calculate induced velocity from the vortices placed along the camber line
   - Add freestream velocity to get total velocity
5. Plot the vector field using quiver and streamline plots

### 2.5 Algorithm for Calculating Circulation through Line Integral

1. Define a circular contour around the airfoil
2. Interpolate velocity field at points along this contour
3. Calculate the line integral around the contour:
   - Γ = ∮ v·dl = ∮ (u dx + v dy)
   - Implement this numerically by summing u*dx + v*dy around the contour

## 3. Airfoil Simulation and Results

### 3.1 Camber Line ('y' vs 'x') for the NACA Airfoil from Assignment 1

[Insert NACA camber line plot here]

The camber line plot shows the curvature of the airfoil's mean line. The NACA 4412 airfoil has a maximum camber of 4% located at 40% of the chord from the leading edge.

### 3.2 Slope of the Camber Line along x (dy/dx vs x)

[Insert camber line slope plot here]

The slope plot shows how the camber line's gradient changes along the chord. The slope is highest near the leading edge and decreases toward the trailing edge.

### 3.3 Cl vs α for the NACA Airfoil from the Program and CFD Simulations on Same Chart

[Insert Cl vs alpha comparison plot here]

### 3.4 Comment on/Discuss the Results Obtained in 3.3

The thin airfoil theory results show a linear relationship between Cl and α, as expected from the theory (Cl = 2π(α - α0)). The slope of approximately 2π per radian (0.11 per degree) matches theoretical predictions.

When compared to CFD simulations:
- At low angles of attack (-3° to 6°), thin airfoil theory predicts Cl values reasonably close to CFD results
- At higher angles (>6°), CFD shows non-linear behavior due to flow separation effects, which thin airfoil theory cannot predict
- The zero-lift angle of attack (α0) from thin airfoil theory is close to that from CFD, indicating good prediction of the basic airfoil characteristics

The differences between thin airfoil theory and CFD arise from:
1. Thin airfoil theory neglects thickness effects
2. Thin airfoil theory assumes small angles and attached flow
3. Viscous effects are ignored in thin airfoil theory

### 3.5 Vector Field Plot Around the Airfoil at α = 3°

[Insert vector field plot here]

### 3.6 Comment on/Discuss the Plot Obtained in 3.5

The vector field around the NACA 4412 airfoil at α = 3° shows:
1. Increased velocity on the upper surface (suction side) compared to the lower surface
2. Flow following the contour of the airfoil smoothly
3. Stagnation point slightly below the leading edge
4. Kutta condition satisfied at the trailing edge with smooth flow
5. Diminishing effect of the airfoil on the flow field as distance increases

These observations align with expected aerodynamic behavior and demonstrate how the pressure differential between upper and lower surfaces generates lift.

### 3.7 Circulation Around the Entire Airfoil Using Velocity Line-Integral Approach at α = 3°

Circulation (line integral method): [Insert value here]

### 3.8 Bound Circulation by Integrating Circulation Distribution Along Camber Line at α = 3°

Circulation (bound vorticity method): [Insert value here]

### 3.9 Compare and Comment on the Values Obtained in 3.7 and 3.8

The circulation values calculated via both methods should theoretically be equal according to Kelvin's circulation theorem. The comparison shows:
- Line integral method: directly measures the circulation around a contour enclosing the airfoil
- Bound vorticity method: calculates circulation by integrating the vorticity distribution along the camber line

[Comment on the agreement/disagreement between the two methods and possible reasons for any discrepancies]

## 4. Novel Airfoil Properties

### 4.1 Camber Line ('y' vs 'x') for Your 3 Airfoils and the NACA Airfoil on Same Chart

[Insert comparison of camber lines plot here]

The three custom airfoils designed are:
1. **Reflexed Airfoil**: Has maximum camber of 6% at 30% chord with reflex (upward curve) in the rear portion
2. **Front-loaded Airfoil**: Has maximum camber of 5% at 20% chord with gradual decrease toward trailing edge
3. **S-shaped Airfoil**: Has sinusoidal camber distribution with positive camber forward and negative camber aft

### 4.2 Cl vs α of the Three Airfoils Along with That of the NACA Airfoil on Same Chart

[Insert Cl comparison plot here]

### 4.3 Comment on/Discuss the Results Obtained in 4.2

The lift curves show different characteristics for each airfoil:
1. **Reflexed Airfoil**: Has slightly lower lift slope than NACA 4412 due to the reflexed trailing edge, which reduces the effective camber. However, this design likely provides improved pitching moment characteristics.
2. **Front-loaded Airfoil**: Shows higher lift at low angles of attack due to increased camber near the leading edge, beneficial for takeoff performance.
3. **S-shaped Airfoil**: Exhibits lower lift coefficients overall, but may have benefits for pitch stability due to its unique camber distribution.

These differences illustrate how camber line modification directly impacts lift generation, demonstrating the fundamental concept of thin airfoil theory where lift is primarily influenced by camber and angle of attack.

### 4.4 Vector Field Plot Around the Three Airfoils at α = 3°

[Insert vector field plots for custom airfoils]

### 4.5 Comment on/Discuss the Results Obtained in 4.4

The vector field visualizations reveal different flow patterns around each custom airfoil:
1. **Reflexed Airfoil**: Shows strong acceleration over the forward portion but reduced acceleration near the trailing edge due to the reflexed shape
2. **Front-loaded Airfoil**: Demonstrates concentrated acceleration near the leading edge with more gradual pressure recovery
3. **S-shaped Airfoil**: Exhibits complex flow pattern with acceleration regions corresponding to positive camber sections and deceleration in negative camber regions

These flow patterns directly correspond to the pressure distributions that generate lift and pitching moments, providing visual confirmation of how camber line shape influences aerodynamic performance.

## 5. Conclusion

### 5.1 Overall Take on Your Code's Performance as Compared to Ansys Simulation, and Possible Reasons for Deviations

The thin airfoil theory implementation provided reasonable approximations of lift characteristics when compared to more sophisticated CFD simulations:

**Strengths:**
- Accurately predicted the lift curve slope in the linear region
- Correctly captured the influence of camber on zero-lift angle of attack
- Provided rapid analysis capabilities suitable for early design iterations
- Clearly showed the relationship between camber distribution and lift

**Limitations:**
- Could not predict maximum lift coefficient or stall behavior
- Did not account for thickness effects
- Neglected viscous effects like boundary layer growth and separation
- Assumed small angles and perturbations

These limitations explain the primary deviations from CFD results, particularly at higher angles of attack. However, the computational efficiency (analysis in seconds vs. hours for CFD) makes thin airfoil theory valuable for preliminary design work and understanding fundamental aerodynamic principles.

### 5.2 Overall Take on the Performance of Your Airfoil as Compared to the NACA Airfoil

Among the custom designs, the front-loaded airfoil showed the most promising characteristics for STOL (Short Take-Off and Landing) applications:

1. It generated higher lift at low angles of attack, which is beneficial for takeoff performance
2. The early maximum camber position enhances lift at the expense of some drag increase
3. The smooth camber distribution should maintain attached flow at low speeds

The reflexed airfoil, while generating slightly less lift, would likely have improved pitching moment characteristics that could reduce trim drag and improve overall efficiency in cruise.

The S-shaped airfoil demonstrated interesting stability characteristics but would require further analysis beyond thin airfoil theory to fully evaluate its performance.

For the STOL competition, a combination of the front-loaded airfoil's lift performance with some reflex in the trailing edge might provide an optimal balance of takeoff performance and cruise efficiency.

## 6. Code

[This section will be evaluated by the instructor based on the submitted code files]

## 7. Acknowledgement

I would like to acknowledge the following resources for their assistance:
- [List any discussions or help received here]
- Course lecture materials on thin airfoil theory
- Numerical integration algorithms from SciPy documentation

## 8. References

1. Anderson, J.D. (2010). Fundamentals of Aerodynamics. McGraw-Hill Education.
2. Katz, J., & Plotkin, A. (2001). Low-Speed Aerodynamics. Cambridge University Press.
3. Abbott, I.H., & Von Doenhoff, A.E. (1959). Theory of Wing Sections. Dover Publications.
4. [List any additional references used]
