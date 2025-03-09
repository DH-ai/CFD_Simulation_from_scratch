# Thin Airfoil Theory Implementation

This package implements thin airfoil theory for analyzing airfoil performance and contains tools for comparing results with CFD simulations.

## Project Files

1. `thin_airfoil_theory.py` - Main implementation of thin airfoil theory
2. `naca_analysis.py` - Script to analyze NACA airfoils and compare with CFD
3. `custom_analysis.py` - Script to analyze three custom airfoil designs
4. `requirements.txt` - Required dependencies
5. `thinairfoilanalyser.py` – GUI program for interactive thin airfoil analysis.

## Installation

1. Make sure you have Python 3.8 or newer installed
2. Create and activate a virtual environment (recommended):
   ```
   python -m venv venv
   
   # On Windows:
   venv\Scripts\activate
   
   # On macOS/Linux:
   source venv/bin/activate
   ```
3. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

### Running the Main Program

To run the interactive program:

```
python thin_airfoil_theory.py
```

This will ask you to input airfoil specifications and analyze the airfoil.

### Running the GUI
To run the GUI

```
python thinairfoilanalyser.py
```

This will ask for parameters and live update the simulations.

### NACA Airfoil Analysis

To analyze a NACA airfoil and compare with CFD results:

```
python naca_analysis.py
```

Note: You may need to modify the CFD data in `naca_analysis.py` to match your actual CFD results from Assignment 1.

### Custom Airfoil Analysis

To analyze the three custom airfoil designs:

```
python custom_analysis.py
```

## Features

- Camber line generation for NACA and custom airfoils
- Calculation of camber line slope
- Lift coefficient prediction based on thin airfoil theory
- Vector field plotting around airfoils
- Circulation calculation using line integral and bound vorticity approaches
- Comparative analysis of multiple airfoil designs

## Modifying Parameters

- To modify the NACA airfoil parameters, edit the `m` and `p` values in `naca_analysis.py`
- To create your own custom airfoils, edit the custom airfoil functions in `custom_analysis.py`
- To change the CFD comparison data, replace the placeholder data in `naca_analysis.py`

## Output Files

The analysis scripts will generate the following output files:

- `naca_camber_line.png` - Plot of NACA airfoil camber line
- `naca_camber_slope.png` - Plot of NACA airfoil camber line slope
- `naca_cl_vs_alpha.png` - Comparison of Cl vs alpha curves
- `naca_vector_field.png` - Vector field around NACA airfoil
- `custom_camber_comparison.png` - Comparison of all airfoil camber lines
- `custom_cl_comparison.png` - Comparison of all Cl vs alpha curves
- `custom1_vector_field.png` - Vector field around custom airfoil 1
- `custom2_vector_field.png` - Vector field around custom airfoil 2
- `custom3_vector_field.png` - Vector field around custom airfoil 3
