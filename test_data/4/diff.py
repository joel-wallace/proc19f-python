import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.colors as mcolors

# Set global font family, font size, and linewidth using rcParams
line_width = 0.5  # Set a smaller line width for the plot
spine_width = 1.0

plt.rcParams['font.family'] = 'Arial'  # Keep font as Arial
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = line_width  # Global line width for axes
plt.rcParams['lines.linewidth'] = line_width  # Global line width for plot lines

# Load diffusion spectrum data
diff_data = np.loadtxt('diff_data.csv', delimiter=',')
x_diff = diff_data[:, 0]  # Chemical shift (¹H ppm)
y_diff = diff_data[:, 1:]  # Intensity values

# Load error bar data (G² vs ln(I/I₀))
errorbar_data = np.loadtxt('errorbar_data.csv', delimiter=',')
x_err = errorbar_data[:, 0]  # G² values
y_err = errorbar_data[:, 1]  # ln(I/I₀) values
y_err_values = errorbar_data[:, 2]  # Error in ln(I/I₀)

# Perform linear regression (Best-fit line)
slope, intercept, r_value, p_value, std_err = linregress(x_err, y_err)

# Generate best-fit line values
x_fit = np.linspace(min(x_err), max(x_err), 100)  # Smooth range for plotting
y_fit = slope * x_fit + intercept

# Create figure with two subplots (4 inch by 2 inch)
fig, axs = plt.subplots(1, 2, figsize=(4, 2))

# Create a custom color map for the spectrum plot (saturation gradient of #889E73)
color_map = mcolors.LinearSegmentedColormap.from_list(
    "green_gradient", ["#cde8b3", "#889E73"], N=y_diff.shape[1]
)

# Plot diffusion spectrum (Chemical shift vs Intensity) in the first subplot
for i in range(y_diff.shape[1]):  # Iterate over all columns in y_diff (now includes all columns)
    axs[0].plot(x_diff, y_diff[:, i], color=color_map(i / y_diff.shape[1]))  # Green gradient
axs[0].invert_xaxis()  # Reverse x-axis
axs[0].set_xlim([-63, -62])  # Zoom in to -62 to -63 ppm
axs[0].set_xlabel(r'$\delta_{\mathrm{F}}$ (ppm)')  # Use \mathrm to ensure F is non-italic

# Remove all spines except the bottom one for the first subplot
axs[0].spines['left'].set_color('none')
axs[0].spines['right'].set_color('none')  # Hide the right spine
axs[0].spines['top'].set_color('none')  # Hide the top spine
axs[0].spines['bottom'].set_linewidth(spine_width)  # Set bottom spine width to match
axs[0].tick_params(axis='y', which='both', length=0)  # Remove y-axis ticks
axs[0].set_ylabel('')  # Remove y-axis label
axs[0].set_yticks([])
axs[0].invert_xaxis()
# Plot data with error bars and best-fit line in the second subplot
# Removed the title from the second subplot
# Create a color map for the error bar points (same gradient as the spectrum plot)
error_point_colors = [color_map(i / y_diff.shape[1]) for i in range(len(y_err))]

# Plot the data points with the gradient color and reduced dot size
for i in range(len(y_err)):
    axs[1].errorbar(x_err[i], y_err[i], yerr=y_err_values[i], fmt='o', color=error_point_colors[i], capsize=3, markersize=4)  # Reduced dot size to 4

# Plot the best-fit line with the specified color (FF7F3E - orange)
axs[1].plot(x_fit, y_fit, color='#FF7F3E', label=f"Fit: y = {slope:.4f}x + {intercept:.4f}")  # Best-fit line in orange

# Use LaTeX formatting for labels
axs[1].set_xlabel(r'G$^{2}$ (T$^{2}$ m$^{-2}$)')  # Subscript and superscript formatting for G²
axs[1].set_ylabel(r'ln(I/I$_{0}$)')  # Subscript for I₀ using LaTeX formatting

# Remove all spines except the bottom one for the second subplot
axs[1].spines['left'].set_linewidth(spine_width)
axs[1].spines['right'].set_color('none')  # Hide the right spine
axs[1].spines['top'].set_color('none')  # Hide the top spine
axs[1].spines['bottom'].set_linewidth(spine_width)  # Set bottom spine width to match

# Save the figure as a PNG file
plt.tight_layout()
plt.savefig('diffusion_spectrum_and_fit.png', format='png', dpi=300)

# Show the figure with both plots
plt.show()

# Calculate diffusion coefficient (D) using MATLAB's prefactor formula
gamma = 2.68E+08  # Gyromagnetic ratio (rad/s/T)
delta = 0.004  # Gradient pulse duration (s)
DELTA = 0.1  # Diffusion time (s)
sigma = 0.9  # Shape factor (depends on gradient shape)

prefactor = gamma**2 * delta**2 * sigma**2 * (DELTA - delta/3)

D = -slope / prefactor  # Diffusion coefficient (m²/s)
D_err = std_err / prefactor  # Error in D

print(f"Diffusion Coefficient (D): {D:.4e} m²/s ± {D_err:.4e}")

# Estimate hydrodynamic radius (rH) using Stokes-Einstein equation
T = 298  # Temperature (K)
n = 0.89E-3  # Viscosity of water (Pa·s)
kb = 1.3806488e-23  # Boltzmann constant (J/K)

rH = kb * T / (6 * np.pi * n * D)  # Hydrodynamic radius (m)
rH_err = rH * (D_err / D)

print(f"Hydrodynamic Radius (rH): {rH*1E9:.2f} nm ± {rH_err*1E9:.2f} nm")

# Estimate molecular weight (MW) assuming a folded protein
MW = rH * 1E9 * 10  # Convert to kDa
MW_err = rH_err * 1E9 * 10

print(f"Estimated Molecular Weight: {MW:.2f} kDa ± {MW_err:.2f} kDa")
