import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt

# Analytical Method Function
def balance_rotating_masses_single_plane(masses, radii, angles):
    """
    Balances rotating masses in a single plane using the analytical method.
    """
    # Convert angles to radians
    angles_rad = np.radians(angles)
    
    # Calculate the x and y components of centrifugal forces
    Fx = np.sum([m * r * np.cos(a) for m, r, a in zip(masses, radii, angles_rad)])
    Fy = np.sum([m * r * np.sin(a) for m, r, a in zip(masses, radii, angles_rad)])
    
    # Resultant force
    F_resultant = np.sqrt(Fx**2 + Fy**2)
    theta_resultant = np.arctan2(Fy, Fx)  # Resultant angle in radians
    
    # Balancing mass and position
    balancing_mass = F_resultant / radii[-1]  # Using the last radius for simplicity
    balancing_angle = np.degrees(theta_resultant) + 180  # Opposite direction
    
    return {
        "balancing_mass": balancing_mass,
        "balancing_angle": balancing_angle,
        "Fx": Fx,
        "Fy": Fy,
    }

# Graph Plotting Function
def plot_balancing(masses, radii, angles, result):
    """
    Plots the graphical representation of balancing rotating masses in a single plane.
    """
    angles_rad = np.radians(angles)
    x_components = [m * r * np.cos(a) for m, r, a in zip(masses, radii, angles_rad)]
    y_components = [m * r * np.sin(a) for m, r, a in zip(masses, radii, angles_rad)]
    
    # Determine scaling factors for arrow size
    max_magnitude = max(np.sqrt(np.array(x_components)**2 + np.array(y_components)**2)) + 1
    arrow_scale = 0.05  # Scale factor for arrows
    head_width = max_magnitude * arrow_scale
    head_length = max_magnitude * arrow_scale * 0.8

    # Plot original vectors
    plt.figure(figsize=(8, 8))
    for i, (x, y) in enumerate(zip(x_components, y_components)):
        plt.arrow(0, 0, x, y, head_width=head_width, head_length=head_length, 
                  length_includes_head=True, label=f"Mass {i+1}")
        # Annotate the angle and notation
        angle_text = f'{angles[i]}°'
        plt.text(x * 1.1, y * 1.1, f'M{i+1} ({angle_text})', fontsize=12, color='blue', ha='center')

    # Plot resultant vector
    plt.arrow(0, 0, result["Fx"], result["Fy"], color="red", head_width=head_width, head_length=head_length,
              length_includes_head=True, label="Resultant")
    plt.text(result["Fx"] * 1.1, result["Fy"] * 1.1, f'R ({result["Fx"]:.3f}, {result["Fy"]:.3f})', fontsize=12, color='red', ha='center')

    # Plot balancing vector
    balancing_angle_rad = np.radians(result["balancing_angle"])
    bx = result["balancing_mass"] * radii[-1] * np.cos(balancing_angle_rad)
    by = result["balancing_mass"] * radii[-1] * np.sin(balancing_angle_rad)
    plt.arrow(0, 0, bx, by, color="green", head_width=head_width, head_length=head_length, 
              length_includes_head=True, label="Balancing Mass")
    plt.text(bx * 1.1, by * 1.1, f'B ({result["balancing_mass"]:.2f} kg)', fontsize=12, color='green', ha='center')

    # Customize plot
    plt.axhline(0, color="black", linewidth=0.5)
    plt.axvline(0, color="black", linewidth=0.5)
    plt.grid(True)
    plt.axis("equal")
    plt.legend()
    plt.title("Balancing of Rotating Masses (Single Plane)")
    plt.xlabel("Fx")
    plt.ylabel("Fy")
    plt.show()

# GUI Code
def get_gui_input():
    """
    Initializes a GUI for user input and calls the calculation and plotting functions.
    """
    def on_calculate():
        try:
            # Get user inputs
            masses = list(map(float, mass_entry.get().strip().split(',')))
            radii = list(map(float, radius_entry.get().strip().split(',')))
            angles = list(map(float, angle_entry.get().strip().split(',')))
            
            # Ensure equal number of inputs
            if not (len(masses) == len(radii) == len(angles)):
                raise ValueError("Lengths of masses, radii, and angles must be the same.")

            # Perform calculations
            result = balance_rotating_masses_single_plane(masses, radii, angles)
            
            # Display results
            results_text.set(
                f"Balancing Mass: {result['balancing_mass']:.3f} kg\n"
                f"Balancing Angle: {result['balancing_angle']:.3f}°\n"
                f"Resultant Force (Fx): {result['Fx']:.3f} N\n"
                f"Resultant Force (Fy): {result['Fy']:.3f} N"
            )
            
            # Plot results
            plot_balancing(masses, radii, angles, result)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # Create main window
    root = tk.Tk()
    root.title("Balancing Rotating Masses (Single Plane)")

    # Input frame
    input_frame = ttk.Frame(root, padding="10")
    input_frame.grid(row=0, column=0, sticky="EW")

    # Input fields
    ttk.Label(input_frame, text="Masses (kg, comma-separated):").grid(row=0, column=0, sticky="W")
    mass_entry = ttk.Entry(input_frame, width=50)
    mass_entry.grid(row=0, column=1, sticky="EW")

    ttk.Label(input_frame, text="Radii (m, comma-separated):").grid(row=1, column=0, sticky="W")
    radius_entry = ttk.Entry(input_frame, width=50)
    radius_entry.grid(row=1, column=1, sticky="EW")

    ttk.Label(input_frame, text="Angles (degrees, comma-separated):").grid(row=2, column=0, sticky="W")
    angle_entry = ttk.Entry(input_frame, width=50)
    angle_entry.grid(row=2, column=1, sticky="EW")

    # Result display
    ttk.Label(input_frame, text="Results:").grid(row=3, column=0, sticky="W", pady=(10, 0))
    results_text = tk.StringVar()
    results_label = ttk.Label(input_frame, textvariable=results_text, foreground="blue", justify="left")
    results_label.grid(row=3, column=1, sticky="W", pady=(10, 0))

    # Buttons
    button_frame = ttk.Frame(root, padding="10")
    button_frame.grid(row=1, column=0, sticky="EW")
    ttk.Button(button_frame, text="Calculate & Plot", command=on_calculate).grid(row=0, column=0, padx=5, pady=5)
    ttk.Button(button_frame, text="Exit", command=root.destroy).grid(row=0, column=1, padx=5, pady=5)

    # Run the GUI
    root.mainloop()

# Launch GUI
get_gui_input()
