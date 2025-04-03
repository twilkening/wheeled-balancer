import numpy as np
import matplotlib.pyplot as plt

# Load the data from the file
data = np.loadtxt("src\simout.txt", skiprows=1, delimiter=',')  # Adjust delimiter if necessary (e.g., delimiter=',')

# Assuming the file has two columns: time and theta
time = data[:, 0]  # First column: time
theta = data[:, 1]  # Second column: theta
omega = data[:, 2]  # Third column: omega (angular velocity)

# # Plot the data
# plt.figure(figsize=(10, 6))
# plt.plot(time, theta, label="Theta (rad)", color="blue")
# plt.xlabel("Time (s)")
# plt.ylabel("Theta (rad)")
# plt.title("Inverted Pendulum Simulation")
# plt.legend()
# plt.grid()
# plt.show()

# Plot both theta and angular velocity
plt.figure(figsize=(10, 6))
plt.plot(time, theta, label="Theta (rad)", color="blue")
plt.plot(time, omega, label="Angular Velocity (rad/s)", color="red")
plt.xlabel("Time (s)")
plt.ylabel("Values")
plt.title("Inverted Pendulum Simulation")
plt.legend()
plt.grid()
plt.show()