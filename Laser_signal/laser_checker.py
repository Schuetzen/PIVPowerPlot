import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

# Parameters (in milliseconds)
exposure = 40  # Exposure time
dur = 500
freq = 10
dt = 3
laser_on = 2

flash = 1  # Example values
exp_flash = exposure - flash - laser_on
laser_flash = laser_on - flash
loop = 1000 / freq - flash - dt - exposure

# Initialize time and pin states
cameraPin_state = []
laserPin_state = []
time_points = []
exposure_bars = []

# Function to calculate states
def calculate_states():
    time = 0
    while time < dur:
        yield time, 1, 0  # Camera on
        time += flash
        yield time, 0, 0  # Camera off, delay for exp_flash
        time += exp_flash
        yield time, 0, 1  # Laser on
        time += dt
        yield time, 1, 1  # Camera on, delay for flash
        time += flash
        yield time, 0, 1  # Camera off, delay for laser_flash
        time += laser_flash
        yield time, 0, 0  # Laser off, delay for loop
        time += loop

# Update function for the animation
def update(frame):
    time, camera_state, laser_state = frame
    time_points.append(time)
    cameraPin_state.append(camera_state)
    laserPin_state.append(laser_state)

    line1.set_data(time_points, cameraPin_state)
    line2.set_data(time_points, laserPin_state)

    # Manage exposure bars
    if camera_state == 1 and (len(time_points) == 1 or cameraPin_state[-2] == 0):
        bar = Rectangle((time, 0), exposure, 0.5, color='red', alpha=0.4)
        ax.add_patch(bar)
        exposure_bars.append(bar)

    return [line1, line2] + exposure_bars

# Set up the figure for plotting
fig, ax = plt.subplots(figsize=(10, 8))
line1, = ax.step([], [], where='post', label='Camera Pin', linestyle='--')
line2, = ax.step([], [], where='post', label='Laser Pin')

# Add a dummy rectangle for exposure time legend
dummy_exposure_bar = Rectangle((0, 0), 0.5, 0.5, color='red', alpha=0.4, label='Exposure Time')
ax.add_patch(dummy_exposure_bar)  # Add dummy bar to be able to create legend

ax.set_xlim(-1, dur)
ax.set_ylim(-0.1, 1.1)
ax.set_yticks([0, 1])  # Set y-axis to only show 0 and 1

ax.set_xlabel('Time (ms)')
ax.set_ylabel('State')
ax.set_title('Simulation of Camera and Laser Pin States Over Time')

# Create the legend
ax.legend(handles=[line1, line2, dummy_exposure_bar], loc='upper left')

ax.grid(False)

# Create the animation
ani = FuncAnimation(fig, update, frames=calculate_states(), blit=False, repeat=False)

plt.show()
