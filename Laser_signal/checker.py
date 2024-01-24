import matplotlib.pyplot as plt
import numpy as np

# Parameters (in milliseconds)
expossure = 40  # Exposure time
dur = 500
freq = 4
dt = 3
laser_on = 2

flash = 1  # Example values
exp_flash = expossure - flash - laser_on
laser_flash = laser_on - flash
loop = 1000/freq - flash - dt - expossure
# dur = 5000  # Duration for which we want to simulate

# Initialize time and pin states
time = 0
cameraPin_state = []
laserPin_state = []
time_points = []

while time < dur:
    # Camera on
    cameraPin_state.append(1)
    laserPin_state.append(0)
    time_points.append(time)
    time += flash

    # Camera off, delay for exp_flash
    cameraPin_state.append(0)
    laserPin_state.append(0)
    time_points.append(time)
    time += exp_flash

    # Laser on
    cameraPin_state.append(0)
    laserPin_state.append(1)
    time_points.append(time)
    time += dt

    # Camera on, delay for flash
    cameraPin_state.append(1)
    laserPin_state.append(1)
    time_points.append(time)
    time += flash

    # Camera off, delay for laser_flash
    cameraPin_state.append(0)
    laserPin_state.append(1)
    time_points.append(time)
    time += laser_flash

    # Laser off, delay for loop
    cameraPin_state.append(0)
    laserPin_state.append(0)
    time_points.append(time)
    time += loop

# Plotting
plt.figure(figsize=(10, 4))
plt.yticks([0, 1])

plt.step(time_points, cameraPin_state, where='post', label='Camera Pin')
plt.step(time_points, laserPin_state, where='post', label='Laser Pin')
plt.xlabel('Time (ms)')
plt.ylabel('State')
plt.title('Simulation of Camera and Laser Pin States Over Time')
plt.legend(loc='upper left')
plt.grid(True)
plt.show()
