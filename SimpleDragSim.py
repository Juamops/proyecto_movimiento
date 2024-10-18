import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.stats import norm
from matplotlib.animation import FuncAnimation

# Simulation Characteristics
timestep = 0.05
simulation_time = 20

# Rock Characteristics in SI units
density = 2990
radius = 0.4
volume = (4/3) * pi * radius**3
mass = volume * density
surface_area = 4 * pi * radius**2
characteristic_length = volume / surface_area

# Environment characteristics in SI units
air_density = 1.182014
kinematic_viscosity = 15.52 * 10**(-6)
dynamic_viscosity = 18.37 * 10**(-6)
drag = 0.147
g = 9.81

x0 = 0
vx0 = 20
ax0 = -(drag * vx0**2) / mass

y0 = 0
vy0 = 20
ay0 = -((drag * vy0**2) / mass) - g

x_coord = [x0]
vx = [vx0]
ax = [ax0]

y_coord = [y0]
vy = [vy0]
ay = [ay0]
time = [0]

for i in range(1, round(simulation_time/timestep)):
    t = round((i * timestep) + timestep, 2)
    time.append(t)

    axn = -(drag * vx[i-1] * abs(vx[i-1])) / mass
    vxn = (ax[i-1] * timestep) + vx[i-1]
    xn = (0.5 * ax[i-1] * timestep**2) + (vx[i-1] * timestep) + (x_coord[i - 1])

    ayn = -((drag * vy[i-1] * abs(vy[i-1])) / mass) - g
    vyn = (ay[i-1] * timestep) + vy[i-1]
    yn = (0.5 * ay[i-1] * timestep**2) + (vy[i-1] * timestep) + (y_coord[i - 1])

    ax.append(axn)
    vx.append(vxn)
    x_coord.append(xn)

    ay.append(ayn)
    vy.append(vyn)
    y_coord.append(yn)

# plt.scatter(x_coord, y_coord)
# plt.show()

tau = mass / (drag * vx0)
vx_real = [vx0 / (1 + (t / tau)) for t in time]
difference = [estimate - real for estimate, real in zip(vx, vx_real)]

average_velocity = np.average([(x**2 + y**2)**0.5 for x, y in zip(vx, vy)])

rmse = (sum([diff**2 for diff in difference]) / len(difference)) ** 0.5
mae = (sum([abs(diff) for diff in difference]) / len(difference))
mape = sum([abs(diff) / real for diff, real in zip(difference, vx_real)]) / len(difference)

reynolds = (air_density * average_velocity * characteristic_length) / dynamic_viscosity

print(f'Reynolds Number: {reynolds}')
print(f'Average Velocity: {average_velocity}')
print(f'Root Mean Square Error: {rmse}')
print(f'Mean Absolute Error: {mae}')
print(f'Mean Absolute Percentual Error: {mape * 100:.2f}%')



# Crear la figura para la animación
# fig, ax = plt.subplots()
# scat = ax.scatter(x_coord[0], y_coord[0], c="b", s=5, label=f'vy0={vy0}, vx0={vx0} m/s')
# ax.set(xlim=[0, 60], ylim=[-20, 40], xlabel='X (m)', ylabel='Y (m)')
# ax.legend()
#
#
# def update(frame):
#     # for each frame, update the data stored on each artist.
#     x = x_coord[:frame]
#     y = y_coord[:frame]
#     # update the scatter plot:
#     data = np.stack([x, y]).T
#     scat.set_offsets(data)
#     # update the line plot:
#     return (scat,)


# ani = FuncAnimation(fig=fig, func=update, frames=len(x_coord), interval=timestep * 10)
# plt.show()

#plt.plot(time, y_coord, label='position')
plt.plot(time, vx, label='velocity (euler)')
plt.plot(time, vx_real, '--', label='velocity (real)')
#plt.plot(time, ax, label='acceleration')
#plt.axhline(y=0, color='black')
plt.legend()
plt.show()