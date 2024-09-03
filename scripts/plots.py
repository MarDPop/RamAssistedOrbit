import matplotlib.pyplot as plt
import subprocess
import numpy as np

rad2deg = 180.0/np.pi

def quaternion_to_matrix(q):
    """
    Converts a quaternion (w, x, y, z) to a 3x3 rotation matrix.
    :param q: Quaternion (numpy array or list).
    :return: 3x3 rotation matrix.
    """
    x, y, z, w = q
    return np.array([[1 - 2 * (y**2 + z**2), 2 * (x*y - w*z), 2 * (x*z + w*y)],
                     [2 * (x*y + w*z), 1 - 2 * (x**2 + z**2), 2 * (y*z - w*x)],
                     [2 * (x*z - w*y), 2 * (y*z + w*x), 1 - 2 * (x**2 + y**2)]])

# subprocess.run(["../bin/sim.exe"], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

track_x = []
track_y = []
track_t = []
track_u = []
track_v = []
track_speed = []
with open("../output/track.dat","r") as trackfile:
    lines = trackfile.readlines()

    for line in lines:
        arr = line.split()
        track_t.append(float(arr[0]))
        track_x.append(float(arr[1]))
        track_y.append(float(arr[2]))
        track_u.append(float(arr[3]))
        track_v.append(float(arr[4]))
        track_speed.append(np.sqrt(track_u[-1]*track_u[-1] + track_v[-1]*track_v[-1]))
        
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(track_x,track_y)
ax1.set_aspect('equal')
ax2.plot(track_t,track_speed)

#%%
rad2deg = 57.12
pos = []
vel = []
quat = []
angVel = []
mass = []
times = []
mach = []
flightAngle = []
gForce = []
altitude = []
elevator = []
thrust = []
with open("../output/result.dat","r") as resultfile:
    lines = resultfile.readlines()

    for line in lines:
        arr = line.split()
        times.append(float(arr[0]))
        pos.append(np.array([float(arr[1]),float(arr[2]),float(arr[3])]))
        vel.append(np.array([float(arr[4]),float(arr[5]),float(arr[6])]))
        quat.append(np.array([float(arr[7]),float(arr[8]),float(arr[9]),float(arr[10])]))
        angVel.append(np.array([float(arr[11]),float(arr[12]),float(arr[13])]))
        mass.append(float(arr[14]))
        mach.append(float(arr[15]))
        flightAngle.append(float(arr[16]))
        gForce.append(float(arr[17]))
        altitude.append(float(arr[18]))
     
        if len(arr) > 20:
            elevator.append(float(arr[19])*rad2deg)
            thrust.append(float(arr[20]))
        else:
            elevator.append(0)
            thrust.append(0)
AoA = []
speed = []
pitch = []
for i in range(len(times)):
    rotm = quaternion_to_matrix(quat[i])
    up = pos[i] / np.linalg.norm(pos[i])
    z = np.dot(rotm[:,0],up)
    p = np.arcsin(z)
    pitch.append(p*rad2deg)
    s = np.linalg.norm(vel[i])
    speed.append(s)
    velz = np.dot(vel[i],up)
    angleOfAttack = p - np.arcsin(velz/s)
    AoA.append(angleOfAttack*rad2deg)
    
pitch[0] = 0
        
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
ax1.plot(times,altitude, label='Altitude')
ax1.set_ylabel('meters')
ax11 = ax1.twinx()
ax11.plot(times, elevator, color = 'tab:red', label='elevator deflection')
ax11.set_ylabel('Degrees')
ax1.set_title("Altitude")
ax1.grid()
ax1.legend();
ax11.legend();

ax2.plot(times, mach)
ax22 = ax2.twinx()
ax22.plot(times, speed, color = 'tab:red')
ax22.set_ylabel('Speed')
ax2.set_title("Mach")
ax2.grid()
ax2.legend()

ax3.plot(times, mass)
ax33 = ax3.twinx()
ax33.plot(times, thrust, color = 'tab:red', label='thrust')
ax3.set_title("Mass")
ax3.grid()
ax33.legend()

ax4.plot(times, pitch, label='pitch')
ax4.set_title("Pitch")
ax44 = ax4.twinx()
ax44.plot(times, AoA, color = 'tab:red', label='AoA')
ax4.grid()
ax4.legend()
ax44.legend()

plt.show()