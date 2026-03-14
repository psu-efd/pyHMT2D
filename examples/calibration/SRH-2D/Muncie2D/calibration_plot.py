
import numpy as np
import math
import matplotlib

import matplotlib.pyplot as plt
from matplotlib import animation

plt.rc('text', usetex=False)  # allow the use of Latex for math expressions and equations
plt.rc('font', family='serif')  # specify the default font family to be "serif"


def plot_comparison():
    """ Plot simulation vs. measurement

    Returns
    -------

    """

    #load the measurement data for wse
    data_m_wse = np.genfromtxt("sampled_wse.csv", dtype=None, delimiter=",",
                               encoding = "utf-8", skip_header=1)

    name_m_wse = []
    x_m_wse = []
    y_m_wse = []
    wse_m = []

    for pointI in range(len(data_m_wse)):
        name_m_wse.append(data_m_wse[pointI][0])
        x_m_wse.append(data_m_wse[pointI][1])
        y_m_wse.append(data_m_wse[pointI][2])
        wse_m.append(data_m_wse[pointI][3])


    #load the simulation data for wse
    data_s_wse = np.genfromtxt("PointSimulation_wse.csv", dtype=None, delimiter=",",
                               encoding = "utf-8",skip_header=1)

    name_s_wse = []
    x_s_wse = []
    y_s_wse = []
    wse_s = []

    for pointI in range(len(data_s_wse)):
        name_s_wse.append(data_s_wse[pointI][0])
        x_s_wse.append(data_s_wse[pointI][1])
        y_s_wse.append(data_s_wse[pointI][2])
        wse_s.append(data_s_wse[pointI][3])


    #plot
    fig, axs = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

    axs[0].scatter(wse_m, wse_s, marker='o', facecolor='none', edgecolor='blue', s=30)

    #draw a diagonal line
    axs[0].plot([930,960],[930,960], color = 'k', linewidth = 1)

    axs[0].set_xlabel("Measured WSE (ft)", fontsize=16)
    axs[0].set_ylabel("Simulated WSE (ft)", fontsize=16)

    axs[0].set_xlim(930, 960)
    axs[0].set_ylim(930, 960)

    axs[0].set_aspect(1)

    # change the fontsize of major and minor ticks label
    axs[0].tick_params(axis='both', which='major', labelsize=12)
    axs[0].tick_params(axis='both', which='minor', labelsize=10)

    #plot velocity vector comparison

    #load the measurement data for velocity
    data_m_vel = np.genfromtxt("sampled_velocity.csv", dtype=None, delimiter=",",
                               encoding = "utf-8", skip_header=1)

    name_m_vel = []
    x_m_vel = []
    y_m_vel = []
    angle_m = []
    mag_m = []
    u_m = []
    v_m = []

    for pointI in range(len(data_m_vel)):
        name_m_vel.append(data_m_vel[pointI][0])
        x_m_vel.append(data_m_vel[pointI][1])
        y_m_vel.append(data_m_vel[pointI][2])
        angle_m.append(data_m_vel[pointI][3])
        mag_m.append(data_m_vel[pointI][4])

        u_m.append(mag_m[pointI] * np.sin(np.deg2rad(angle_m[pointI])))
        v_m.append(mag_m[pointI] * np.cos(np.deg2rad(angle_m[pointI])))

    #load the simulation data for velocity
    data_s_vel = np.genfromtxt("PointSimulation_velocity.csv", dtype=None, delimiter=",",
                               encoding = "utf-8", skip_header=1)

    name_s_vel = []
    x_s_vel = []
    y_s_vel = []
    angle_s = []
    mag_s = []
    u_s = []
    v_s = []

    for pointI in range(len(data_s_vel)):
        name_s_vel.append(data_s_vel[pointI][0])
        x_s_vel.append(data_s_vel[pointI][1])
        y_s_vel.append(data_s_vel[pointI][2])
        angle_s.append(data_s_vel[pointI][3])
        mag_s.append(data_s_vel[pointI][4])

        u_s.append(mag_s[pointI] * np.sin(np.deg2rad(angle_s[pointI])))
        v_s.append(mag_s[pointI] * np.cos(np.deg2rad(angle_s[pointI])))

    #Plotted in an array on a [0,1]x[0,1]canvas.

    #number of point vectors
    numPoints = len(data_s_vel)

    #The array has nColumns columns.
    nColumns = 3
    nRows =  math.ceil(numPoints/nColumns)

    dx = 1.0/nColumns
    dy = 1.0/nRows

    x_vel = []
    y_vel = []

    #the coordinates of each vector
    for pointI in range(numPoints):
        iRow = math.floor(pointI/nColumns)
        iColumn = pointI % nColumns

        x_vel.append(iColumn * dx + dx/2)
        y_vel.append(iRow * dy + dy / 2)


    #plot measured and simulation velocity
    axs[1].quiver(x_vel, y_vel, u_m, v_m, scale=50, color='r')

    q = axs[1].quiver(x_vel, y_vel, u_s, v_s, scale=50, color='b',
                      linestyle='dashed',facecolor='none', linewidth=1,
                        width=0.0001, headwidth=200, headlength=300)

    #add point names
    for pointI in range(numPoints):
        #add some offset to the position of text so it does not cover the vector
        #offset in the opposite direction of velocity
        offsetx = 0.05 * np.sin(np.deg2rad(angle_s[pointI]))
        offsety = 0.05 * np.cos(np.deg2rad(angle_s[pointI]))
        axs[1].text(x_vel[pointI]-offsetx, y_vel[pointI]-offsety, name_s_vel[pointI])

    #add a key
    axs[1].quiverkey(q, X=0.1, Y=0.9, U = 5, label='5 ft/s', color='k')

    axs[1].set_xlim(0, 1)
    axs[1].set_ylim(0, 1)

    axs[1].set_aspect(1)

    # change the fontsize of major and minor ticks label
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    #axs[1].axis('off')

    #save figure
    plt.savefig("Muncie_2d_calibration_comparison.png", dpi=300, bbox_inches='tight', pad_inches=0)

    plt.show()

if __name__ == "__main__":

    #plot WSE and velocity vector comparison between measurement and simulation
    plot_comparison()

    print("All done!")