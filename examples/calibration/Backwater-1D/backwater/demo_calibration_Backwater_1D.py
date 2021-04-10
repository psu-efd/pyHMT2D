import numpy as np
from scipy import interpolate


import pyHMT2D

from matplotlib import pyplot as plt

plt.rc('text', usetex=False)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def demo_backwater_1d_model_data():
    #load Backwater_1D_Data configuration data
    my_backwater_1D_data = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Data("calibration.json")

    #create a Backwater_1D_Model object
    my_backwater_1D_model = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Model()

    #set the simulation case in the Backwater_1D_Model object
    my_backwater_1D_model.set_simulation_case(my_backwater_1D_data)

    #run the Backwater_1D_Model model
    my_backwater_1D_model.run_model()

    #sample some results as manufactured solution
    xmin = np.min(my_backwater_1D_data.gridx)
    xmax = np.max(my_backwater_1D_data.gridx)

    #make some sampling points
    sampling_points = np.linspace(xmin+100,xmax-100,4)

    #build the WSE interpolator
    f_WSE = interpolate.interp1d(my_backwater_1D_data.gridx, my_backwater_1D_data.WSE)

    #bed elevation
    bedElev = my_backwater_1D_data.WSE - my_backwater_1D_data.waterDepth

    #sample the WSE
    sampled_WSE = f_WSE(sampling_points)

    #sampled the x-velocity
    sampled_U = my_backwater_1D_data.specificDischarge/sampled_WSE

    #output the true solution and sampled WSE and velocity
    fid = open("sampled_wse.csv", "w")
    fid.write("Name,x,y,wse\n")
    for pointI in range(len(sampling_points)):
       fid.write("Point%d, %f, 0, %f\n" % (pointI,sampling_points[pointI],sampled_WSE[pointI]))
    fid.close()

    fid = open("sampled_velocity.csv", "w")
    fid.write("Name,x,y,angle,velocity\n")
    for pointI in range(len(sampling_points)):
       fid.write("Point%d, %f, 0, 90, %f\n" % (pointI,sampling_points[pointI],sampled_U[pointI]))
    fid.close()

    #plot result
    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    negX = -(my_backwater_1D_data.gridx - my_backwater_1D_data.startx) / 1000
    #print("negX, WSE = ", negX, my_backwater_1D_data.WSE)

    wse_plot = ax.plot(negX[::50], my_backwater_1D_data.WSE[::50], 'g-', label='WSE')
    bedElev_plot = ax.plot(negX[::50], bedElev[::50], 'k-', label='Bed')

    #plot the sampled WSE
    ax.scatter(-(sampling_points-my_backwater_1D_data.startx)/1000, sampled_WSE, marker='o', facecolor='none',
               edgecolor='blue', s=20, linewidths=1)

    ax.set_xlabel('x (km)', fontsize=16)
    ax.set_ylabel('Elevation (m)', color='g', fontsize=16)

    #plot Manning's n
    n_plot = ax2.plot(negX[::50], my_backwater_1D_data.gridManningN[::50], 'b--', label='Manning\'s $n$')
    ax2.set_ylabel('Manning\'s $n$', color = 'b', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=12)

    #plt.title('Backwater curve')

    ax.set_xlim(-11, 1)
    ax.set_ylim(-0.3, 3)
    ax2.set_ylim(0.03,0.06)

    # add all plots' label
    lns = wse_plot + bedElev_plot + n_plot
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='center left', fontsize=14, frameon=False)

    fig.tight_layout()

    #save figure
    plt.savefig("backwater_1d_calibration_example.png", dpi=300, bbox_inches='tight', pad_inches=0)

    plt.show()


def demo_calibrator():
    """ Demo the Calibrator class

    Returns
    -------

    """

    my_calibrator = pyHMT2D.Calibration.Calibrator("calibration.json")

    my_calibrator.calibrate()


if __name__ == "__main__":

    #run the Backwater_1D model and sample some data for calibration later (manufactured solution)
    #demo_backwater_1d_model_data()

    #run the calibration
    demo_calibrator()

    print("All done!")