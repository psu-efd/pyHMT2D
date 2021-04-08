import numpy as np
import os

import matplotlib.pyplot as plt

import pyHMT2D

plt.rc('text', usetex=False)   #whether allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def create_exit_rating_curve(my_srh_2d_data):
    """ Create a rating curve

    Parameters
    ----------
    my_srh_2d_data

    Returns
    -------

    """

    #output boundary profile
    my_srh_2d_data.srhgeom_obj.output_nodeString_line_coordinates(nodeStringID=1,
                                    nodeStringFileName="exit_boundary_profile.dat", dir='')

    #output boundary Manning's n profile
    my_srh_2d_data.output_boundary_manning_n_profile(nodeStringID=1,
                                    nodeStringFileName="exit_boundary_ManningN.dat", dir='')

    #read in the station and elevation of the boundary profile
    station_profile, x, y, z = np.loadtxt('exit_boundary_profile.dat', delimiter=',', unpack=True, skiprows=1)

    #read in the station and Manning's n of the boundary profile
    station_ManningN, ManningN = np.loadtxt('exit_boundary_ManningN.dat', delimiter=',', unpack=True, skiprows=1)

    #create a rating curve for a boundary
    overboard = 10.0
    #ManningN = 0.04   #need to check this. How to deal with composite Manning n?
    slope = 0.00064
    number_of_rc_points = 30
    units = 'EN'

    stage, Q, Area, Pwet, station_resample, z_resample = pyHMT2D.Misc.generate_rating_curve_based_on_Mannings_equation(
                                                         station_profile, z, overboard, station_ManningN, ManningN,
                                                         slope, number_of_rc_points, nResample=201,units=units)

    #plot all: cross-section profile, stage vs. Q, stage vs. area, stage vs. wetted perimeter

    fig, axs = plt.subplots(2, 2, constrained_layout=True)
    #fig.tight_layout()

    if units == 'EN':
        Station_text = 'Station (ft)'
        Elevation_text = 'Elevation (ft)'
        Discharge_text = 'Discharge (cfs)'
        Area_text = 'Area (ft$^2$)'
        PWet_text = 'Wetted perimeter (ft)'
        Stage_text = 'Stage (ft)'
    else:
        Station_text = 'Station (m)'
        Elevation_text = 'Elevation (m)'
        Discharge_text = 'Discharge (m$^3$)'
        Area_text = 'Area (m$^2$)'
        PWet_text = 'Wetted perimeter (m)'
        Stage_text = 'Stage (m)'

    #plot the cross-section profile
    sub1_ax2 = axs[0, 0].twinx()
    axs[0, 0].plot(station_resample, z_resample, label='Profile')
    #axs[0, 0].scatter(station_profile, z, c='red', s=5, alpha=0.3)
    sub1_ax2.plot(station_ManningN, ManningN, linestyle='--', linewidth=1.5, color='m', label='Manning\'s $n$')
    axs[0, 0].plot(np.nan, linestyle='--', linewidth=1.5, color='m', label='Manning\'s $n$')  # Make an agent in ax
    #axs[0, 0].set_title('Cross section profile')
    axs[0, 0].set(xlabel=Station_text, ylabel=Elevation_text)
    sub1_ax2.set(ylabel="Manning's $n$")
    sub1_ax2.set_ylim(0.02, 0.08)
    axs[0, 0].legend(loc="lower left", fontsize=10,frameon=False)

    #plot rating curve stage vs. Q
    axs[0, 1].plot(Q, stage)
    #axs[0, 1].set_title('Stage vs. Q')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set(xlabel=Discharge_text, ylabel=Stage_text)

    #plot stage vs. area
    axs[1, 0].plot(Area, stage)
    #axs[1, 0].set_title('Area vs. Stage')
    axs[1, 0].set(xlabel=Area_text, ylabel=Stage_text)

    #plot stage vs. wetted perimenter
    axs[1, 1].plot(Pwet, stage)
    #axs[1, 1].set_title('Stage vs. Wetted perimeter')
    axs[1, 1].set(xlabel=PWet_text, ylabel=Stage_text)

    # show legend, set its location, font size, and turn off the frame
    #plt.legend(loc='lower right', fontsize=14, frameon=False)

    #save the plot to file
    plt.savefig("exit_rating_curve.png", dpi=300, bbox_inches='tight', pad_inches=0)

    plt.show()

    # save the rating curve to file
    fid = open("exit_rating_curve.dat", "w")

    fid.write("RATING_CURVE\n")
    fid.write("// \n")
    fid.write("// Rating curve at the exit \n")
    fid.write("// Q(cfs)  WSE(ft)\n")
    fid.write("// \n")

    for i in range(len(Q)):
        fid.write("%f %f \n" % (Q[i], stage[i]))

    fid.close()

def convert_results_to_vtk():
    """ Read the resutls and output to VTK

    Returns
    -------

    """

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center (currently it is all nodal)
    bNodal = False

    #my_srh_2d_data.readSRHXMDFFile("Muncie_XMDFC.h5", bNodal)

    #outputXMDFDataToVTK(bNodal, timeStep=-1,lastTimeStep=False, dir=''):
    #my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='')

    # User specified SRH-2D result in SRH (point) or SRHC (cell center) format
    #srhFileName = 'backwater_curve_SRHC24.dat'

    # whehter it is cell center data or at point (need to be checked by user)
    #bCellData = True

    # Read SRH-2D result file
    #resultVarNames, resultData = my_srh_2d_data.readSRHFile(srhFileName)
    # print(resultVarNames, resultData)

    # output SRH-2D result to VTK
    #srhName = os.path.splitext(srhFileName)[0]
    #vtkFileName = srhName + ".vtk"
    #print("vtk file name = ", vtkFileName)

    #my_srh_2d_data.outputVTK(vtkFileName, resultVarNames, resultData, bCellData)


def main():
    """ Testing SRH_2D_Data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie.srhhydro")

    #output the flat 2D mesh
    #my_srh_2d_data.output_flat_mesh_to_vtk('Duck_pond_flat_mesh.vtk')

    #create a rating curve
    create_exit_rating_curve(my_srh_2d_data)

    #convert the results to VTK format
    #convert_results_to_vtk()


    print("All done!")


if __name__ == "__main__":
    main()