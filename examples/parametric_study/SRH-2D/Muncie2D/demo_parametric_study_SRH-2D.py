import numpy as np
import csv
import vtk
from matplotlib import pyplot as plt

from pyHMT2D import Parametric_Study

from pyHMT2D import gVerbose

plt.rc('text', usetex=False)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_searchspace(x, parameter_names, title):
    v1 = 0
    v2 = 2

    fig, ax = plt.subplots()
    plt.plot(np.array(x)[:, v1], np.array(x)[:, v2], 'bo', label='samples')
    plt.plot(np.array(x)[:, v1], np.array(x)[:, v2], 'bo', markersize=10, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel(parameter_names[v1])
    #ax.set_xlim([-5, 10])
    ax.set_ylabel(parameter_names[v2])
    #ax.set_ylim([0, 15])
    plt.title(title)

    plt.show()

def demo_parametric_study():
    """ Demo the parametric study class

    Returns
    -------

    """

    my_parametric_study = Parametric_Study("parametric_study.json")

    samples, parameter_names = my_parametric_study.create_all_cases()

    plot_searchspace(samples, parameter_names, "Samples")


if __name__ == "__main__":

    #run the calibration
    demo_parametric_study()

    print("All done!")