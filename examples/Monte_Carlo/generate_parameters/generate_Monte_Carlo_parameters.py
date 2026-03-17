"""
Generate parameters for Monte Carlo simulations. 

It demonstrates how to generate parameters for Monte Carlo simulations with a truncated normal distribution for Manning's n.
"""

import numpy as np
from scipy.stats import truncnorm
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

def generate_ManningN_from_distribution(n_min, n_max, n_mean, n_std, nSamples):
    """
    Generate Manning's n values from a distribution. Currently the distribution is a truncated normal distribution.

    Parameters
    ----------
    n_min : float
        minimum value of Manning's n
    n_max : float
        maximum value of Manning's n
    n_mean : float
        mean value of Manning's n
    n_std : float
        std of Manning's n
    nSamples : int
        number of samples to generate

    Returns
    -------

    """

    #convert the user specified min, mean, max, and std to [a, b] on standard normal distribution
    a, b = (n_min - n_mean) / n_std, (n_max - n_mean) / n_std
    #print("a,b =", a, b)

    #call scipy's truncnorm to generate the samples (on the standard normal distribution
    samples = truncnorm(a=a, b=b, scale=1.0).rvs(size=nSamples)

    #convert the random values back to the normal distribution scale
    samples = samples * n_std + n_mean

    #print("samples = ", samples)
    print("min, max, and mean of samples = ", np.min(samples), np.max(samples), np.mean(samples))

    #make a plot to show the sample distribution
    plt.hist(samples, 11, facecolor='gray', alpha=0.75)

    # set the limit for the x and y axes
    plt.xlim([n_min, n_max])
    #plt.ylim([0, 45])

    # set x and y axes label and font size
    plt.xlabel('Manning\'s n', fontsize=16)
    plt.ylabel('Count', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=12)

    # set axis label format
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}'))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    # show title and set font size
    plt.title('Histogram of sampled Manning\'s n values', fontsize=16)

    # save the plot to a file
    plt.savefig('sampledManningN_histogram.png', dpi=300, bbox_inches='tight')

    # show legend, set its location, font size, and turn off the frame
    #plt.legend(loc='lower left', fontsize=14, frameon=False)
    #plt.show()

    return samples

def sample_ManningN():
    """

    Returns
    -------

    """

    print("Generating Manning's n values from a truncated normal distribution ...")

    #define parameters
    n_min = 0.03
    n_mean = 0.04
    n_max = 0.05
    n_std = 0.005
    nSamples = 100

    #generate random samples for Manning's n
    samples = generate_ManningN_from_distribution(n_min, n_max, n_mean, n_std, nSamples)

    #print the first 5 samples
    print("first 5 samples = ", samples[:5])

    #save the sampled Manning's n values for record
    date_time = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")    #take the current date and time
    fileName_ManningN = "sampledManningN_"+date_time+".dat"
    np.savetxt(fileName_ManningN, samples, delimiter=',')

    return samples

if __name__ == "__main__":
    #sample Manning's n values
    samples = sample_ManningN()
   
    print("All done!")

