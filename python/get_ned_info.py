#!/usr/bin/env python

# Get the NED info for the detected dwarf candidates

# Python 2/3 compatibility
from __future__ import print_function # to use print() as a function in Python 2

try:
    input = raw_input # use 'input' function in both Python 2 and 3
except NameError:
    pass

# std lib
#from getpass import getpass

# 3rd party
#import pandas as pd
import numpy as np
#import pylab as plt
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from astropy import utils, io, convolution, stats, wcs
#from astropy.visualization import make_lupton_rgb
from astropy import units as u
#from astropy.stats import median_absolute_deviation as mad
from astropy.table import Table
#from scipy.interpolate import interp1d
#from scipy import arange, array, exp
#%matplotlib inline
#import healpy as hp
#from photutils import find_peaks, data_properties
#import logging
import socket
import os
import sys
import time
from astropy.io import fits

from astroquery.ned import Ned
from astropy import coordinates

# Data Lab
#from dl import authClient as ac, queryClient as qc, storeClient as sc, helpers

# Quiet the Astropy warnings
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


# Do a single NED query 
def getned(ra,dec,radius=0.1):
    
    try: 
        co = coordinates.SkyCoord(ra=ra, dec=dec,
                                  unit=(u.deg, u.deg), frame='fk4')
        result_table = Ned.query_region(co, radius=radius * u.deg)
    except:
        result_table = None
    
    return result_table


if __name__ == "__main__":

    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    t0 = time.time()

    print(sys.argv)

    # Not enough inputs
    n = len(sys.argv)
    if n < 2:
        print("Syntax - get_ned_info.py catalogfile outfile")
        sys.exit()

    # Inputs
    infile = sys.argv[1]
    outfile = sys.argv[2]

    # Load the input file
    cat = Table(fits.getdata(infile,1))
    ncat = len(cat)

    # Add columns
    cat['NED_NOBJ'] = 0
    cat['NED_NAME'] = '                                                     '
    cat['NED_RA'] = 999999.0
    cat['NED_DEC'] = 999999.0
    cat['NED_TYPE'] = '                                                     '
    cat['NED_REDSHIFT'] = 0.0

    # Loop over the peaks and download data
    for i in range(ncat):
        print(str(i+1)+' '+str(cat[i]['RA'])+' '+str(cat[i]['DEC']))
        ned_info = getned(cat[i]['RA'],cat[i]['DEC'],0.1)
        if ned_info is not None:
            # Sort by distance from queried position
            ned_info.sort('Distance (arcmin)')
            # Take closest one
            ned_info1 = ned_info[0]
            print(ned_info1)
            # Stuff information into the table
            cat[i]['NED_NOBJ'] = len(ned_info)
            cat[i]['NED_NAME'] = ned_info1['Object Name']
            cat[i]['NED_RA'] = ned_info1['RA(deg)']
            cat[i]['NED_DEC'] = ned_info1['DEC(deg)']
            cat[i]['NED_TYPE'] = ned_info1['Type']
            cat[i]['NED_REDSHIFT'] = ned_info1['Redshift']


    # Save the output file
    cat.write(outfile,'fits')
