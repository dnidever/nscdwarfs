#!/usr/bin/env python

# 1) pull data from object table for 2x2 deg
#       star/galaxy separation, color cut
#    pull coverage information from coverage table



# 2) color cut in g-r

# 3) make density map

# 4) correct for covering fraction, set missing data to average value

# 5) use luminosity function to correct for depth
#     try to to calculate this empircally, density vs. depth
#     or query Besancon at that position and figure out luminosity function

# 6) mexican hat filter

# 7) run peak finder, store peaks
#      ra, dec, significance, density above the background, shape parameter

# 8) loop over peaks, do narrow cone search (size depends on shape)
#     everything from catalog around select portion of sky

# 9) visualization to look at CMD and possibly an image

# 10) machine learning to weed out the contaminants, large galaxies, galaxy clusters, etc.

# Python 2/3 compatibility
from __future__ import print_function # to use print() as a function in Python 2

try:
    input = raw_input # use 'input' function in both Python 2 and 3
except NameError:
    pass

# std lib
from getpass import getpass

# 3rd party
import pandas as pd
import numpy as np
import pylab as plt
import matplotlib
from astropy import utils, io, convolution, stats, wcs
from astropy.visualization import make_lupton_rgb
from astropy import units as u
from astropy.stats import median_absolute_deviation as mad
#%matplotlib inline
import healpy as hp
from photutils import find_peaks
import logging
import socket
import os
import sys
import time

# Data Lab
from dl import authClient as ac, queryClient as qc, storeClient as sc, helpers

#Simple Image Access (SIA) service
from pyvo.dal import sia
DEF_ACCESS_URL = "http://datalab.noao.edu/sia"
svc = sia.SIAService(DEF_ACCESS_URL)

# Quiet the Astropy warnings
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

# A function to retrieve data from a point on the sky
def getData (ra,dec,radius=1.0,columns='*'):

    query_template = \
    """SELECT {0:s} FROM nsc_dr1.object
       WHERE q3c_radial_query(ra,dec,{1:f},{2:f},{3:f})"""

    query = query_template.format(columns,ra,dec,radius)
    print(query)
    
    try:
        result = qc.query(token,sql=query) # by default the result is a CSV formatted string
    except Exception as e:
        print(e.message)
        
    df = helpers.convert(result,'pandas')
    
    return df

# A function to retrieve data from a point on the sky
def getBlueStarData (ra,dec,radius=1.0,columns='ra,dec,gmag,rmag',colcutlo=None,colcuthi=None,classcut=None,fwhmcut=None,errcut=None):

    query_template = \
    """SELECT {0:s} FROM nsc_dr1.object
       WHERE q3c_radial_query(ra,dec,{1:f},{2:f},{3:f})"""
    #   (gmag-rmag)>({4:f}) and (gmag-rmag)<{5:f} and class_star>{6:f} and
    #   fwhm<{7:f} and gerr<{8:f} and rerr<{8:f}

    query = query_template.format(columns,ra,dec,radius,colcutlo,colcuthi,classcut,fwhmcut,errcut)
    if colcutlo is not None: query+=" and (gmag-rmag)>("+"{0:f}".format(colcutlo)+")"
    if colcuthi is not None: query+=" and (gmag-rmag)<"+"{0:f}".format(colcuthi)
    if classcut is not None: query+=" and class_star>"+"{0:f}".format(classcut)
    if fwhmcut is not None: query+=" and fwhm<"+"{0:f}".format(fwhmcut)
    if errcut is not None: query+=" and gerr<"+"{0:f}".format(errcut)
    if errcut is not None: query+=" and rerr<"+"{0:f}".format(errcut)
    print(query)
    
    try:
        result = qc.query(token,sql=query) # by default the result is a CSV formatted string
    except Exception as e:
        print(e.message)
        
    df = helpers.convert(result,'pandas')
    
    return df

# A function to retrieve data from a point on the sky
def getCovData (ra,dec,radius=1.0,columns='*'):

    query_template = \
    """SELECT {0:s} FROM nsc_dr1.coverage
       WHERE q3c_radial_query(ra,dec,{1:f},{2:f},{3:f})"""

    query = query_template.format(columns,ra,dec,radius)
    print(query)
    
    try:
        result = qc.query(token,sql=query) # by default the result is a CSV formatted string
    except Exception as e:
        print(e.message)
        
    df = helpers.convert(result,'pandas')
    
    return df

# A Mexican-hat convolution filter
def dwarf_filter (ra,dec,fwhm_small=2.0,fwhm_big=20):

    # Based on Koposov et al. (2008).
    # Code by Ken Mighell and Mike Fitzpatrick.
    # Minor edits by Rorbert Nikutta.
    
    x, y = ra, dec

    # Information about declination (y) [degrees]
    ymean = (y.min() + y.max()) / 2.0
    ydiff_arcmin = (y.max() - y.min()) * 60.0 # convert from degrees to arcmin

    # Information about right ascension (x) [degrees in time]:
    xdiff = x.max() - x.min() # angular separation [degrees (time)] 
    xmean = (x.min() + x.max())/2.0

    # convert from degrees in time to separation in angular degrees:
    xdiff_angular = (x.max() - x.min()) * np.cos(ymean*(np.pi/180.0))

    # convert from degress to arcmin
    xdiff_angular_arcmin = xdiff_angular * 60.0 

    # Get the number of one-arcmin pixels in the X and Y directions:
    nx = np.rint (xdiff_angular_arcmin).astype('int')
    ny = np.rint (ydiff_arcmin).astype('int')
    
    # Create a two-dimensional histogram of the raw counts:
    Counts, xedges, yedges  = np.histogram2d (x, y, (nx,ny) )
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    raw_hist = np.rot90(Counts).copy() # hack around Pythonic weirdness

    # Make the small and big Gaussian kernels with a standard deviation
    # of the given FWHM in arcmin^2 pixels.
    kernel_small = convolution.Gaussian2DKernel(fwhm_small/2.35,factor=1)
    kernel_big = convolution.Gaussian2DKernel(fwhm_big/2.35,factor=1)

    # Compute the differential convolution kernels.
    conv_big = convolution.convolve (raw_hist, kernel_big)
    conv_small = convolution.convolve (raw_hist, kernel_small)
    conv_delta = conv_small - conv_big
    delta = conv_delta.copy()

    # Compute statistics and the floor
    mean = np.mean (delta, dtype='float64')
    sigma = np.std (delta, dtype='float64')
    sigmaRaw = np.std(raw_hist,dtype='float64')
    median = np.median (delta)                       # not used
    floor = mean 

    # Clip to specified limits.
    clipped = delta.copy()
    clipped[ delta < floor ] = floor

    # Return the computed fields.
    return raw_hist, extent, delta, clipped, sigma

# A little function to download the deepest stacked images
#   adapted from R. Nikutta
def download_deepest_image(ra,dec,fov=0.1,band='g'):
    imgTable = svc.search((ra,dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).votable.to_table()
    print("The full image list contains", len(imgTable), "entries")
    
    sel0 = imgTable['obs_bandpass'].astype(str)==band
    sel = sel0 & ((imgTable['proctype'].astype(str)=='Stacked') & (imgTable['prodtype'].astype(str)=='image')) # basic selection
    Table = imgTable[sel] # select
    if (len(Table)>0):
        row = Table[np.argmax(Table['exptime'].data.data.astype('float'))] # pick image with longest exposure time
        url = row['access_url'].decode() # get the download URL
        print ('downloading deepest image...')
        image = io.fits.getdata(utils.data.download_file(url,cache=True,show_progress=False,timeout=120))

    else:
        print ('No image available.')
        image=None
        
    return image

# Multi panel image plotter
def plot_images(images,geo=None,panelsize=4,bands=list('gri'),cmap=matplotlib.cm.gray_r):
    n = len(images)
    if geo is None: geo = (n,1)
        
    fig = plt.figure(figsize=(geo[0]*panelsize,geo[1]*panelsize))
    for j,img in enumerate(images):
        ax = fig.add_subplot(geo[1],geo[0],j+1)
        if img is not None:
            print(img.min(),img.max())
            vmin = np.median(img)-2*np.std(img)
            vmax = np.median(img)+2*np.std(img)
            ax.imshow(img,origin='lower',interpolation='none',cmap=cmap,norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title('%s band' % bands[j])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

# Create a WCS for a tangent plane projection in our region
def get_wcs(ra,dec,image,fov=1.,unit='deg',projection=("RA---TAN","DEC--TAN")):
    npix = image.shape[0]
    crpix = npix/2 + 1
    cdelt = fov/float(npix)
    w = wcs.WCS(naxis=2)
    w.wcs.cunit = (unit,unit)
    w.wcs.crpix = (crpix,crpix)
    w.wcs.cdelt = np.array((-cdelt,cdelt))
    w.wcs.ctype = projection
    w.wcs.crval = (ra,dec) #coords.ra.to(unit).value, coords.dec.to(unit).value)
    return w

def plotpanel(axid,x,y,title='',xlim=(-1,2),ylim=(25.2,14)):
    ax = fig.addxs_subplot(axid)
    ax.scatter(x,y,marker='.',s=10, alpha=0.8)
    ax.set_xlabel(x.name)
    ax.set_ylabel(y.name)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title)


if __name__ == "__main__":

    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Version
    verdir = "v1"
    if len(sys.argv) > 4:
       version = sys.argv[4]
       verdir = version if version.endswith('/') else version+"/"

    # on thing/hulk use
    if (host == "thing") | (host == "hulk"):
        dir = "/dl1/users/dnidever/nsc/dwarfs/"+verdir
    # on gp09 use
    if (host == "gp09") | (host == "gp08") | (host == "gp07") | (host == "gp06") | (host == "gp05"):
        dir = "/net/dl1/users/dnidever/nsc/dwarfs/"+verdir
    if (host == "NideverMacBookPro"):
        dir = "/Users/nidever/datalab/nsc/dwarfs/"+verdir
        
    t0 = time.time()

    print(sys.argv)

    # Not enough inputs
    n = len(sys.argv)
    if n < 4:
        print("Syntax - ns_dwarfs.py ra dec radius version")
        sys.exit()

    # Inputs
    ra0 = np.float(sys.argv[1])
    dec0 = np.float(sys.argv[2])
    radius = np.float(sys.argv[3])
    
    #ra0 = 185.43
    #dec0 = -31.99
    #radius = 1.0

    # Output name
    outbase = str(ra0)+'_'+str(dec0)+'_'+str(radius)+'_'+str(version)

    # Check if the output already exists
    if os.path.exists(dir+outbase+'_peaks.fits'):
        print(dir+outbase+'_peaks.fits already exists')
        sys.exit()
    
    # Set up logging to screen and logfile
    #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    #logfile = tmpdir+"/"+base+".log"
    logfile = dir+outbase+".log"
    #fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, fileName))
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    #rootLogger.setLevel(logging.NOTSET)
    rootLogger.setLevel(logging.INFO)

    rootLogger.info("Searching for overdensities at RA="+str(ra0)+" DEC="+str(dec0)+" Radius="+str(radius)+" on host="+host)
    #rootLogger.info("  Temporary directory is: "+tmpdir)


    # Either get token for anonymous user
    token = ac.login('anonymous')

    # Authenticated users please uncomment the next line
    #token = ac.login(input("Enter user name: "),getpass("Enter password: "))

    #df0 = getData(ra0,dec0,radius=0.01)
    #print(str(len(df0))+' objects found')
    df = getBlueStarData(ra0,dec0,radius=radius,colcutlo=-0.2,colcuthi=0.8,classcut=0.6,fwhmcut=1.5,errcut=0.1)
    rootLogger.info(str(len(df))+' objects found')
    if len(df) == 0:
        rootLogger.info('No data')
        sys.exit()
    
    # Get coverage information
    cov = getCovData(ra0,dec0,radius=radius,columns='ra,dec,pix,pix128,gcoverage,gdepth,rcoverage,rdepth')
    rootLogger.info(str(len(cov))+' coverage pixels returned')

    # Create the healpix map
    #NSIDE = 1024  #4096
    #map = np.zeros(hp.nside2npix(NSIDE),dtype='float')
    #map[:] = hp.UNSEEN   # all unseen/masked to start
    #map = hp.ma(np.zeros(hp.nside2npix(NSIDE),dtype='float'))
    #map.mask = True
    #map[:] = hp.UNSEEN
    #objpix = hp.pixelfunc.ang2pix(NSIDE,df['ra'],df['dec'],lonlat=True)


    # Make healpix density map
    #npix = hp.nside2npix(NSIDE)
    #hist,bin_edges = np.histogram(objpix,bins=npix,range=[0,npix-1])
    #ind, = np.where(hist > 0)
    #map[ind] = hist[ind]
    #map[ind].mask = False

    # Downgrade the coverage map resolution
    #gcov_map2 = hp.pixelfunc.ud_grade(gcov_map,1024)
    #rcov_map2 = hp.pixelfunc.ud_grade(rcov_map,1024)

    #%%time
    ## 20 and 2 arcmin
    #smap1 = hp.sphtfunc.smoothing(map,fwhm=(20./60.)*(3.14159/180.),iter=1,lmax=2*NSIDE)
    #bmap1 = hp.sphtfunc.smoothing(map,fwhm=(120./60.)*(3.14159/180.),iter=1,lmax=2*NSIDE)
    #smmap1 = smap1-bmap1
    ##smmap2 = hp.sphtfunc.smoothing(map2,fwhm=(20./60.)*(3.14159/180.))

    #%%time
    small_k, big_k = 2., 20.  # kernel sizes in arcminutes
    raw, extent, delta, clipped, dsigma = dwarf_filter(df['ra'],df['dec'],fwhm_small=small_k,fwhm_big=big_k)

    #fig, ax = plt.subplots(figsize=(8,8))
    #im = plt.imshow(clipped)
    #plt.xlabel('pixel')
    #plt.ylabel('pixel')
    #plt.colorbar(label='relative spatial density after convolution');
    #plt.savefig(dir+outbase+'_skymap.png')
    
    # find peaks
    small_k = 2.0
    mn, med, std = stats.sigma_clipped_stats(clipped,sigma=3.0,iters=5)
    nsigma = 2.5
    tbl = find_peaks(clipped,med+nsigma,box_size=small_k*2)

    # add ra & dec positions of peaks found
    a, b = extent[:2]
    xvec = np.arange(a,b,(b-a)/clipped.shape[1])
    a, b = extent[2:]
    yvec = np.arange(a,b,(b-a)/clipped.shape[0])

    tbl['ra'] = xvec[tbl['x_peak']]
    tbl['dec'] = yvec[-tbl['y_peak']-1]
    #print(tbl)
    rootLogger.info(str(tbl))

    fig, ax = plt.subplots(figsize=(8,8))
    im = plt.imshow(clipped)
    ax.scatter(tbl['x_peak'],tbl['y_peak'],marker='s',s=tbl['peak_value']*40,c='none',edgecolors='w',lw=3) # keeps writing to previous ax
    plt.xlabel('pixel')
    plt.ylabel('pixel')
    plt.title('%d overdensities, RA=%f  DEC=%f' % (len(tbl), ra0, dec0))
    plt.colorbar(label='relative spatial density after convolution');
    plt.savefig(dir+outbase+'_skymap_peaks.png')

    
    #ecs = ['w','y'] # colors of box frames
    #ax.scatter(tbl['x_peak'],tbl['y_peak'],marker='s',s=tbl['peak_value']*40,c='none',edgecolors=ecs,lw=3) # keeps writing to previous ax
    #ax.scatter(tbl['x_peak'],tbl['y_peak'],marker='s',s=tbl['peak_value']*40,c='none',edgecolors='w',lw=3) # keeps writing to previous ax
    #fig  # repeats (the updated) figure
    #plt.savefig(dir+outbase+'_skymap_peak.png')

    # Loop over the peaks and download data
    for i in range(len(tbl)):
        tbl0 = tbl[i]
        rootLogger.info("  i="+str(i+1)+" RA="+str(tbl0['ra'])+" DEC="+str(tbl0['dec']))
        
        # Convert X/Y to RA/DEC
        df0 = getBlueStarData(tbl0['ra'],tbl0['dec'],radius=0.1,classcut=0.6,fwhmcut=1.5)
        print(str(len(df))+' objects found')

        #fig = plt.figure(figsize=(8,8))
        #plotpanel(121,df0['gmag']-df0['rmag'],df0['gmag'],'Overdensity CMD, %d objects RA=%f  DEC=%f' % (len(df0), tbl0['ra'], tbl0['dec']))
        #plotpanel(122,R1['g_r'],R1['gmag'],'yellow box, %d objects' % len(R1))

        fig =  plt.figure(figsize=(8,8))
        plt.scatter(df0['gmag']-df0['rmag'],df0['gmag'],marker='.',s=10, alpha=0.8)
        plt.xlabel('g-r')
        plt.ylabel('g')
        xlim=(-1,2)
        ylim=(25.2,14)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.title('Overdensity CMD, %d objects RA=%f  DEC=%f' % (len(df0), tbl0['ra'], tbl0['dec']))
        plt.savefig(dir+outbase+'_cmd'+str(i+1)+'.png')
        
        
    # Save the table
    rootLogger.info("Saving info to "+dir+outbase+'_peaks.fits')
    tbl.write(dir+outbase+'_peaks.fits')
