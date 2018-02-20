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
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy import arange, array, exp
#%matplotlib inline
import healpy as hp
from photutils import find_peaks, data_properties
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


# Extrapolation function
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

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
        
    df = helpers.convert(result,'table')   # pandas
    
    return df

# A function to retrieve data from a point on the sky
def getDataCuts (ra,dec,radius=1.0,columns='ra,dec,gmag,rmag',colcutlo=None,colcuthi=None,classcut=None,fwhmcut=None,errcut=None):

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
        
    df = helpers.convert(result,'table')  # pandas
    
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
    if len(sys.argv) > 2:
       version = sys.argv[3]
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
    if n < 2:
        print("Syntax - ns_dwarfs_hpix.py hpix radius version")
        sys.exit()

    # Inputs
    hpix = np.int(sys.argv[1])
    radius = np.float(sys.argv[2])

    # Convert healpix to ra/dec
    nside = np.int(64)
    ra0, dec0 = hp.pix2ang(nside,hpix,lonlat=True)

    # Output name
    outbase = str(hpix)+'_'+str(radius)+'_'+str(version)
    outdir = dir+str(hpix/1000)+'/'+str(hpix)+'/'
    if not os.path.exists(outdir):     # make output directory if necessary
        os.makedirs(outdir)
    outfile = outdir+outbase+'_peaks.fits'
    donefile = outdir+outbase+'.done'

    # Check if the "done" file already exists
    if os.path.exists(donefile):
        print('This healpix was already done')
        sys.exit()

    # Set up logging to screen and logfile
    #logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    #logfile = tmpdir+"/"+base+".log"
    logfile = outdir+outbase+".log"
    # Delete file if it exists
    if os.path.exists(logfile):
        os.remove(logfile)
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
    bcatall = getDataCuts(ra0,dec0,radius=radius,colcutlo=-0.2,colcuthi=0.8,classcut=0.6,fwhmcut=1.5,errcut=0.1)
    rootLogger.info(str(len(bcatall))+' objects found')
    if len(bcatall) == 0:
        rootLogger.info('No data')
        # Create done file
        f = open(donefile,'w')
        f.write(host)
        f.close()
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

    # Run dwarf filter
    small_k, big_k = 2., 20.  # kernel sizes in arcminutes
    try:
        raw, extent, delta, clipped, dsigma = dwarf_filter(bcatall['ra'],bcatall['dec'],fwhm_small=small_k,fwhm_big=big_k)
    except:
        rootLogger.info('Problems with dwarf filter')
        # Create done file
        f = open(donefile,'w')
        f.write(host)
        f.close()
        sys.exit()

    # find peaks
    small_k = 2.0
    mn, med, std = stats.sigma_clipped_stats(clipped,sigma=3.0,iters=5)
    nsigma = 2.5
    tbl = find_peaks(clipped,med+nsigma,box_size=small_k*2)
    rootLogger.info(str(len(tbl))+' peaks found')
    if len(tbl) == 0:
        # Create done file
        f = open(donefile,'w')
        f.write(host)
        f.close()
        sys.exit()

    # add ra & dec positions of peaks found
    a, b = extent[:2]
    xvec = np.arange(a,b,(b-a)/clipped.shape[1])
    a, b = extent[2:]
    yvec = np.arange(a,b,(b-a)/clipped.shape[0])

    tbl['ra'] = xvec[tbl['x_peak']]
    tbl['dec'] = yvec[-tbl['y_peak']-1]
    #print(tbl)
    rootLogger.info(str(tbl))

    # Make sky map of overdensities
    fig, ax = plt.subplots(figsize=(8,8))
    im = plt.imshow(clipped)
    ax.scatter(tbl['x_peak'],tbl['y_peak'],marker='s',s=tbl['peak_value']*40,c='none',edgecolors='w',lw=3) # keeps writing to previous ax
    plt.xlabel('pixel')
    plt.ylabel('pixel')
    plt.title('%d overdensities, RA=%f  DEC=%f' % (len(tbl), ra0, dec0))
    plt.colorbar(label='relative spatial density after convolution');
    plt.savefig(outdir+outbase+'_skymap_peaks.png')

    # Create output table
    dt = np.dtype([('id',np.str_,20),('hpix',int),('num',int),('x_peak',int),('y_peak',int),('peak_value',float),('ra_peak',float),('dec_peak',float),
                   ('x_centroid',float),('y_centroid',float),('ra',float),('dec',float),('asemi',float),('bsemi',float),('theta',float),
                   ('back_color',float),('back_mag',float),('back_fwhm',float),('back_class_star',float),('back_gdepth',float),('back_rdepth',float),
                   ('nobj',int),('color',float),('mag',float),('fwhm',float),('class_star',float),('nblue',int),('blue_mag',float),
                   ('gdepth',float),('rdepth',float),
                   ('pmra',float),('pmdec',float),('pmraerr',float),('pmdecerr',float),('nexp',int),('deltamjd',float)])
    peaks = None
    npeaks = len(tbl)
    if npeaks > 0:
        peaks = np.zeros(npeaks,dtype=dt)
        peaks['hpix'] = hpix
        peaks['num'] = np.arange(npeaks)+1
        tempid = np.char.add(peaks['hpix'].astype(str),'.')
        peaks['id'] = np.char.add(tempid,peaks['num'].astype(str))
        # Copying info in tbl
        peaks['x_peak'] = tbl['x_peak']
        peaks['y_peak'] = tbl['y_peak']
        peaks['peak_value'] = tbl['peak_value']
        peaks['ra_peak'] = tbl['ra']
        peaks['dec_peak'] = tbl['dec']
        # Add the "background" values
        columns = 'ra,dec,gmag,rmag,fwhm,class_star'
        catall = getDataCuts(ra0,dec0,radius=radius,columns=columns)
        rootLogger.info(str(len(catall))+' total objects found')
        gdcol = (catall['gmag'] < 50) & (catall['rmag'] < 50)
        if np.sum(gdcol) > 0:
            peaks['back_color'] = np.median(catall['gmag'][gdcol]-catall['rmag'][gdcol])
        else:
            peaks['back_color'] = 999999.
        gdmag = (catall['gmag'] < 50)
        if np.sum(gdmag) > 0:
            peaks['back_mag'] = np.median(catall['gmag'][gdmag])
        else:
            peaks['back_mag'] = 999999.
        peaks['back_fwhm'] = np.median(catall['fwhm'])
        peaks['back_class_star'] = np.median(catall['class_star'])
        # Depth information
        gdgdepth = (cov['gdepth'] < 50) & (cov['gdepth'] > 0)
        if np.sum(gdgdepth) > 0:
            peaks['back_gdepth'] = np.median(cov['gdepth'][gdgdepth])
        else:
            peaks['back_gdepth'] = 999999.
        gdrdepth = (cov['rdepth'] < 50) & (cov['rdepth'] > 0)
        if np.sum(gdrdepth) > 0:
            peaks['back_rdepth'] = np.median(cov['rdepth'][gdrdepth])
        else:
            peaks['back_rdepth'] = 999999.
        

    # Loop over the peaks and download data
    for i in range(npeaks):
        peaks0 = peaks[i]
        rootLogger.info("  i="+str(i+1)+" RA="+str(peaks0['ra_peak'])+" DEC="+str(peaks0['dec_peak']))
        
        # Convert X/Y to RA/DEC
        #cat0 = getDataCuts(peaks0['ra'],peaks0['dec'],radius=0.1,classcut=0.6,fwhmcut=1.5)
        columns = 'ra,dec,pmra,pmdec,pmraerr,pmdecerr,ndet,deltamjd,gmag,rmag,fwhm,class_star'
        cat0 = getDataCuts(peaks0['ra_peak'],peaks0['dec_peak'],radius=0.1,columns=columns)
        print(str(len(cat0))+' objects found')

        # Measure the morphology around the overdensity
        shp = clipped.shape
        x0 = np.int(np.floor(peaks0['x_peak']-10))
        if x0 < 0: x0=0
        x1 = np.int(np.ceil(peaks0['x_peak']+10))
        if x1 > (shp[1]-1): x1=(shp[1]-1)   # X is 2nd dimension
        y0 = np.int(np.floor(peaks0['y_peak']-10))
        if y0 < 0: y0=0
        y1 = np.int(np.ceil(peaks0['y_peak']+10))
        if y1 > (shp[0]-1): y1=(shp[0]-1)   # Y is 1st dimension
        clipped0 = clipped[y0:y1+1,x0:x1+1]
        props = data_properties(clipped0)
        pcolumns = ['id', 'xcentroid', 'ycentroid', 'semimajor_axis_sigma','semiminor_axis_sigma', 'orientation']
        #  semi axes in pixels and orientation in radians
        # 1 pixel is 1 armcin, good unit to use
        pcat = props.to_table(columns=pcolumns)
        peaks[i]['x_centroid'] = props['xcentroid'].value+x0
        peaks[i]['y_centroid'] = props['ycentroid'].value+y0
        peaks[i]['asemi'] = props['semimajor_axis_sigma'].value  # pixel=arcmin
        peaks[i]['bsemi'] = props['semiminor_axis_sigma'].value  # pixel=arcmin
        peaks[i]['theta'] = np.rad2deg(props['orientation'].value)
        # ra & dec positions of centroid
        #   use xvec/yvec defined above
        xf = extrap1d(interp1d(np.arange(len(xvec)),xvec))  # function to interpolate x
        peaks[i]['ra'] = xf([peaks[i]['x_centroid']])
        yf = extrap1d(interp1d(np.arange(len(yvec)),yvec))  # function to interpolate y
        peaks[i]['dec'] = np.float( yf([len(yvec)-peaks[i]['x_centroid']-1]) )
        #peaks[i]['ra'] = xvec[peaks[i]['x_centroid']]
        #peaks[i]['dec'] = yvec[-peaks[i]['y_centroid']-1]

        # Add the median values for this peak
        gdcol = (cat0['gmag'] < 50) & (cat0['rmag'] < 50)
        if np.sum(gdcol) > 0:
            peaks[i]['color'] = np.median(cat0['gmag'][gdcol]-cat0['rmag'][gdcol])
        else:
            peaks[i]['color'] = 999999.
        gdmag = (cat0['gmag'] < 50)
        if np.sum(gdmag) > 0:
            peaks[i]['mag'] = np.median(cat0['gmag'][gdmag])
        else:
            peaks[i]['mag'] = 999999.
        peaks[i]['fwhm'] = np.median(cat0['fwhm'])
        peaks[i]['class_star'] = np.median(cat0['class_star'])
        peaks[i]['nobj'] = len(cat0)
        gdpm = ((np.abs(cat0['pmra']) < 1e5) & np.isfinite(cat0['pmra']) &
                (np.abs(cat0['pmdec']) < 1e5) & np.isfinite(cat0['pmdec']))
        if np.sum(gdpm) > 0:
            peaks[i]['pmra'] = np.median(cat0['pmra'][gdpm])
            peaks[i]['pmdec'] = np.median(cat0['pmdec'][gdpm])
            peaks[i]['pmraerr'] = np.median(cat0['pmraerr'][gdpm]) / np.sqrt(np.sum(gdpm))
            peaks[i]['pmdecerr'] = np.median(cat0['pmdecerr'][gdpm]) / np.sqrt(np.sum(gdpm))
            peaks[i]['nexp'] = np.median(cat0['ndet'][gdpm])
            peaks[i]['deltamjd'] = np.median(cat0['deltamjd'][gdpm])
        else:
            peaks[i]['pmra'] = 999999.
            peaks[i]['pmdec'] = 999999.
            peaks[i]['pmraerr'] = 999999.
            peaks[i]['pmdecerr'] = 999999.
            peaks[i]['nexp'] = np.median(cat0['ndet'])
            peaks[i]['deltamjd'] = np.median(cat0['deltamjd'])
        # Select only very blue stars
        gdblue = ((cat0['gmag']-cat0['rmag']) < 0.1) & (cat0['gmag'] < 50) & (cat0['rmag'] < 50)
        peaks[i]['nblue'] = np.sum(gdblue)
        if np.sum(gdblue) > 0:
            peaks[i]['blue_mag'] = np.median(cat0['gmag'][gdblue])
        else:
            peaks[i]['blue_mag'] = 999999.
        # Get coverage information
        cov0 = getCovData(peaks0['ra'],peaks0['dec'],radius=0.1,columns='ra,dec,pix,pix128,gcoverage,gdepth,rcoverage,rdepth')
        gdgdepth = (cov0['gdepth'] < 50) & (cov0['gdepth'] > 0)
        if np.sum(gdgdepth) > 0:
            peaks[i]['gdepth'] = np.median(cov0['gdepth'][gdgdepth])
        else:
            peaks[i]['gdepth'] = 999999.
        gdrdepth = (cov0['rdepth'] < 50) & (cov0['rdepth'] > 0)
        if np.sum(gdrdepth) > 0:
            peaks[i]['rdepth'] = np.median(cov0['rdepth'][gdrdepth])
        else:
            peaks[i]['rdepth'] = 999999.

        # Make sky map with THIS overdensity highlighted
        fig, ax = plt.subplots(figsize=(8,8))
        im = plt.imshow(clipped)
        ax.scatter(peaks['x_peak'],peaks['y_peak'],marker='s',s=peaks['peak_value']*40,c='none',edgecolors='w',lw=3) # keeps writing to previous ax
        ax.scatter(peaks0['x_peak'],peaks0['y_peak'],marker='s',s=peaks0['peak_value']*40,c='none',edgecolors='y',lw=3)
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('%d overdensities, RA=%f  DEC=%f' % (npeaks, ra0, dec0))
        plt.colorbar(label='relative spatial density after convolution');
        plt.savefig(outdir+outbase+'_skymap_peak'+str(i+1)+'.png')

        # Make CMD for this overdensity
        fig =  plt.figure(figsize=(8,8))
        plt.scatter(cat0['gmag']-cat0['rmag'],cat0['gmag'],marker='.',s=10, alpha=0.8)
        plt.xlabel('g-r')
        plt.ylabel('g')
        xlim=(-1,2)
        ylim=(25.2,14)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.title('Overdensity CMD, %d objects RA=%f  DEC=%f' % (len(cat0), peaks0['ra'], peaks0['dec']))
        plt.savefig(outdir+outbase+'_cmd'+str(i+1)+'.png')

        # Make CMD for this overdensity ONLY GALAXIES
        gals = cat0['class_star'] < 0.1
        fig =  plt.figure(figsize=(8,8))
        plt.scatter(cat0[gals]['gmag']-cat0[gals]['rmag'],cat0[gals]['gmag'],marker='.',s=10, alpha=0.8)
        plt.xlabel('g-r')
        plt.ylabel('g')
        xlim=(-1,2)
        ylim=(25.2,14)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.title('Overdensity CMD Galaxies, %d objects RA=%f  DEC=%f' % (len(cat0), peaks0['ra'], peaks0['dec']))
        plt.savefig(outdir+outbase+'_cmd'+str(i+1)+'_gals.png')

        plt.close('all')  # close all figures
        
    # Save the table
    rootLogger.info("Saving info to "+outdir+outbase+'_peaks.fits')
    #peaks.write(outdir+outbase+'_peaks.fits')
    Table(peaks).write(outfile)

    # Create done file
    f = open(donefile,'w')
    f.write(host)
    f.close()
