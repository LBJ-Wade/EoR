#!/usr/bin/python
# Removal of handle stuff as async is defunct.
def fov_to_cellsize(fov_deg, size):
    import numpy as np
    rmax = np.sin(fov_deg/2.0*(np.pi/180.0))
    inc = rmax / (0.5 * size)
    return np.arcsin(inc)*((180.0*3600.0)/np.pi)

def casa_make_images(msname, fov_deg, size, rootname):
    from USER_DEFINED import runtime_variables as rv

    default(clean)
    cellsize = fov_to_cellsize(fov_deg, size)
    print 'cellsize = ',cellsize
    handle = clean(
        vis=msname,
        imagename=rootname,
        niter=rv.user_image['niter'], #cw 0,
        gridmode=rv.user_image['gridmode'], #cw 'widefield',
        wprojplanes=rv.user_image['wprojplanes'],#cw 256,
        imsize=[size,size],
        cell=[cellsize,cellsize],
	    uvrange=rv.user_image['uvrange'], #cw '0.01~0.8klambda', # added Emma 15/12
	    weighting=rv.user_image['weighting'] #cw 'uniform'
        )
    print niter
    print gridmode
    print rv.user_image['wprojplanes']
    print uvrange
    print weighting

    return handle

#===============================================================================

if __name__ == '__main__':
    import sys
    import time
    import os
    import shutil
    from USER_DEFINED import runtime_variables as rv

    # This script should be run with:
    # $> casapy --nogui --nologger --log2term -c make_images_simple.py [start channel] [end channel]

    if len(sys.argv) < 3:
        print len(sys.argv)
        print 'Usage: '
        print '  casapy --nogui --nologger --log2term -c make_images_simple.py [start channel] [end channel]'
        print ''
        sys.exit(2)

    #===========================================================================
    size      = rv.user_image['num_pixels_side'] #cw  image dimensions
    fov       = rv.user_image['fov_deg'] #cw deg.
    #===========================================================================

    visDir    = '../Vis'
    #cw Nb. visDir is a variable in Emma's, maybe it should become another global
    #cw (if do this then make sure to replace /Vis in runOSKAR.py)

    end_channel   = int(sys.argv[-1])
    start_channel = int(sys.argv[-2])

    # List of measurement sets in the vis folder (within channel range)
    listdir_ = os.listdir(visDir)
    listdir_.sort()
    msList = []
    for i in range(0, len(listdir_)):
        path_ = os.path.join(visDir, listdir_[i])
        if os.path.isdir(path_) and path_[-2:] == 'ms':
            channelID = int(path_[-6:-3])
            if channelID >= start_channel and channelID <= end_channel:
                #if listdir_[i][0] == 'n':
            	msList.append(listdir_[i])
        	print '%s/%s' % (visDir, msList[-1][0:])

    for i in range(0, len(msList)):
        ms  = '%s/%s' % (visDir, msList[i])
        img = '%s_%04i_%.1f' % (msList[i][0:-3], size, fov)

        print ''
        print '************************************************************'
        print '*  - ms     = %s' % (ms)
        print '*  - images = %s.{image|psf|...}' % (img)
        print '*  - FoV    = %f deg.' % (fov)
        print '*  - size   = %i' % (size)
        print '************************************************************'
        print ''

        casa_make_images(ms, fov, size, img)

    # Convert CASA images to fits.
    for f in os.listdir('.'):
	print "fname: %s"%f
        if (f.endswith('.image') and not os.path.isfile(f+'.fits')):
            exportfits(imagename=f, fitsimage=f+'.fits')
        if (f.endswith('.psf') and not os.path.isfile(f+'.fits')):
            exportfits(imagename=f, fitsimage=f+'.fits')
        if (not f.endswith('.fits')):
            print '>> removing CASA image: %s' % (f)
            if os.path.isdir(f):
                        shutil.rmtree(f)
            else:
                  	os.remove(f)
