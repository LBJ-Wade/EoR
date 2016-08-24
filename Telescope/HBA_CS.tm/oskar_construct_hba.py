#!/usr/bin/python

from lofar_meta_data_utils import *;
from math import *;
import matplotlib.pyplot as plt;
import matplotlib;
import numpy;
import os;

telescope_name = 'LOFAR_HBA_CS_RS';

filenames = [
'CS001',
'CS002',
'CS003',
'CS004',
'CS005',
'CS006',
'CS007',
'CS011',
'CS013',
'CS017',
'CS021',
'CS024',
'CS026',
'CS028',
'CS030',
'CS031',
'CS032',
'CS101',
'CS103',
'CS201',
'CS301',
'CS302',
'CS401',
'CS501',
'RS106',
'RS205',
'RS208',
'RS305',
'RS306',
'RS310',
'RS307',
'RS406',
'RS407',
'RS409',
'RS503',
'RS508',
'RS509'
]

def make_tile_data(dirname):
    # Load the layout.txt file.
    coords = numpy.loadtxt("layout.txt");
    
    # Get the x and y coordinates of two reference positions.
    ref_x = coords[0,0];
    ref_y = coords[0,1];
    first_x = coords[1,0];
    first_y = coords[1,1];
    
    # Get the angle between the first two tile centres.
    dir_x = first_x - ref_x;
    dir_y = first_y - ref_y;
    angle = atan2(dir_y, dir_x);
    angles.write("%f\n" % angle);
    
    # Transform the default element coordinates to match the tile orientation.
    rot_mat = numpy.array([[cos(angle),-sin(angle)],[sin(angle),cos(angle)]]);
    num_elements = tile_dim * tile_dim;
    coord = numpy.zeros([num_elements, 2]); # x and y element coordinates
    for i in range(0, num_elements):
        element_coord = numpy.array([[element_pos_x[i]], [element_pos_y[i]]]);
        coord[i] = numpy.transpose(numpy.dot(rot_mat, element_coord));
        
    # Write out the coordinates to the tile layout.
    dirname = "tile";
    os.mkdir(dirname);
    os.chdir(dirname);
    numpy.savetxt("layout.txt", coord);
    os.chdir('..')
    
    # Plot (debug).
    #matplotlib.axes.Axes.set_aspect(plt.gca(), 1)
    #plt.plot(coords[:,0], coords[:,1], 'ro')
    #plt.plot(coord[:,0], coord[:,1], 'bo')
    #plt.show()


def rotate_data(data, matrix, ref_pos):
    return numpy.dot((data - ref_pos), matrix);

# Load the station data as a multi-level dictionary.
num_stations = len(filenames);
station_dict = load_data();

# Make telescope directory, and enter it.
os.mkdir(telescope_name);
os.chdir(telescope_name);

# Write the top-level layout file.
fp = open("layout_ecef.txt", "w");
for i in range(0, num_stations):
    f = filenames[i];
    st = station_dict[f];
    #print f
    if f.startswith("CS"):
        fp.write("%f %f %f\n" % (st['hba0_c'][0], st['hba0_c'][1], st['hba0_c'][2]));
        fp.write("%f %f %f\n" % (st['hba1_c'][0], st['hba1_c'][1], st['hba1_c'][2]));
        #print st['hba0_c']
        #print st['hba1_c']
    else:
        fp.write("%f %f %f\n" % (st['hba_c'][0], st['hba_c'][1], st['hba_c'][2]));
        #print st['hba_c']
fp.close();

# Generate default 4x4 tile, 1.25 m separation.
tile_dim = 4;
t = numpy.linspace(-1.875, 1.875, tile_dim);
[element_pos_x, element_pos_y] = numpy.meshgrid(t, t);
element_pos_x = numpy.reshape(element_pos_x, -1, 1);
element_pos_y = numpy.reshape(element_pos_y, -1, 1);

# Create file for station orientation angles.
angles = open("angles.txt", "w");

# Loop over stations.
for i in range(0, num_stations):
    f = filenames[i];
    st = station_dict[f];
    
    # Extract both HBA0 and HBA1 as separate stations from CS* stations.
    if f.startswith("CS"):
        # Create directory for first HBA sub-station.
        dirname = "{:s}_HBA0".format(f);
        os.mkdir(dirname);
        os.chdir(dirname);
        
        # Write tile coordinates for HBA0.
        dat = rotate_data(st['hba0'], st['hba0_c_rot'], st['hba0_c']).reshape((-1, 3));
        # Get every other row (since identical positions are given for each dipole).
        dat = dat[::2, :];
        numpy.savetxt("layout.txt", dat);
        
        # Make tile directory and data.
        make_tile_data(dirname);
        os.chdir('..')

        # Create directory for second HBA sub-station.
        dirname = "{:s}_HBA1".format(f);
        os.mkdir(dirname);
        os.chdir(dirname);
        
        # Write tile coordinates for HBA1.
        dat = rotate_data(st['hba1'], st['hba1_c_rot'], st['hba1_c']).reshape((-1, 3));
        # Get every other row (since identical positions are given for each dipole).
        dat = dat[::2, :];
        numpy.savetxt("layout.txt", dat);
        
        # Make tile directory and data.
        make_tile_data(dirname);
        os.chdir('..')
    else:
        # Create directory for HBA station.
        dirname = "{:s}_HBA".format(f);
        os.mkdir(dirname);
        os.chdir(dirname);
        
        # Write tile coordinates for HBA.
        dat = rotate_data(st['hba'], st['hba_c_rot'], st['hba_c']).reshape((-1, 3));
        # Get every other row (since identical positions are given for each dipole).
        dat = dat[::2, :];
        numpy.savetxt("layout.txt", dat);
        
        # Make tile directory and data.
        make_tile_data(dirname);
        os.chdir('..')

angles.close();

