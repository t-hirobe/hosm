#!/usr/local/app/python/2.7.5/bin/python
import struct
import sys
import math
import argparse

int_size = 4
float_size = 4

def read_int(file):
    return struct.unpack('i', file.read(int_size))[0]

def read_float(file):
    return struct.unpack('f', file.read(float_size))[0]

def read_grid_data_int(file, mx, my):
    data = file.read((mx*my)*int_size)
    return np.array(struct.unpack('{0:d}i'.format(mx*my), data)).reshape((mx, my))

parser = argparse.ArgumentParser(description='')
parser.add_argument('-f', '--file', metavar='file_name',
                    help='specify wave.bin file')
parser.add_argument('-x', '--kx', metavar='kx', type=int,
                    help='')
parser.add_argument('-y', '--ky', metavar='ky', type=int,
                    help='')
parser.add_argument('-i', '--info', action='store_true',
                    help='')
args  = vars(parser.parse_args())
filename  = args['file']
show_info = args['info']
kx        = args['kx']
ky        = args['ky']


file = open(filename, "rb")

# read header
level = read_int(file)
nx    = read_int(file)
ny    = read_int(file)
lx    = read_float(file)
ly    = read_float(file)
dkx   = read_float(file)
dky   = read_float(file)
dt    = read_float(file)
endTime  = read_float(file)
outputStep  = read_int(file)
grav = read_float(file)

mkx = int(nx/(level+1))
mky = int(ny/(level+1))
mx = int(2*mkx+1)
my = int(2*mky+1)

#Courant condition
cmax = math.sqrt(grav/max(dkx, dky))
dx = lx / nx
cfl = dt * cmax / dx


kx_list = read_grid_data_int(file, mx, my)
ky_list = read_grid_data_int(file, mx, my)

if show_info:
    print "M         = {0}".format(level)
    print " nx, ny   = {0}, {1}".format(nx, ny)
    print " ly, ly   = {0}, {1}".format(lx, ly)
    print "dky, dky  = {0}, {1}".format(dkx, dky)
    print "-mkx:mkx  = {0} : {1}".format(-mkx, mkx)
    print "-mky:mky  = {0} : {1}".format(-mky, mky)
    print "dt        = {0}".format(dt)
    print "end time  = {0}".format(endTime)
    print "cfl       = {0}".format(cfl)
    exit()

# one time data
#one_time_data_len = 4*4*(mx*my)+4

pos = (kx + mkx)*my + ky + mky
time = 0.0
print "# kx, ky = {0}, {1}".format(kx_list[pos], ky_list[pos])
while time <= endTime:
    try:
        time = struct.unpack('f', file.read(float_size))[0]
        file.seek(4*pos, 1)
        eta_h = struct.unpack('f', file.read(float_size))[0]
        file.seek(4*(mx*my - pos - 1), 1)
        file.seek(4*pos, 1)
        eta_p = struct.unpack('f', file.read(float_size))[0]
        file.seek(4*(mx*my - pos - 1), 1)
        file.seek(4*pos, 1)
        phi_h = struct.unpack('f', file.read(float_size))[0]
        file.seek(4*(mx*my - pos - 1), 1)
        file.seek(4*pos, 1)
        phi_p = struct.unpack('f', file.read(float_size))[0]
        file.seek(4*(mx*my - pos - 1), 1)

        print "{time:.8f} {eta_h: .5e} {eta_p: .5f} {phi_h: .5e} {phi_p: .5f}".format(
            time = time,
            eta_h = eta_h,
            eta_p = eta_p,
            phi_h = phi_h,
            phi_p = phi_p
            )
        pass
    except IOError:
        break
    except struct.error:
        break
    pass
