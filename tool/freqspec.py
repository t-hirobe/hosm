#!/usr/local/app/python/2.7.5/bin/python
import struct
import sys
import math
import argparse
import numpy as np

int_size = 4
float_size = 4

def read_int(file):
    return struct.unpack('i', file.read(int_size))[0]

def read_float(file):
    return struct.unpack('f', file.read(float_size))[0]

def read_grid_data_int(file, mx, my):
    data = file.read((mx*my)*int_size)
    return np.array(struct.unpack('{0:d}i'.format(mx*my), data)).reshape((mx, my))

def read_grid_data_float(file, mx, my):
    data = file.read((mx*my)*float_size)
    return np.array(struct.unpack('{0:d}f'.format(mx*my), data)).reshape((mx, my))



#if __name__ == '__main__':


parser = argparse.ArgumentParser(description='')
parser.add_argument('-f', '--file', metavar='file_name',
                    help='specify wave.bin file')
parser.add_argument('-s', '--step', metavar='step', type=int,
                    help='specify step number')
parser.add_argument('-i', '--info', action='store_true',
                    help='')
args  = vars(parser.parse_args())
filename  = args['file']
show_info = args['info']
seek_num  = args['step']


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
outputStep = read_int(file)
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

k_list = ((kx_list*dkx)**2 + (ky_list*dky)**2)**0.5

if show_info:
    print "M         = {0}".format(level)
    print " nx, ny   = {0}, {1}".format(nx, ny)
    print " lx, ly   = {0}, {1}".format(lx, ly)
    print "dky, dky  = {0}, {1}".format(dkx, dky)
    print "-mkx:mkx  = {0} : {1}".format(-mkx, mkx)
    print "-mky:mky  = {0} : {1}".format(-mky, mky)
    print "dt        = {0}".format(dt)
    print "end time  = {0}".format(endTime)
    print "cfl       = {0}".format(cfl)
    exit()

# one time data
one_time_data_len = 4*4*(mx*my)+4
file.seek(one_time_data_len * seek_num, 1)
time = read_float(file)

keta_h = read_grid_data_float(file, mx, my)
keta_p = read_grid_data_float(file, mx, my)
kphi_h = read_grid_data_float(file, mx, my)
kphi_p = read_grid_data_float(file, mx, my)
#keta = keta_h * np.exp(keta_p * 1.0j) * 2.0 * math.pi
#kphi = kphi_h * np.exp(kphi_p * 1.0j) * 2.0 * math.pi
keta = keta_h * np.exp(keta_p * 1.0j)
kphi = kphi_h * np.exp(kphi_p * 1.0j)

grav = 1.0
w_list = (grav*k_list)**0.5
k_list[mkx][mky] = 1
w_list[mkx][mky] = 1
b_list = (w_list/2.0/k_list)**0.5 * keta + 1j*(k_list/2.0/w_list)**0.5 * kphi
b_list[mkx][mky] = 0

dw = 1/20.0
vals = w_list * (abs(b_list)**2)
max_w = np.max(w_list)
output = [[] for i in range(int(max_w/dw + 2))]
for ix in range(0, mx):
    for iy in range(0, my):
        w = w_list[ix][iy]
        w_id = int(w/dw+0.5)
        output[w_id].append( vals[ix][iy] / dw )
        pass
    pass

print "# Time = {0}".format(time)
for iw in range(0, int(max_w/dw+1)):
    print iw*dw, np.sum(np.array(output[iw])), len(output[iw])
    pass
