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

def rearrange_spectrumgrid(data, mkx, mky, mx, my):
    tmp = np.array([0+0j]*mx*my).reshape((mx, my))
    for ix in range(0, mx):
        for iy in range(0, my):
            if mkx < ix  and mky < iy : tmp[ix][iy] = data[ix-mkx-1][iy-mky-1]
            if mkx < ix  and iy <= mky: tmp[ix][iy] = data[ix-mkx-1][iy+mky]
            if ix <= mkx and mky < iy : tmp[ix][iy] = data[ix+mkx][iy-mky-1]
            if ix <= mkx and iy <= mky: tmp[ix][iy] = data[ix+mkx][iy+mky]
            pass
        pass
    return tmp


def show_curve(kx_list, ky_list, keta, lx, ly, mx, my):
    mkx = (mx-1)/2
    mky = (my-1)/2
    keta_kdx = rearrange_spectrumgrid(keta * kx_list**2 * -1, mkx, mky, mx, my)
    keta_kdy = rearrange_spectrumgrid(keta * ky_list**2 * -1, mkx, mky, mx, my)
    eta_dx = np.fft.ifft2(keta_kdx)*mx*my
    eta_dy = np.fft.ifft2(keta_kdy)*mx*my

    print "# x y ddx ddy"
    for ix in range(0, mx):
        for iy in range(0, my):
            print "{x:.5f} {y:.5f} {dx:.5e} {dy:.5e}".format(
                x = ix * lx/mx,
                y = iy * ly/my,
                dx = eta_dx[ix][iy].real,
                dy = eta_dy[ix][iy].real
                )
            pass
        print ""
        pass
    pass

def show_gradient(kx_list, ky_list, keta, lx, ly, mx, my):
    mkx = (mx-1)/2
    mky = (my-1)/2
    keta_kdx = rearrange_spectrumgrid(keta * kx_list * 1j, mkx, mky, mx, my)
    keta_kdy = rearrange_spectrumgrid(keta * ky_list * 1j, mkx, mky, mx, my)
    eta_dx = np.fft.ifft2(keta_kdx)*mx*my
    eta_dy = np.fft.ifft2(keta_kdy)*mx*my

    print "# x y dx dy"
    for ix in range(0, mx):
        for iy in range(0, my):
            print "{x:.5f} {y:.5f} {dx:.5e} {dy:.5e}".format(
                x = ix * lx/mx,
                y = iy * ly/my,
                dx = eta_dx[ix][iy].real,
                dy = eta_dy[ix][iy].real
                )
            pass
        print ""
        pass
    pass


def show_space(keta, kphi, lx, ly, mx, my):
    mkx = (mx-1)/2
    mky = (my-1)/2
    keta_tmp = rearrange_spectrumgrid(keta, mkx, mky, mx, my)
    kphi_tmp = rearrange_spectrumgrid(kphi, mkx, mky, mx, my)
    eta = np.fft.ifft2(keta_tmp)*mx*my
    phi = np.fft.ifft2(kphi_tmp)*mx*my

    print "# x y eta phi"
    for ix in range(0, mx):
        for iy in range(0, my):
            print "{x:.5f} {y:.5f} {eta:.5e} {phi:.5e}".format(
                x = ix * lx/mx,
                y = iy * ly/my,
                eta = eta[ix][iy].real,
                phi = phi[ix][iy].real
                )
            pass
        print ""
        pass
    pass

def show_space2(kx_list, ky_list, keta, kphi, lx, ly, mx, my):
    print "# x y eta phi"
    for ix in range(0, mx):
        sys.stderr.write("{0}\n".format(ix))
        for iy in range(0, my):
            x = ix * lx/mx
            y = iy * ly/my
            h = np.sum(keta * np.exp(1j * (kx_list*x*2*math.pi/lx + ky_list*y*2*math.pi/ly))).real
            print "{x:.5f} {y:.5f} {eta:.5e}".format(
                x = x,
                y = y,
                eta = h
                )
            pass
        print ""
        sys.stdout.flush()
        pass
    pass

def show_spectrum(kx_list, ky_list, keta, kphi, mx, my):
    print "# kx ky abs(eta) arg(eta) abs(phi) arg(phi)"
    for ix in range(0, mx):
        for iy in range(0, my):
            print "{kx:5d} {ky:5d} {eta_h: .5e} {eta_p: .5f} {phi_h: .5e} {phi_p: .5f}".format(
                kx = kx_list[ix][iy],
                ky = ky_list[ix][iy],
                eta_h = np.abs(keta[ix][iy]),
                eta_p = np.angle(keta[ix][iy]),
                phi_h = np.abs(kphi[ix][iy]),
                phi_p = np.angle(kphi[ix][iy])
                )
            pass
        print ""
        pass
    pass


#if __name__ == '__main__':


parser = argparse.ArgumentParser(description='')
parser.add_argument('-f', '--file', metavar='file_name',
                    help='specify wave.bin file')
parser.add_argument('-s', '--step', metavar='step', type=int,
                    help='specify step number')
parser.add_argument('-i', '--info', action='store_true',
                    help='')
parser.add_argument('-d', '--display', choices='hsgc',
                    help='display option, height, spectrum, gradient or curveture')
args  = vars(parser.parse_args())
filename  = args['file']
show_info = args['info']
seek_num  = args['step']
display  = args['display']


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
keta = keta_h * np.exp(keta_p * 1.0j)
kphi = kphi_h * np.exp(kphi_p * 1.0j)

print "# Time = {0}".format(time)

#if   display == "h":  show_space2(kx_list, ky_list, keta, kphi, lx, ly, mx, my)
if   display == "h":  show_space(keta, kphi, lx, ly, mx, my)
elif display == "s":  show_spectrum(kx_list, ky_list, keta, kphi, mx, my)
elif display == "g":  show_gradient(kx_list, ky_list, keta, lx, ly, mx, my)
elif display == "c":  show_curve(kx_list, ky_list, keta, lx, ly, mx, my)
