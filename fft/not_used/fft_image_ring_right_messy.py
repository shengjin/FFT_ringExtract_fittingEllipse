#!/usr/bin/env python2.7
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos

debug = False
#debug = False
# calculate the Brightness Temperature
#calc_TB = False
calc_TB = True


######### THE MAIN PROGRAM ############

#### THE Spectral Index got from ALMA
spectral_index = 2.72256367113

## 
v_mod = cc*10/1.0
v_obs = 343500000000.0
#av4    343500000000.0
#0	337244575000.0
#1	339182075000.0
#2	347744575000.0
#3	349744575000.0

u_wav = 870.0
#av4 870.0
#0 888.947
#1 883.869
#2 862.105
#3 857.175
v_wav = u_wav


print v_mod, v_obs

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage()
dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out
ImageJyppix_scaled = image.imageJyppix[:,:,n_image]/dpc/dpc
print np.sum(ImageJyppix_scaled, axis=None)*1000.0

# scale the flux using Fv2/Fv1 = (v2/v1)**spectral_index
ImageJyppix_scaled = ImageJyppix_scaled*(v_obs/v_mod)**spectral_index 
print np.sum(ImageJyppix_scaled, axis=None)*1000.0
ImageJyTotal = np.sum(ImageJyppix_scaled, axis=None)

m,n = ImageJyppix_scaled.shape
print "dimension check (m,n):", m, n, "\n"


###################### ###################### 
####### increase the flux in the center
###################### ###################### 

### find the max of the scaled Flux
centerflux = ImageJyppix_scaled.max()
print "centerflux:"
print centerflux, "\n"

# parameters for the center ellipse
cter_au = 2.0 # in AU
# assuming sizepix_x = sizepix_y
e_a_grid = cter_au*au/image.sizepix_x  # long semi-axis
incli = math.radians(46.72)
e_b_grid = e_a_grid * cos(incli)     # short semi-axis
print "e_a_grid, e_b_grid"
print e_a_grid, e_b_grid, "\n"

e_h = 0
e_k = 0
rotang = math.radians(90-138.02)

#####################
# set the grids in pixel
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
fft_dx = pix_x[1] - pix_x[0]
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)
fft_dy = pix_y[1] - pix_y[0]
pix_x_shifted = np.zeros(NN+1)
pix_x_shifted[1:NN+1:1] = pix_x
pix_x_shifted[0] = pix_x[0] - fft_dx
pix_y_shifted = pix_x_shifted

print "fft_dx,fft_dy: "
print fft_dx,fft_dy, "\n"

def outer_boundary(number):
    if (number < 0):
        number = number + 1.5
    else:
        number = number - 0.5
    return number
# for a circle ring
out_au = 50
int_au = 40
out_au_grid = out_au*au/image.sizepix_x  # long semi-axis
int_au_grid = int_au*au/image.sizepix_x  # long semi-axis
print out_au_grid
print int_au_grid

for i in range(m):
    for j in range(n):
        i_cter = outer_boundary(pix_x[i])
        j_cter = outer_boundary(pix_y[j])
        if (i_cter**2+j_cter**2>=int_au_grid**2) and (i_cter**2+j_cter**2<=out_au_grid**2):
            ImageJyppix_scaled[i,j] = 0.1
        else:
            ImageJyppix_scaled[i,j] = 0

ImageJyppix_scaled_shifted = np.zeros([m+1,n+1], dtype=float)

for i in range(1,m+1):
    for j in range(1,n+1):
        ImageJyppix_scaled_shifted[i,j] = ImageJyppix_scaled[i-1,j-1] 

plt.contourf(pix_x_shifted+0.5, pix_y_shifted+0.5, ImageJyppix_scaled_shifted, cmap=cm.gray)
plt.colorbar()
grids = "%s%s%s%s" % (nx+1, " x ", ny+1, " grids")
plt.title(grids)
plt.savefig('ori.png')
plt.clf()

 
###################### ###################### 
############ info printting
###################### ###################### 

print "pixel-size (x) in cm :"
print "     ", image.sizepix_x
print "pixel-size (y) in cm :"
print "     ", image.sizepix_y
print "distance from the source (cm), PC x dpc(nPCs) :"
print "     ", dpc*3.0856775e18
print "radian in x :"
print "     ", image.sizepix_x/dpc/3.0856775e18
print "radian in y :"
print "     ", image.sizepix_y/dpc/3.0856775e18


###################### ###################### 
####### set teh u_max,v_max for FFT
###################### ###################### 

u_max = u_wav/1e6/(image.sizepix_x/dpc/3.0856775e18)
v_max = v_wav/1e6/(image.sizepix_y/dpc/3.0856775e18)
print "u_Max (lamda_obs / radian_x) :"
print "     ", u_wav/1e6/(image.sizepix_x/dpc/3.0856775e18)
print "y_Max (lamda_obs / radian_y) :"
print "     ", v_wav/1e6/(image.sizepix_y/dpc/3.0856775e18), "\n"


###################### ###################### 
# Calcualte the Brightness Temperature
###################### ###################### 

if calc_TB:
    # calculate pixel size in cm^2
    pixel_size = image.sizepix_x * image.sizepix_y
    # calculate steradian 
    pixel_2_ster = pixel_size/140.0/140.0/pc/pc
    # convert Jy/pixel to Jy/sr
    ImageJyppix_scaled_Jy2sr = ImageJyppix_scaled/pixel_2_ster
    # calculate Brightness temperature
    # fomula1: T_B = Jy/sr * lamda(in meter)**2 / 2760
    # fomula2: T_B = S(flux density) * 0.32x10^23 * lamda**2 / theta**2   < mks units >
    # fomula3: T_B = S(mJy/beam) * 1.36 * lamda**2 / theta**2  ; < cm, seconds of arc >
    #           theta is the HPBW
    #Tb(K) = 13.5958 * (lambda/mm)^2 * (arcsec/bmaj)*(arcsec/bmin)*(flux/Jy) 
    #bmin=12.2842*1000 # mm
    #bmax=15238.4*1000 # mm
    #theta_bmin=(u_wav/1000/bmin)*206264.806247 # 1 radian = 206264.806247 arcsecond
    #theta_bmax=(u_wav/1000/bmax)*206264.806247
    #print "theta_bmin, theta_bmax, u_wav"
    #print theta_bmin, theta_bmax, u_wav, "\n"
    theta = image.sizepix_x/140/pc*206264.806247
    theta_bmin=theta
    theta_bmax=theta
    print theta
    #pixperbm = arcsec2perbm / arcsec2perpix
    T_B = 13.5958 * ImageJyppix_scaled * (u_wav/1000)**2.0 * (1.0/theta_bmin) * (1.0/theta_bmax)
    #T_B = ImageJyppix_scaled_Jy2sr * (u_wav*1e-6)**2.0 / 2760.0
    #T_B = np.array(T_B).reshape(m*n)


def Robs2Rreal(posang,incli):
    """
    calc the ratio of R_obs to R_real after a 'posang' degree 
         rotation of a disk with 'incli' degree inclination.
    NOTE: this ration is in the x-horizon direction
    """
    import math
    if (posang < 90) and (posang > 0):
        ang1 = 90.0 - posang
    elif (posang < 180) and (posang > 90):
        ang1 = posang - 90
    else:
        print "posang must between [0:180]"
        quit()
    ang1 = math.radians(ang1)
    incli = math.radians(incli)
    leng1 = math.sin(ang1)
    leng2 = math.cos(ang1)*math.cos(incli)
    print ang1
    print leng1, leng2
    ratio = math.sqrt(leng1**2.0+leng2**2.0)
    return ratio


# Output the flux in Jy/pixel for debug
if debug:
    image_output = np.zeros(m*n)
    for i in range(m):
        for j in range(n):
            image_output[n*i+j]=ImageJyppix_scaled[i,j]
    np.savetxt('ImageJyppix_scaled', np.transpose(image_output))


#####################
# save the mid-panel line for T_B plotting
ratio=Robs2Rreal(138.02,46.72)
print ratio
n_mid=(NN-1)/2
x_to_au = (pix_x+0.5)*image.sizepix_x/ratio/au
np.savetxt('T_B.out', np.transpose([T_B[n_mid,:], x_to_au]))
#####################

###################### ###################### 
######### FFT 
###################### ###################### 

print "Start FFT, this may take some time..."

FFT_image = np.fft.fft2(ImageJyppix_scaled_shifted)

FreqRows = np.fft.fftfreq(FFT_image.shape[0], fft_dx)
FreqCols = np.fft.fftfreq(FFT_image.shape[1], fft_dy)

FFT_image = np.fft.fftshift(FFT_image)   # Shift zero freq to center
FreqRows = np.fft.fftshift(FreqRows)   # Shift zero freq to center
FreqCols = np.fft.fftshift(FreqCols)   # Shift zero freq to center
#ImageJyppix_scaled = np.fft.ifft2(FFT_image)



#####################
# set the grids in AU
# change
im_x = (pix_x+0.5)*image.sizepix_x/au
im_y = (pix_x+0.5)*image.sizepix_y/au


#########################
#### Write the FFT of the image
#### Write the FFT of the image
u_array1 = np.linspace(-u_max/2.0, u_max/2.0, m, endpoint=False)
u_array2 = np.roll(u_array1, -1)
u_array2[m-1] = -u_array1[0]
u_array = np.zeros(m+1, dtype=float)
u_array[1:m+1:1] = u_array2
u_array[0] = u_array[1]-(u_array[2]-u_array[1])
v_array = u_array


# save u,v grid for interpolation 
np.savetxt('./u.grid', np.transpose(u_array))
np.savetxt('./v.grid', np.transpose(v_array))
np.savetxt('./rows.grid', np.transpose(FreqRows))
np.savetxt('./cols.grid', np.transpose(FreqCols))


fft_out = np.zeros(((m+1)*(n+1),4),  dtype=float)
# 0:u, 1:v, 2:R, 3:I
for i in range(m+1):
    for j in range(n+1):
        fft_out[i*(n+1)+j,0] = u_array[i]
        fft_out[i*(n+1)+j,1] = v_array[j] 
        fft_out[i*(n+1)+j,2] = np.real(FFT_image[i,j])*((-1)**(i+j))
        fft_out[i*(n+1)+j,3] = np.imag(FFT_image[i,j])*((-1)**(i+j))

uv_mid = np.zeros((NN+1,2),  dtype=float)
n_mid=(NN+1-1)/2
for i in range(NN+1):
    uv_mid[i,0] = fft_out[n_mid*(NN+1)+i,2] 
    uv_mid[i,1] = fft_out[n_mid*(NN+1)+i,3] 
np.savetxt('./u_v_mid.dat', np.transpose([uv_mid[:,0],uv_mid[:,1]]))

print "saving uvRI*.dat, this may take some time..."
#np.savetxt('./uvRI_model_av4.dat', np.transpose([fft_out[:,2],fft_out[:,3]]))
#np.savetxt('./uvRI_model_av4.dat', np.transpose([fft_out[:,0],fft_out[:,1],fft_out[:,2],fft_out[:,3]]))


############################
##### Plot 

plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.pcolormesh(pix_x+0.5, pix_y+0.5, ImageJyppix_scaled, cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("Janskey/pixel")
plt.title("Flux density")
plt.savefig('flux_density.png')

plt.clf()

plt.xlabel('Distance [AU]')
plt.ylabel('Distance [AU]')
plt.pcolormesh(im_x, im_y, T_B, cmap='RdBu', norm=LogNorm(), vmin=0.1, vmax=T_B.max())
plt.ylim(-180,180)
plt.xlim(-180,180)
cbar1=plt.colorbar()
cbar1.set_label("Temperature [ K ]")
plt.title("Brightness Temmperature")
plt.savefig('brightness_temperature.png')

plt.clf()

fft_real = fft_out[:,2]    #T_B = np.array(T_B).reshape(m*n)
fft_real = np.array(fft_real).reshape(m+1,n+1)
plt.xlabel('u')
plt.ylabel('v')
plt.ylim(-500,500)
plt.xlim(-500,500)
plt.pcolormesh(u_array, v_array, fft_real)
#, norm=LogNorm(), vmin=fft_out[:,2].max()*1e-6, vmax=fft_out[:,2].max())
cbar2=plt.colorbar()
cbar2.set_label("Real")
plt.grid()
plt.savefig('image_fft_real_uv.png')

plt.clf()

fft_imag = fft_out[:,3]    #T_B = np.array(T_B).reshape(m*n)
fft_imag = np.array(fft_imag).reshape(m+1,n+1)
plt.xlabel('u')
plt.ylabel('v')
plt.ylim(-500,500)
plt.xlim(-500,500)
plt.pcolormesh(u_array, v_array, fft_imag)
plt.grid()
#, norm=LogNorm(), vmin=fft_out[:,3].max()*1e-6, vmax=fft_out[:,3].max())
cbar2=plt.colorbar()
cbar2.set_label("Imag")
plt.savefig('image_fft_imag_uv.png')
quit()


