# 
# Script to azimuthally average brightness temperature of an image:
#

import numpy as np
import matplotlib.pyplot as plt

# Get file names
image_name = 'obs.fits'
output_text = 'Tb_vs_R.txt'
image_center_RA = '04:31:38.425'
image_center_Dec = '+018.13.57.241'
disk_PA = '138.02'
#disk_incl = 0
disk_incl = 46.72
distance = 140.0
smooth_beam = 0.0335
#smooth_beam = 0.031
smooth_rms  = 0.000186679

#radius_wanted = 40.0
#radius_wanted = 68.0
#radius_wanted = 32.9
radius_wanted = 74.0


# Convert to fits file if necessary
if image_name.endswith('.fits'):
    output_name = image_name[:-5] + '.image'
    importfits(fitsimage=image_name,imagename=output_name,overwrite=True)
    # Smooth image to circular beam (for azimuthal averaging)
    imsmooth(imagename=output_name, outfile=image_name[:-5]+'_smooth.image',
             kernel='gauss',
             major=str(smooth_beam)+'arcsec',
             minor=str(smooth_beam)+'arcsec',
             pa='0.0deg',targetres=True,overwrite=True)
    ia.open(image_name[:-5]+'_smooth.image')
    
else:
    # Smooth image to circular beam (for azimuthal averaging)
    imsmooth(imagename=image_name, outfile=image_name[:-6]+'_smooth.image',
             kernel='gauss',
             major=str(smooth_beam)+'arcsec',
             minor=str(smooth_beam)+'arcsec',
             pa='0.0deg',targetres=True,overwrite=True)
    ia.open(image_name[:-6]+'_smooth.image')


# Get beam sizes
unit_bmaj = ia.summary()['restoringbeam']['major']['unit']
unit_bmin = ia.summary()['restoringbeam']['minor']['unit']
value_bmaj = ia.summary()['restoringbeam']['major']['value']
value_bmin = ia.summary()['restoringbeam']['minor']['value']

print "unit_bmaj", unit_bmaj
print "unit_bmin", unit_bmin
print "value_bmaj", value_bmaj
print "value_bmin", value_bmin, "\n"


if unit_bmaj == 'arcsec' and unit_bmin == 'arcsec':
    a_maj = value_bmaj
    a_min = value_bmin
else:
    print "Units of major/minor axis are *not* arcseconds!"

# Get image size
size_x = ia.shape()[0]
size_y = ia.shape()[1]
print "size_x", size_x
print "size_y", size_y, "\n"

# Get pixel size
pix_size_x = abs(ia.summary()['incr'][0])*180./pi*3600.
pix_size_y = abs(ia.summary()['incr'][1])*180./pi*3600.
print "pix_size_x", pix_size_x
print "pix_size_y", pix_size_y, "\n"

# Pixel area
pix_area = pix_size_x * pix_size_y

# Define steps radially
step_size = a_min/4.0    # Double the Nyquist rate
inner_rad = np.arange(step_size,2.0,step_size) # keeped for future use in case want to output all the rad
outer_rad = np.arange(inner_rad[1],2.0,step_size)
print "step_size", step_size, "\n"
#print "inner_rad", inner_rad
#print "outer_rad", outer_rad


rad_wanted = radius_wanted/distance
for i in range(np.size(outer_rad)):
    if (outer_rad[i] < rad_wanted) and (outer_rad[i+1] > rad_wanted):
        n_rad_wanted = i+1
        break
print rad_wanted
print n_rad_wanted, outer_rad[i], outer_rad[i+1]

#radius_wanted = (inner_rad-step_size/2.)*distance


major_large = str(outer_rad[n_rad_wanted+1])
minor_large = str(outer_rad[n_rad_wanted+1] * np.cos(disk_incl*pi/180.))
major_small = str(outer_rad[n_rad_wanted])
minor_small = str(outer_rad[n_rad_wanted] * np.cos(disk_incl*pi/180.))
ellipse_region_large = "ellipse [[" + image_center_RA + ", " + image_center_Dec + \
                       "], [" + major_large + "arcsec, " + minor_large + "arcsec], "+ disk_PA +"deg]"
ellipse_region_small = "ellipse [[" + image_center_RA + ", " + image_center_Dec + \
                       "], [" + major_small + "arcsec, " + minor_small + "arcsec], "+ disk_PA +"deg]"
# Get image coordinates and set to region
csys=ia.coordsys()
print "csys", csys, "\n"
rg.setcoordinates(csys.torecord())
print "torecord", csys.torecord(), "\n"
print "rg.setcoordinates(csys.torecord())", rg.setcoordinates(csys.torecord()), "\n"
a=rg.fromtext(ellipse_region_large, shape=[size_x, size_y,    1, 1])
b=rg.fromtext(ellipse_region_small, shape=[size_x, size_y,    1, 1])
bcomp=rg.complement(b)
regs={'1':a, '2':bcomp}
inter=rg.intersection(regs)
ib=ia.subimage(outfile='test_Snu.image', region=inter, overwrite=True)
ib.done()
imstat_Snu = imstat('test_Snu.image')

ib_l=ia.subimage(outfile='test_Snu_l.image', region=ellipse_region_large, overwrite=True)
ib_l.done()

# Number of independent beams in annuli:
npoints = imstat_Snu['npts'][0]
print "npoints", npoints
area_of_annuli = npoints * pix_area
print "area_of_annuli", area_of_annuli
beam_correction = ( area_of_annuli / (pi*smooth_beam**2.0) )
print "beam_correction", beam_correction

mean_Snu=imstat_Snu['mean'][0]
sigma_Snu=smooth_rms/(beam_correction)**0.5
print mean_Snu, sigma_Snu

intensity=imval(imagename='test_Snu.image', region=ellipse_region_large)
intensity_data = intensity['data'][:]
intensity_mask = intensity['mask'][:]
intensity_coords = intensity['coords'][:]
print intensity_data.shape
print intensity_coords.shape

mask_true = np.where(intensity_mask == True)
fin_data = intensity_data[mask_true]
ra_all=intensity_coords[:,:,0]
dec_all=intensity_coords[:,:,1]
fin_ra = ra_all[mask_true]
fin_dec = dec_all[mask_true]

np.savetxt
np.savetxt('JyPixel_Ra_Dec.out', np.transpose([fin_data,fin_ra,fin_dec]))

