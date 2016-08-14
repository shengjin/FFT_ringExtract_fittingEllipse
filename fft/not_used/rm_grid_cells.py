#! /usr/bin/env python2

import numpy as np

  
u_grid = np.genfromtxt('./u.grid.middle', dtype=float)
n_u = u_grid.shape[0]
n_data = n_u*n_u

R_mod = np.zeros(n_data,dtype=float)

print "Reading files"

f = open('I.dat', 'r')
for i in range(n_data):
    if(i%max(1,int(0.01*n_data))==0):
        print "%s%s" % (100*i/n_data, "% done ...")
    tmp = float(f.readline())
    R_mod[i] = tmp

# reshape the R,I model
R_mod = R_mod.reshape(n_u,n_u)
print "R_mod reshaped to: "
print "     ", R_mod.shape

n_u_new = int((n_u-1)/2)+1
u_grid_new = np.zeros(n_u_new, dtype=float)
R_new = np.zeros([n_u_new, n_u_new], dtype=float)

print n_u_new, n_u

j = 0


for i in range(n_u):
    if(i%2==0):
        u_grid_new[j] = u_grid[i]
        j = j + 1

print "sift data"
for i in range(n_u_new):
    ## print a progress bar
    if(i%max(1,int(0.01*n_u_new))==0):
        print "%s%s" % (100*i/n_u_new, "% done ...")
    for j in range(n_u_new):
        R_new[i,j] = R_mod[i*2,j*2]

print R_new[50,20], R_mod[100,40]

R_new = R_new.reshape(n_u_new*n_u_new)

np.savetxt('v.grid.new', np.transpose(u_grid_new))
np.savetxt('I.new', np.transpose(R_new))

