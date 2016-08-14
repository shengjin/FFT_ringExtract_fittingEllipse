#!/bin/bash
awk '{print $1}' RI_model.dat > ../../R.dat 
awk '{print $2}' RI_model.dat > ../../I.dat 
cp v.grid ../../v.grid
cp u.grid ../../u.grid
