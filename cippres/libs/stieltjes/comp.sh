#!/bin/bash

#ifort -c -C -traceback stieltjes.f90 imaging.f90 interp.f90 pythag.f sort2.f90 tql2.f 
#ifort -C -traceback stieltjes.o imaging.o interp.o pythag.o sort2.o tql2.o -o stieltjes 

ifort -c stieltjes.f90 imaging.f90 interp.f90 pythag.f sort2.f90 tql2.f 
ifort  stieltjes.o imaging.o interp.o pythag.o sort2.o tql2.o -o stieltjes 

#ifort -c stieltjes.f90 imaging.f90 interp.f90 pythag.f sort2.f90 tql2.f -C -traceback
#ifort stieltjes.o imaging.o interp.o pythag.o sort2.o tql2.o -o stieltjes -C -traceback
