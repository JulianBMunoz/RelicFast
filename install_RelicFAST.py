#!/bin/python

import os 

print('Cloning RelicFast repository...') 
os.system('git clone --branch nick --recurse-submodules git@github.com:JulianBMunoz/RelicFast.git')
os.chdir('./RelicFast')
#os.system('git submodule init')
#os.system('git submodule update')

print('Fixing length_transfer_camb parameter in common.h file...') 
reading_file = open("./include/common.h", "r")
new_file_content = ""
for line in reading_file:
  stripped_line = line.strip()
  new_line = stripped_line.replace("length_transfer_camb 700", "length_transfer_camb 595")
  new_line = stripped_line.replace("boltzmann_tag  _CLASS_", "boltzmann_tag  _CAMB_")
  new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./include/common.h') 
writing_file = open("./include/common.h", "w")
writing_file.write(new_file_content)
writing_file.close()

print('Setting axionCAMB as current Boltzmann solver...') 
os.system('cp -r ./axionCAMB ./CAMB_Current') 
os.chdir('./CAMB_Current')

print('Setting fortran as compiler for axionCAMB code...') 
reading_file = open("./Makefile", "r")
new_file_content = ""
for line in reading_file:
  stripped_line = line.strip()
  new_line = stripped_line.replace("F90C     = ifort", "F90C     = gfortran")
  new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./Makefile') 
writing_file = open("./Makefile", "w")
writing_file.write(new_file_content)
writing_file.close()

print('Compiling axionCAMB...') 
os.system('make all') 

print('Removing openmp from RelicFast compilation...') 
os.chdir('..') 
reading_file = open("./Makefile", "r")
new_file_content = ""
for line in reading_file:
    new_line = line.replace("parallel = -fopenmp", "")
    new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./Makefile') 
writing_file = open("./Makefile", "w")
writing_file.write(new_file_content)
writing_file.close()

print('Compiling RelicFast...') 
os.system('make') 

