#!/bin/python

import os 

print('Cloning RelicFast repository...') 
os.system('git clone --branch nick git@github.com:JulianBMunoz/RelicFast.git')
os.chdir('./RelicFast')

print('Setting axionCAMB as Boltzmann Solver in common.h file...') 
reading_file = open("./include/common.h", "r")
new_file_content = ""
for line in reading_file:
  stripped_line = line.strip()
  new_line = stripped_line.replace(
        "boltzmann_tag  _CLASS_", "boltzmann_tag  _AXIONCAMB_"
    ).replace(
        "boltzmann_tag  _CAMB_", "boltzmann_tag  _AXIONCAMB_"
    )
  new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./include/common.h') 
writing_file = open("./include/common.h", "w")
writing_file.write(new_file_content)
writing_file.close()

print('Compiling CLASS Boltzmann solver...')
os.system('cp -r ./CLASS ./CLASS_Current')          
os.chdir('./CLASS_Current')
reading_file = open("./Makefile", "r")
new_file_content = ""
for line in reading_file:
  stripped_line = line.strip()
  new_line = stripped_line.replace(
    "CC       = gcc", "CC       = gcc-10"
  ).replace(
    "OMPFLAG   = -openmp", ""
  )
  new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./Makefile')
writing_file = open("./Makefile", "w")
writing_file.write(new_file_content)
writing_file.close()
os.system('make all')
os.chdir('../')

print('Compiling CAMB Boltzmann solver...')
os.system('cp -r ./CAMB ./CAMB_Current') 
os.chdir('./CAMB_Current')
print('Setting fortran as compiler for CAMB code...') 
reading_file = open("./Makefile", "r")
new_file_content = ""
for line in reading_file:
  stripped_line = line.strip()
  new_line = stripped_line.replace(
    "F90C     = ifort", "F90C     = gfortran"
  )
  new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./Makefile') 
writing_file = open("./Makefile", "w")
writing_file.write(new_file_content)
writing_file.close()
os.system('make all') 
os.chdir('../') 

print('Compiling axionCAMB Boltzmann solver...')
os.system('cp -r ./axionCAMB ./axionCAMB_Current') 
os.chdir('./axionCAMB_Current')
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
os.system('make all') 
os.chdir('../')

print('Removing openmp flag from RelicFast compilation...') 
reading_file = open("./Makefile", "r")
new_file_content = ""
for line in reading_file:
    new_line = line.replace(
        "parallel = -fopenmp", ""
    ).replace(
        "CC=gcc", "CC=gcc-10"
    )
    new_file_content += new_line +"\n"
reading_file.close()
os.system('rm ./Makefile') 
writing_file = open("./Makefile", "w")
writing_file.write(new_file_content)
writing_file.close()

print('Compiling RelicFast...') 
os.system('make') 

