import os

os.system('g++ src/buckling_beam.cpp -larmadillo')
os.system('g++ src/quantum_dots_one_electron.cpp -larmadillo')
os.system('g++ src/quantum_dots_two_electrons.cpp -larmadillo')

if os.path.isfile('a.out'):
    os.system('./a.out')
    os.system('rm a.out')
