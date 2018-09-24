import os, sys
import matplotlib.pyplot as plt
from math import log10

dirname = os.path.dirname(sys.argv[0])
if len(dirname) > 0:
    dirname += '/'

# Running C++ program
os.popen('c++ tridiagonal.cpp')
os.system('./' + dirname + 'a.out > out/tridiagonal_output.dat')
#os.popen('c++ LU.cpp')
#os.system('./' + dirname + 'a.out > out/LU_output.dat')
# os.popen()

# Preparing plots
plt.xlabel('x')
plt.ylabel('u(x)')

# Making plots
for n in [10, 100, 1000]:
    infile = open('out/tridiagonal_n=%g.dat' % n, 'r')

    x = []
    analytic = []
    numeric  = []

    for line in infile:
        numbers = line.split()
        x.append(eval(numbers[0]))
        analytic.append(eval(numbers[1]))
        numeric.append(eval(numbers[2]))

    plt.plot(x, analytic, label='Analytic')
    plt.plot(x, numeric, label='Numeric')
    plt.title('Calculation of u(x) for n=%g' % n)
    plt.legend()
    plt.savefig('plots/tridiagonal_n=%g.png' % n)
    plt.clf()

    infile.close()

# Collecting data for table
infile = open('out/tridiagonal_output.dat', 'r')

epsilons = []
loghs    = []

for i in range(1, 8):
    infile.readline()
    infile.readline()
    epsilons.append(eval(infile.readline().strip().split()[-1]))
    infile.readline()
    loghs.append(eval(infile.readline().strip().split()[-1]))
    infile.readline()

infile.close()

# Making table
outfile = open('out/project1d_table.dat', 'w')



outfile.write('\\begin{tabular}{ c c c } ')
outfile.write(' '*9 + '& $log(h)$ & $\epsilon_{max}$ \\\\ \n')
outfile.write('\\hline')
for i, e, l in zip(range(1, 8), epsilons, loghs):
    outfile.write(' $n=10^%1.0f$ & %5.2f & %8.5f \\\\ \n' % (i, l, e))
outfile.write('\\end{tabular}\n')

outfile.close()
