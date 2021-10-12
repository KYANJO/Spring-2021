import matplotlib.pyplot as plt
import numpy as np
import sys

for filenum in range(len(sys.argv)-1):

    filename = str(sys.argv[filenum+1])
    print "Reading data from file:", filename
    n, err, time = np.loadtxt(filename, delimiter='\t', unpack=True,skiprows=0)
    plt.loglog(n,time,'-',linewidth=2,label=filename);

plt.xlabel('# of points')
plt.ylabel('time (s)')
plt.title('Timing')
plt.legend()
plt.grid(True)
plt.savefig("timing_mpi.png");
plt.show()

