import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
filename = sys.argv[1] #file name passed as argument

timing_data = pd.read_csv(filename)

print(timing_data)

# --- PLOT ONLY TOTAL, COMPUTATION AND COMMUNICATION TIMES
#timing_plot = timing_data.plot(x="nproc",y=['elapsed_time','time_comput','time_comm'],loglog="True",style='-o')

# --- PLOT TIME FOR COMMUNICATION ROUTINES
timing_plot = timing_data.plot(x="nproc",y=timing_data.columns[4:],loglog="True",style='-o')

# --- PLOT EVERYTHING
#timing_plot = timing_data.plot(x="nproc",loglog="True",style='-o')


plt.ylabel('time')
plt.title('Timing of a parallel std_dev program')
plt.grid(True)

fig = timing_plot.get_figure()
fig.savefig("timing.png")
plt.show()




