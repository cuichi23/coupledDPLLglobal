import datetime
import numpy as np
import matplotlib.pyplot as plt
import csv

''' All plots in latex mode '''
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams['agg.path.chunksize'] = 10000

''' STYLEPACKS '''
titlefont = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 10,
        }

labelfont = {
        'family' : 'sans-serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
        }

annotationfont = {
        'family' : 'monospace',
        'color'  : (0, 0.27, 0.08),
        'weight' : 'normal',
        'size'   : 16,
        }


data_set_two = 'False';
method       = 'C';                                                             # choose method B (with R(0)->0 or R(0)->1) or C (K+dK with last K's final state)
plot_adiabtic= 0;

# load data from param-sweep
with open('sweep_results.txt', 'rt') as f:                                      # _N10000_Tsim25000
    reader = csv.reader(f, delimiter=',', skipinitialspace=True)
    # to skip header: next(reader, None)
    lineData = list()
    cols = next(reader)
    print(cols)

    for col in cols:
    # Create a list in lineData for each column of data.
        lineData.append(list())


    for line in reader:
        for i in range(0, len(lineData)):
        # Copy the data from the line into the correct columns.
            lineData[i].append(line[i])

    data = dict()

    for i in range(0, len(cols)):
    # Create each key in the dict with the data in its column.
        data[cols[i]] = lineData[i]

Kstart_data=np.array(data['Kstart'], dtype=float)
lastRt_data=np.array(data['lastR'], dtype=float)
meanRt_data=np.array(data['meanR'], dtype=float)

# load data from param-sweep
if data_set_two == 'True':
    with open('sweep_results_highR.txt', 'rt') as f:
        reader = csv.reader(f, delimiter=',', skipinitialspace=True)
        # to skip header: next(reader, None)
        lineData = list()
        cols = next(reader)
        #print(cols)

        for col in cols:
        # Create a list in lineData for each column of data.
            lineData.append(list())


        for line in reader:
            for i in range(0, len(lineData)):
            # Copy the data from the line into the correct columns.
                lineData[i].append(line[i])

        data1 = dict()

        for i in range(0, len(cols)):
        # Create each key in the dict with the data in its column.
            data1[cols[i]] = lineData[i]

    Kstart_data1=np.array(data1['Kstart'], dtype=float)
    lastRt_data1=np.array(data1['lastR'], dtype=float)
    meanRt_data1=np.array(data1['meanR'], dtype=float)
#print(Kstart_data)
#print(lastRt_data)

# load data from adiabatic simulation
if plot_adiabtic == 1:
    dataRt=np.load('orderparam_Kstart0.002_Fc0.05305_tau3.0_2018_10_1.npz')
    dataK =np.load('Kt_Kstart0.002_Fc0.05305_tau3.0_2018_10_1.npz')

    Rt    = dataRt['Rt']
    K     = dataK['Kt']

init_phases  = 'iid'		# specifies initial phases of the PLLs [zero, gaussian, iid]
init_hist    = 'free'		# specifies initial history for simulation [free, synced]
int_frequen  = 'lorentz'	# specifies dist. of intrinsic frequencies [const, gaussian, lorentz]
N            = 10000;  		# number of PLLs in the system
Tsim         = 1500;  		# simulation time in seconds, or Trelax the relaxation time (adiabtic case)
Kstart       = 0.002;		# starting coupling strength in Hz(/V)
Kend         = 0.4;			# maximum coupling strength in Hz(/V), e.g., 1.1/(2pi), 0.25/(2pi)
fc           = 0.05305;		# cut-off frequency in Hz [fc=1/(2*pi*m)], e.g., m=3->fc=0.05305, m=0.1->fc=1.5915
tau          = 3.00;		# transmission delay in seconds        e.g., 3.00, 4.18
m_init_phase = 3.14;		# mean initial phases of PLLs [rad]]
var_init_phi = 6.28;		# variance of gaussian or iid distribution of initial phases
m_int_freqs  = 0.4779;		# mean intrinsic frequency of the PLLs in Hz
v_int_freqs  = 0.0159;		# variance of gaussian or scale of lorentz distributed int. frequencies
sync_freq    = 0.4779;		# initial frequency of the system in case of coupled initial state [Hz]
SamplesPerT  = 125;			# samples per period, determines time-step for iteration of Euler-sheme / Runge-Kutta 4th order
adiabatic    = 1;			# whether or not adiabatic change of variable
ad_rate      = 0.001; 		# change rate Krate = ad_rate / Tsim
downscale    = 10;			# divisor of sample rate when plotting and saving data
plot_phi_freq= 0;			# whether to plot phases and frequencies
plot_out     = 0;			# whether to open plots in matplotlib windows
down         = 10;			# see simulation

if tau == 0:
    dt = 1.0 / ( SamplesPerT * m_int_freqs );
elif tau > 0:
    dt = 1.0 / ( SamplesPerT * sync_freq );
    tausteps = int(round(tau/dt));

Krate = ad_rate / Tsim;                   # the rate at which K(t) changes adiabatically
Treve = 0;                                # after relaxation, the time until K(t) growth reverses
if Krate > 0:
    Treve = (Kend - Kstart) / Krate;
# time-steps and container for histories, coupling strengths and phase time-series
trelsteps = int(round(Tsim/dt));              # relaxation time in steps dt
trevsteps = 0;                                # define time of adiabatic change in steps dt
extrsteps = 0;                                # define time of adiabatic change in steps dt
if adiabatic == 1:
    trevsteps = int(round(2.0*Treve/dt));     # time of adiabatic change in steps dt
    extrsteps = int(round(0.5*Tsim/dt));      # time of adiabatic change in steps dt
itersteps = trelsteps+trevsteps+extrsteps;    # total simulation time in dt

t = np.arange(tausteps+itersteps, dtype=np.float64);
t = t[::down];

now = datetime.datetime.now()
#******************************************************************************************************************
print('Now plot!')

# plt.rcParams['figure.figsize'] = (1.25*30, 1.25*18)
# plt.figure('orderparameter vs K'); plt.clf();
# #plt.plot(2.0 * np.pi * K, Rt);
# plt.plot(2.0*np.pi*K[0:int(tausteps/down)], Rt[0:int(tausteps/down)],'g--');
# plt.plot(2.0*np.pi*K[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)], Rt[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)],'r--');
# plt.plot(2.0*np.pi*K[int((tausteps+trelsteps+0.5*trevsteps)/down):], Rt[int((tausteps+trelsteps+0.5*trevsteps)/down):],'b-');
# #plt.plot([2.0*np.pi*0.23,2.0*np.pi*0.223,2.0*np.pi*0.302,2.0*np.pi*0.27],[0.41,0.38561,0.69434,0.59477],'g*')
# # plot 1/N order parameter
# plt.plot(2.0*np.pi*Kstart_data,(np.zeros(len(Kstart_data))+1)*1/N,'k-');
# # plot data
# plt.plot(2.0*np.pi*Kstart_data, lastRt_data, 'cD--', ms=5)
# plt.plot(2.0*np.pi*Kstart_data, meanRt_data, 'r+--', ms=6)
# # plot data 1
# if data_set_two == 'True':
#     plt.plot(2.0*np.pi*Kstart_data1, lastRt_data1, 'cD-', ms=5)
#     plt.plot(2.0*np.pi*Kstart_data1, meanRt_data1, 'r+-', ms=6)
#
# plt.xlabel(r'$K(t)$ [radHz]', fontdict = labelfont);
# plt.ylabel(r'$R(t)$', fontdict = labelfont);
# plt.savefig('order_vs_Kt_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
# plt.savefig('order_vs_Kt_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)

fig, ax = plt.subplots(figsize=(0.5*30, 0.5*18));
if plot_adiabtic == 1:
    ax.plot(2.0*np.pi*K[0:int(tausteps/down)], Rt[0:int(tausteps/down)],'g--');
    ax.plot(2.0*np.pi*K[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)], Rt[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)],'r--');
    ax.plot(2.0*np.pi*K[int((tausteps+trelsteps+0.5*trevsteps)/down):], Rt[int((tausteps+trelsteps+0.5*trevsteps)/down):],'b-');
# plot 1/N order parameter
ax.plot(2.0*np.pi*Kstart_data,(np.zeros(len(Kstart_data))+1)*1/np.sqrt(N),'k-',label=r'$1 \sqrt{N}$ treshold');
# plot data
if method == 'B':
    ax.plot(2.0*np.pi*Kstart_data, lastRt_data, 'cD--', ms=5, label=r'last $R(t)$, $R_{init}\rightarrow 0$')
    ax.plot(2.0*np.pi*Kstart_data, meanRt_data, 'r+--', ms=6, label=r'mean $R(t)$, $\bar{R}_{init}\rightarrow 0$')
elif method == 'C':
    ax.plot(2.0*np.pi*Kstart_data, lastRt_data, 'cD--', ms=5, label=r'last $R(t)$, $\bar{R}_{init}\rightarrow$ final state last K(n)')
    ax.plot(2.0*np.pi*Kstart_data, meanRt_data, 'r+--', ms=6, label=r'mean $R(t)$, $\bar{R}_{init}\rightarrow$ final state last K(n)')
# plot data 1
if data_set_two == 'True':
    if method == 'B':
        ax.plot(2.0*np.pi*Kstart_data1, lastRt_data1, 'cD-', ms=5, label=r'last $R(t)$, $R_{init}\rightarrow 1$')
        ax.plot(2.0*np.pi*Kstart_data1, meanRt_data1, 'r+-', ms=6, label=r'mean $R(t)$, $\bar{R}_{init}\rightarrow 1$')
    elif method == 'C':
        ax.plot(2.0*np.pi*Kstart_data1, lastRt_data1, 'cD-', ms=5, label=r'last $R(t)$, $\bar{R}_{init}\rightarrow$ final state last K(n)')
        ax.plot(2.0*np.pi*Kstart_data1, meanRt_data1, 'r+-', ms=6, label=r'mean $R(t)$, $\bar{R}_{init}\rightarrow$ final state last K(n)')
ax.legend(loc='upper left', shadow=True, fontsize='large')
ax.set_xlabel(r'$K(t)$ [radHz]', fontdict = labelfont);
ax.set_ylabel(r'$R(t)$', fontdict = labelfont);
plt.savefig('order_vs_Kt_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
plt.savefig('order_vs_Kt_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)

if plot_adiabtic == 1:
    fig, ax1 = plt.subplots();
    ax1.plot(t[int(tausteps/down)]*dt, 0.05, 'ro', ms=2)
    if adiabatic == 1:
        ax1.plot(t[int((tausteps+trelsteps)/down)]*dt, 0.05, 'yo', ms=2)

    ax1.plot(t*dt,Rt,'b-');
    ax1.set_xlabel(r'$t$ $[s]$', fontdict = labelfont);
    ax1.set_ylabel(r'$R(t)$', fontdict = labelfont, color='b');
    ax1.tick_params('y', colors='b');
    ax2 = ax1.twinx();
    ax2.plot(t*dt, 2.0*np.pi*K, 'r-');
    ax2.set_ylabel(r'$K(t)$ [radHz]', fontdict = labelfont, color='r'); #y-axis match the line color
    ax2.tick_params('y', colors='r');

    fig.tight_layout();
    plt.savefig('order_vs_t_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
    plt.savefig('order_vs_t_ADD_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)

plt.draw()
plt.show()
