import global_coupling_largeN as simlarge
import numpy as np
import sys

# ''' MAIN '''
if __name__ == '__main__':
	''' MAIN:

    Variables
    ---------
    N		   : integer
				 number of oscillators

    Note: all frequencies are in Hz, convertion to radHz automatically in the program

	Returns
	-------


	------------------------------------------------------------------------------------
	For adiabtic change cases: going up by starting with equally distribued phases, when
	going down, starting from the perfectly synchronized state.
    '''
# simulate for different value of the delay and the integration time/inertia
tauValues 		= np.linspace(3,3,1);
fcValues  		= np.linspace(0.05305,0.05305,1);
KstartValues	= np.linspace(0.002,0.4,20);

for n in range(len(fcValues)):
	for m in range(len(tauValues)):
		for o in range(len(KstartValues)):

			tau          = tauValues[n];		# transmission delay in seconds        e.g., 3.00, 4.18
			fc           = fcValues[n];			# cut-off frequency in Hz [fc=1/(2*pi*m)], e.g., m=3->fc=0.05305, m=0.1->fc=1.5915
			init_phases  = 'iid'				# specifies initial phases of the PLLs [zero, gaussian, iid]
			init_hist    = 'free'				# specifies initial history for simulation [free, synced, synced_highR]
			int_frequen  = 'lorentz'			# specifies dist. of intrinsic frequencies [const, gaussian, lorentz]
			N            = 25000;  			    # number of PLLs in the system
			Tsim         = 10000;  			    # simulation time in seconds, or Trelax the relaxation time (adiabtic case)
			Kstart       = KstartValues[o];		# starting coupling strength in Hz(/V)
			Kend         = 0.4;					# maximum coupling strength in Hz(/V), e.g., 1.1/(2pi), 0.25/(2pi)
			fc           = fcValues[n];			# cut-off frequency in Hz [fc=1/(2*pi*m)], e.g., m=3->fc=0.05305, m=0.1->fc=1.5915
			tau          = tauValues[m];  		# transmission delay in seconds        e.g., 3.00, 4.18
			m_init_phase = 3.14;				# mean initial phases of PLLs [rad]]
			var_init_phi = 6.28;				# variance of gaussian or iid distribution of initial phases
			m_int_freqs  = 0.4779;				# mean intrinsic frequency of the PLLs in Hz
			v_int_freqs  = 0.0159;				# variance of gaussian or scale of lorentz distributed int. frequencies
			sync_freq    = 0.4779;				# initial frequency of the system in case of coupled initial state [Hz]
			SamplesPerT  = 250;					# samples per period, determines time-step for iteration of Euler-sheme / Runge-Kutta 4th order
			adiabatic    = 0;					# whether or not adiabatic change of variable
			ad_rate      = 0.0005; 				# change rate Krate = ad_rate / Tsim
			downscale    = 10;					# divisor of sample rate when plotting and saving data
			plot_phi_freq= 0;					# whether to plot phases and frequencies
			plot_out     = 0;					# whether to open plots in matplotlib windows
			wrap_hist_tau= 1; 					# whether time-series wrapped into vector of length tau --> evolveOrderReduced vs evolveOrder
			specified_distribution = 1; 		# if this value is non-zero, a specific distribution of intrinsic frequencies is drawn via random seed

			simlarge.main(init_phases, init_hist, int_frequen, N, Tsim, Kstart, Kend, fc, tau, downscale, m_init_phase,
			         var_init_phi, m_int_freqs, v_int_freqs, sync_freq, SamplesPerT, adiabatic, ad_rate, specified_distribution,
					 plot_phi_freq, wrap_hist_tau, plot_out)

# init_phases  = 'iid'                   # specifies initial phases of the PLLs [zero, gaussian, iid]
# init_hist    = 'free'                  # specifies initial history for simulation [free, synced]
# int_frequen  = 'const'               # specifies dist. of intrinsic frequencies [const, gaussian, lorentz]
# N            = 100;                  # number of PLLs in the system
# Tsim         = 50;                    # simulation time in seconds, or Trelax the relaxation time (adiabtic case)
# Kstart       = 0.002;                  # starting coupling strength in Hz(/V)
# Kend         = 5;                    # maximum coupling strength in Hz(/V)
# fc           = 0.1;                # cut-off frequency in Hz [fc=1/(2*pi*m)]
# tau          = 1.0;                   # transmission delay in Hz
# m_init_phase = 3.14;                   # mean initial phases of PLLs [rad]]
# var_init_phi = 0.1;                    # variance of gaussian or iid distribution of initial phases
# m_int_freqs  = 3;                 # mean intrinsic frequency of the PLLs in Hz
# v_int_freqs  = 0.0001;            # variance of gaussian or scale of lorentz distributed int. frequencies
# sync_freq    = 3;                 # initial frequency of the system in case of coupled initial state [Hz]
# SamplesPerT  = 100;                    # samples per period, determines time-step for iteration of Euler-sheme / Runge-Kutta 4th order
# adiabatic    = 1; 					   # whether or not adiabatic change of variable
# ad_rate      = 0.01;				   # change rate Krate = ad_rate / Tsim
# downscale    = 1;					   # divisor of sample rate
# plot_phi_freq= 1;					   # whether to plot phases and frequencies
# plot_out     = 1;                      # whether to open plots in matplotlib windows

# 	''' SIMULATION PARAMETER'''
# init_phases  = str(sys.argv[1])      # specifies initial phases of the PLLs [zero, gaussian, iid]
# init_hist    = str(sys.argv[2])      # specifies initial history for simulation [free, synced]
# int_frequen  = str(sys.argv[3])      # specifies dist. of intrinsic frequencies [const, gaussian, lorentz]
# N            = int(sys.argv[4])      # number of PLLs in the system
# Tsim         = float(sys.argv[5]);   # simulation time in seconds, or Trelax the relaxation time (adiabtic case)
# Kstart       = float(sys.argv[6]);   # starting coupling strength in Hz(/V)
# Kend         = float(sys.argv[7]);   # maximum coupling strength in Hz(/V)
# fc           = float(sys.argv[8]);   # cut-off frequency in Hz
# tau          = float(sys.argv[9]);   # transmission delay in Hz
# m_init_phase = float(sys.argv[10]);  # mean initial phases of PLLs
# var_init_phi = float(sys.argv[11]);  # variance of gaussian or iid distribution of initial phases
# m_int_freqs  = float(sys.argv[12]);  # mean intrinsic frequency of the PLLs in Hz
# v_int_freqs  = float(sys.argv[13]);  # variance of gaussian or scale of lorentz distributed int. frequencies
# sync_freq    = float(sys.argv[14]);  # initial frequency of the system in case of coupled initial state
# SamplesPerT  = int(sys.argv[15]);    # samples per period, determines time-step for iteration of Euler-sheme / Runge-Kutta 4th order
# adiabatic    = float(sys.argv[16]);  # whether or not adiabatic change of variable
# ad_rate      = float(sys.argv[17]);  # change rate Krate = ad_rate / Tsim
# downscale    = int(sys.argv[18]);    # divisor of sample rate
# plot_phi_freq= int(sys.argv[19]);    # whether to plot phases and frequencies
# plot_out     = int(sys.argv[20]);    # whether to open plots in matplotlib windows

#sim.main(init_phases, init_hist, int_frequen, N, Tsim, Kstart, Kend, fc, tau, m_init_phase,
#         var_init_phi, m_int_freqs, v_int_freqs, sync_freq, dt, adiabatic, plot_phi_freq, plot_out)
