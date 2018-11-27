#python related imports
import os
import sys
import time
import os.path
import datetime
import numpy as np
import matplotlib
from pathlib import Path
if not os.environ.get('SGE_ROOT') == None: # this environment variable is set within the queue network, i.e., if it exists, 'Agg' mode to supress output
    print('NOTE: \"matplotlib.use(\'Agg\')\"-mode active, plots are not shown on screen, just saved to results folder!\n')
    matplotlib.use('Agg') #'%pylab inline'
import matplotlib.pyplot as plt
from scipy.stats import cauchy
import gc

''' if plots should only be saved, not shown '''
# matplotlib.use('TkAgg')

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
# C-related imports
cimport numpy as np
cimport cython
# from cython.parallel import parallel, prange
from libc.math cimport round
# from libc.math cimport cos

# MAIN
@cython.cdivision(True)    # enables c-division without 1/0 checking
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def main(str init_phases, str init_hist, str int_frequen, int N, double Tsim, double Kstart, double Kend, double fc, double tau, int downscale,
        double m_init_phase, double var_init_phi, double m_int_freqs, double v_int_freqs, double sync_freq, int SamplesPerPeriod, int adiabatic,
        double ad_rate, int specified_distribution, int plot_phi_freq, int wrap_hist_tau, int custom_Kvalues_vector, Kvalues_vector, int plot_out):
    # for simulations with equal distributions of the intrinsic frequencies choose a value for 'specified_distribution'
    if specified_distribution != 0:
        np.random.seed(specified_distribution);
        print('\nFixed distribution for intrinsic frequencies!')
    else:
        np.random.seed();
    #print('fixed random seed for debugging purposes!')

    #check whether there is a file for saving results of sweeps
    save_file=Path('results/sweep_results.txt')
    if (adiabatic == 0 or adiabatic == 2 and not save_file.is_file()):
        print('sweep-results.txt file created and header written.')
        sweep_results=open('results/sweep_results.txt','w+');
        sweep_results.write('Kstart, Fc, tau, lastR, meanR\n');
        sweep_results.close();

    # since c-division is turned on, seperately check for the variables that can be in the denominator
    if Kstart==0 or Kend==0 or N==0 or Tsim==0 or m_int_freqs==0 or downscale==0:
        print('c-division safety check not passed! Check variables.')
        sys.exit()
    # calculate delay steps and time-step ################################################################
    cdef int down = downscale;          # factor to downscale samplerate
    cdef unsigned int tausteps = 0;     # delay in simulation steps
    cdef double dt = 0;                 # time-step according to given sample frequency
    if tau == 0:
        dt = 1.0 / ( SamplesPerPeriod * m_int_freqs );
        print('Simulating at time-step dt  =', dt, '[s] for N = ', N, ', Kstart =', Kstart,
        '[Hz], Kend = ', Kend, '[Hz], fc = ', fc,'[Hz]')
    elif (tau > 0 and sync_freq != 0):
        dt = 1.0 / ( SamplesPerPeriod * sync_freq );
        print('Simulating a delay-case: tau =', tau, '[s] at time-step dt  =', dt,
        '[s] for N = ', N, ', Kstart =', Kstart,'[Hz], Kend = ', Kend, '[Hz], fc = ', fc,'[Hz]')
        tausteps = int(round(tau/dt));
    else:
        print('\nProblem with delay step setup!\n')

    # vector with intrinsic frequencies of the PLLs
    cdef np.ndarray[np.float64_t, ndim=1] int_freqs = np.zeros(N, dtype=np.float64);
    if   int_frequen == 'const':
        int_freqs[:] = m_int_freqs;
    elif int_frequen == 'gaussian':
        int_freqs[:] = np.random.normal(loc=m_int_freqs, scale=np.sqrt(v_int_freqs), size=N);
    elif int_frequen == 'lorentz':
        int_freqs[:] = cauchy.rvs(loc=m_int_freqs, scale=v_int_freqs, size=N);
    else:
        print('\nPlease specify one of the valid options from [const, gaussian, lorentz]')
    now = datetime.datetime.now()
    np.savez('results/intFreqs_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), intF=int_freqs);

    cdef double Krate = ad_rate / Tsim;                   # the rate at which K(t) changes adiabatically
    cdef double Treve = 0;                                # after relaxation, the time until K(t) growth reverses
    if Krate > 0:
        Treve = (Kend - Kstart) / Krate;
    # time-steps and container for histories, coupling strengths and phase time-series
    cdef long int trelsteps = int(round(Tsim/dt));        # relaxation time in steps dt
    cdef long int trevsteps = 0;                          # define time of adiabatic change in steps dt
    cdef long int extrsteps = 0;                          # define time of adiabatic change in steps dt
    if adiabatic == 1:
        trevsteps = int(round(2.0*Treve/dt));             # time of adiabatic change in steps dt
        extrsteps = int(round(0.5*Tsim/dt));              # time of adiabatic change in steps dt
    elif adiabatic == 2:
        extrsteps = int(round(1.0*Tsim/dt));              # time of adiabatic change in steps dt
    cdef long int itersteps = trelsteps+trevsteps+extrsteps;  # total simulation time in dt
    cdef np.ndarray[np.float64_t, ndim=2] phases;
    if plot_phi_freq == 1:
        phases = np.zeros([N, tausteps+itersteps], dtype=np.float64);
        print('Phases vectors of full simulation length, len(phases[0,:])=',len(phases[0,:]))
    elif plot_phi_freq == 0:
        if wrap_hist_tau == 1:                             # evolveOrderReduced will be called
            phases = np.zeros([N, tausteps], dtype=np.float64);
            print('Phases vectors length only tausteps, len(phases[0,:])=',len(phases[0,:]))
        elif wrap_hist_tau == 0:
            phases = np.zeros([N, tausteps+itersteps], dtype=np.float64);
            print('Phases vectors of full simulation length, len(phases[0,:])=',len(phases[0,:]))

    cdef np.ndarray[np.float64_t, ndim=1] K = np.zeros(tausteps+itersteps, dtype=np.float64);
    if adiabatic == 1:
        print('Adiabatic change with Krate =', Krate,'[1/s] for 2xTreve =', 2*Treve,'[s]')
    elif adiabatic == 2:
        print('Stepwise change while keeping the state preserved.')
    else:
        print('No adiabatic change, sweep mode or constant K case.')
    print('Total simulation time is T =', (tausteps+itersteps)*dt,'[s] with Trelax =', Tsim ,'[s]')

    params = {'N': N, 'Kstart': Kstart, 'Kend': Kend, 'tau': tau, 'fc': fc, 'm_int_freqs': m_int_freqs, 'v_int_freqs': v_int_freqs, 'ad_rate': ad_rate,
                ' m_init_phase':  m_init_phase, 'var_init_phi': var_init_phi, 'SamplesPerPeriod': SamplesPerPeriod, 'down': downscale, 'Tsim': Tsim,}
    np.savez('results/params_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), params=params);

    # carry out simulation #############################################################################################################
    ####################################################################################################################################
    cdef double deltaK = 0.025;
    cdef int numberK = int((Kend-Kstart)/deltaK);
    cdef np.ndarray[np.float64_t, ndim=1] Kvalues;
    if ( custom_Kvalues_vector == 0 and adiabatic == 2 ):
        Kvalues = np.linspace(Kstart,Kend,numberK, dtype=np.float64);
    elif ( custom_Kvalues_vector == 1 and adiabatic == 2 ):
        print('NOTE: working with custom values vector for Kvalues! See line 160 in the code.')
        #Kvalues = np.array([0.08, 0.11, 0.14, 0.15, 0.1525, 0.155, 0.1575, 0.16, 0.165, 0.17, 0.18, 0.2,
        #                            0.23, 0.26, 0.29, 0.32, 0.35, 0.38, 0.41, 0.44, 0.47, 0.5]);
        Kvalues = Kvalues_vector;
    print('Kvalues:', Kvalues);
    cdef int ii;
    t0 = time.time()
    inidata = init(init_phases, init_hist, tausteps, m_init_phase, var_init_phi, specified_distribution,
                    dt, N, sync_freq, Kend, Kstart, int_freqs, phases, K)
    # print('time needed for preparation of simulations: ', (time.time()-t0), ' seconds')
    if plot_phi_freq == 1:
        print('Called evolvePhases.')
        simdata = evolvePhases(Kstart, Krate, tausteps, extrsteps, trelsteps, trevsteps, adiabatic, down,
                    N, dt, fc, Tsim, itersteps, int_freqs, inidata['phases'], inidata['K'])
    elif (plot_phi_freq == 0 and wrap_hist_tau == 0 and adiabatic != 2):
        print('Called evolveOrder.')
        simdata = evolveOrder(Kstart, Krate, tausteps, extrsteps, trelsteps, trevsteps, adiabatic, down,
                    N, dt, fc, Tsim, itersteps, int_freqs, inidata['phases'], inidata['K'])
    elif (plot_phi_freq == 0 and wrap_hist_tau == 1 and adiabatic != 2):
        print('Called evolveOrderReduced.')
        simdata = evolveOrderReduced(Kstart, Krate, tausteps, extrsteps, trelsteps, trevsteps, adiabatic, down,
                    N, dt, fc, Tsim, itersteps, int_freqs, inidata['phases'], inidata['K'])
    elif (plot_phi_freq == 0 and wrap_hist_tau == 1 and adiabatic == 2):
        print('Called evolveOrderReduced_discreteChangeK.')
        # run initial case which will provide the new phase history for the next step K + deltaK
        phases = np.array(inidata['phases']);
        #print('last and first phases (initial): ', phases[0,0], phases[0,tausteps])
        # run the rest of the cases
        for ii in xrange(0,len(Kvalues)):
            simdata = evolveOrderReduced_discreteChangeK(Kvalues[ii], tausteps, extrsteps, trelsteps, trevsteps, adiabatic, down,
                        N, dt, fc, Tsim, itersteps, int_freqs, phases)
            Rt     = np.array(simdata['Rt']);
            phases = np.array(simdata['phases']);
            #print('last and first phases (in loop',ii,'): ', phases[0,0], phases[0,tausteps])
            # save result in file
            sweep_results=open('results/sweep_results.txt','a+');
            sweep_results.write('{0}, {1}, {2}, {3}, {4}\n'.format(Kvalues[ii], fc, tau, Rt[len(Rt)-1], np.mean(Rt[-10*int(1/(m_int_freqs*dt)):])));
            sweep_results.close();
        print('Now going backwards and decreasing K.')
        for ii in xrange(len(Kvalues)-2,-1,-1):
            simdata = evolveOrderReduced_discreteChangeK(Kvalues[ii], tausteps, extrsteps, trelsteps, trevsteps, adiabatic, down,
                        N, dt, fc, Tsim, itersteps, int_freqs, phases)
            Rt     = np.array(simdata['Rt']);
            phases = np.array(simdata['phases']);
            #print('last and first phases (in loop',ii,'): ', phases[0,0], phases[0,tausteps])
            # save result in file
            sweep_results=open('results/sweep_results.txt','a+');
            sweep_results.write('{0}, {1}, {2}, {3}, {4}\n'.format(Kvalues[ii], fc, tau, Rt[len(Rt)-1], np.mean(Rt[-10*int(1/(m_int_freqs*dt)):])));
            sweep_results.close();
        # clean up, quit program
        print('\nFinished, results saved in sweep_results.txt')
        exit()

    print('time needed for preparation and execution of simulations: ', (time.time()-t0), ' seconds')
    ####################################################################################################################################
    ####################################################################################################################################

    # delete unnecessary data and cleanup
    del inidata; gc.collect(); print('inidata deleted');
    # extract the system variables from the simulation data and save
    K = simdata['K'];                   # this is already reduced by downscale!
    np.savez('results/Kt_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), Kt=K);
    if plot_phi_freq == 1:
        phases = simdata['phases'];     # this is already reduced by downscale!
        del simdata; gc.collect();
        print('simdata deleted');
        np.savez('results/phases_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), phases=phases);
        # calculate Kuramoto order parameter
        Rt = calcKuramotoOrderParameter(phases);
        np.savez('results/orderparam_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), Rt=Rt);
        # calculate instantaneous frequencies
        instfreq = np.diff(phases, axis=1)/(down*dt);
        if N==3:
            print('\nlast phases:', phases[0:3,-1]%(2.*np.pi), ', phase-differences:', (phases[1,-1]%(2.*np.pi))-(phases[2,-1]%(2.*np.pi)), ', ', (phases[2,-1]%(2.*np.pi))-(phases[3,-1]%(2.*np.pi)));
            # print('{phases[0,-1], phases[1,-1], phases[2,-1]}', phases[3,-1], phases[1,-1], phases[2,-1])
        else:
            del phases; gc.collect()
            print('phases deleted')
    elif plot_phi_freq == 0:
        Rt = simdata['Rt'];             # this is already reduced by downscale!
        np.savez('results/orderparam_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day), Rt=Rt);
        del simdata; gc.collect();
        print('simdata deleted');

    # create time vector
    cdef np.ndarray[np.float64_t, ndim=1] t = np.arange(tausteps+itersteps, dtype=np.float64);
    t = t[::down]; gc.collect()

    # detect when the order parameter crossed the 1/N treshold, save the time, the direction of crossing into a database! ####################
    # average/smooth the order parameter-curve to prevent picking up too many treshold crosses
    treshold_times = [];
    cdef int j;
    for j in range(len(Rt)-1):
        if( Rt[j]>( 1./np.sqrt(float(N)) - 0.1/np.sqrt(float(N)) ) and Rt[j]<( 1./np.sqrt(float(N)) + 0.1/np.sqrt(float(N)) ) ):
            treshold_times.append(t[j])

    # print('treshold_times:', treshold_times)


    # signal = np.random.rand(1000000)
    # th = signal > 0.5
    # th[1:][th[:-1] & th[1:]] = False
    #
    # treshold = Rt[::10] > ( 1./float(N) - 0.01 );

    # plot results
    cdef int i;
    #print('len(t)', len(t), 'len(Rt)', len(Rt), 'len(K)', len(K))
    #print('tausteps', tausteps,'trelsteps', trelsteps,'trevsteps', trevsteps,'extrsteps',extrsteps)
    if ( N == 3 and plot_phi_freq == 1 ):
        plt.figure('phases'); plt.clf();
        for i in xrange(N):
            plt.plot((t*dt), phases[i,:]);
        plt.plot(t[int(tausteps/down)]*dt, 1.0, 'ro', ms=2)
        if adiabatic == 1:
            plt.plot(t[int((tausteps+trelsteps)/down)]*dt, 1.0, 'yo', ms=2)
            plt.plot(t[int((tausteps+trelsteps+trevsteps)/down)]*dt, 1.0, 'co', ms=2)
        plt.xlabel(r'$t$ $[s]$', fontdict = labelfont);
        plt.ylabel(r'$\phi(t)$', fontdict = labelfont);
        # #plt.savefig('phases_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
        # #plt.savefig('phases_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)

    if plot_phi_freq == 1:
        plt.figure('frequencies'); plt.clf();
        for i in xrange(N):
            plt.plot((t[0:len(t)-1])*dt, instfreq[i,:]);
        plt.plot(t[int(tausteps/down)]*dt, 2.0*np.pi*m_int_freqs, 'ro', ms=2)
        if adiabatic == 1:
            plt.plot(t[int((tausteps+trelsteps)/down)]*dt, 2.0*np.pi*m_int_freqs, 'yo', ms=2)
            plt.plot(t[int((tausteps+trelsteps+trevsteps)/down)]*dt, 2.0*np.pi*m_int_freqs, 'co', ms=2)
        plt.xlabel(r'$t$ $[s]$', fontdict = labelfont);
        plt.ylabel(r'$\dot{\phi}(t)$', fontdict = labelfont);
        plt.savefig('results/freqs_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
        plt.savefig('results/freqs_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)
        del instfreq; gc.collect();
        if plot_out == 0:
            plt.close();

    fig, ax1 = plt.subplots();
    ax1.plot(t[int(tausteps/down)]*dt, 0.05, 'ro', ms=2)
    if adiabatic == 1:
        ax1.plot(t[int((tausteps+trelsteps)/down)]*dt, 0.05, 'yo', ms=2)
    else:
        ax1.set_title(r'average order parameter over last 3periods $\bar{R}(-3T_{\omega})=$%.5f' %np.mean(Rt[-10*int(1/(m_int_freqs*dt)):]), fontdict = titlefont);
        sweep_results=open('results/sweep_results.txt','a+');
        sweep_results.write('{0}, {1}, {2}, {3}, {4}\n'.format(Kstart, fc, tau, Rt[len(Rt)-1], np.mean(Rt[-10*int(1/(m_int_freqs*dt)):])));
        sweep_results.close();
    ax1.plot(t*dt,Rt,'b-');
    ax1.plot(np.array(treshold_times)*dt,np.zeros(len(treshold_times))+( 0.1/np.sqrt(float(N)) ),'ko', ms=2)

    ax1.set_xlabel(r'$t$ $[s]$', fontdict = labelfont);
    ax1.set_ylabel(r'$R(t)$', fontdict = labelfont, color='b');
    ax1.tick_params('y', colors='b');
    ax2 = ax1.twinx();
    ax2.plot(t*dt, 2.0*np.pi*K, 'r-');
    ax2.set_ylabel(r'$K(t)$ [radHz]', fontdict = labelfont, color='r'); #y-axis match the line color
    ax2.tick_params('y', colors='r');
    fig.tight_layout();
    plt.savefig('results/order_vs_t_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
    plt.savefig('results/order_vs_t_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)
    if plot_out == 0:
        plt.close();

    if adiabatic == 1:
        plt.figure('orderparameter vs K'); plt.clf();
        #plt.plot(2.0 * np.pi * K, Rt);
        plt.plot(2.0*np.pi*K[0:int(tausteps/down)], Rt[0:int(tausteps/down)],'g--');
        plt.plot(2.0*np.pi*K[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)], Rt[int(tausteps/down):int((trelsteps+0.5*trevsteps)/down)],'r--');
        plt.plot(2.0*np.pi*K[int((tausteps+trelsteps+0.5*trevsteps)/down):], Rt[int((tausteps+trelsteps+0.5*trevsteps)/down):],'b-');
        plt.xlabel(r'$K(t)$ [radHz]', fontdict = labelfont);
        plt.ylabel(r'$R(t)$', fontdict = labelfont);
        plt.savefig('results/order_vs_Kt_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
        plt.savefig('results/order_vs_Kt_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)
        if plot_out == 0:
            plt.close();

    plt.figure('histogram of intrinsic frequencies'); plt.clf();    # plot a histogram of the intrinsic frequencies of the oscillators
    if adiabatic == 1:
        Km = 0.5*( Kend - Kstart ); width=1.0;
    else:
        Km = Kstart; width=1.0;
    plt.hist(int_freqs[:], bins=np.linspace(m_int_freqs-width*abs(Km), m_int_freqs+width*abs(Km), num=51), rwidth=0.75, normed=True);
    # plot Lorentzian
    def test_func(x, a, b):
        return b / (np.pi * ( ( x - a )**2 + b**2 ))
    fs = np.linspace(m_int_freqs-1, m_int_freqs+1, num=10000)
    plt.plot(fs, test_func(fs, m_int_freqs, v_int_freqs),label='Lorentzian')
    plt.legend(loc='best',fontsize=10)
    plt.xlim(m_int_freqs-width*abs(Km), m_int_freqs+width*abs(Km));
    if int_frequen == 'lorentz':
        plt.title(r'Lorentz dist. int. freq. [Hz] $\bar{f}=$%.3f and scale $\bar{\sigma}_f=$%.4f' %( np.mean(int_freqs), v_int_freqs ), fontdict = titlefont);
    elif int_frequen == 'gaussian':
        plt.title(r'intrinsic frequency [Hz] $\bar{f}=$%.3f and std $\bar{\sigma}_f=$%.4f' %( np.mean(int_freqs), np.std(int_freqs) ), fontdict = titlefont);
    elif int_frequen == 'const':
        plt.title(r'const intrinsic frequency [Hz] $\bar{f}=$%.3f' %( np.mean(int_freqs) ), fontdict = titlefont);
    plt.xlabel(r'$\dot{\phi}(0)$ [radHz]', fontdict = labelfont);
    plt.ylabel(r'histogram', fontdict = labelfont);
    plt.savefig('results/dist_init_freqs_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.pdf'.format(Kstart, fc, tau, now.year, now.month, now.day))
    plt.savefig('results/dist_init_freqs_Kstart{0}_Fc{1}_tau{2}_{3}_{4}_{5}.png'.format(Kstart, fc, tau, now.year, now.month, now.day), dpi=300)
    if plot_out == 0:
        plt.close();

    if plot_out == 1:
        plt.draw();
        plt.show();

    return 0

# function to set the initial history for all PLLs ##############################################################
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def init(str init_phases, str init_hist, int tausteps, double m_init_phase, double var_init_phi, int specified_distribution,
         double dt, int N, double sync_freq, double Kend, double Kstart, np.ndarray[np.float64_t, ndim=1] int_freqs,
         np.ndarray[np.float64_t, ndim=2] phases, np.ndarray[np.float64_t, ndim=1] K):
    # for simulations with equal distributions of the intrinsic frequencies choose a value for 'specified_distribution'
    # here we use specified_distribution*0.1, since we want another constant distribution of random numbers for the initial phases
    if specified_distribution != 0:
        np.random.seed(specified_distribution);
        print('\nFixed distribution for initial phases!')
    else:
        np.random.seed();
    # set initial phases according to 'init_phases'
    if init_phases == 'zero':           # set all initial phases to zero
        phases[:,0] = 0;
    elif init_phases == 'gaussian':     # set all initial phases to values drawn from a gaussian distribution
        phases[:,0] = np.random.normal(loc=m_init_phase, scale=np.sqrt(var_init_phi),
                               size=len(phases[:,0]))%(2.*np.pi);
    elif init_phases == 'iid':          # set all initial phases identical distributed in mean +/- var/2
        # print('\n\nphases:', phases, '\n\n');
        phases[:,0] = np.random.uniform(low=m_init_phase-0.5*var_init_phi,
                               high=m_init_phase+0.5*var_init_phi, size=len(phases[:,0]))%(2.*np.pi);
    # write initial history according to 'init_hist' if there is a delay in the system
    cdef int i, j;
    if tausteps > 0:
        if init_hist == 'free':
            K[0:tausteps] = 0;      # if uncoupled system history, coupling strength zero
            for i in xrange(N):
                for j in xrange(tausteps-1):                                        # xrange(tausteps-1) --> 0, ..., tausteps-2
                    phases[i,j+1] = phases[i,j] + 2.0 * np.pi * int_freqs[i] * dt;  # last entry in phases for tausteps-1, since j+1!
        elif init_hist == 'synced':
            K[0:tausteps] = Kstart; # if synchronized system history, coupling strength > 0
            for i in xrange(N):
                for j in xrange(tausteps-1):                                    # xrange(tausteps-1) --> 0, ..., tausteps-2
                    phases[i,j+1] = phases[i,j] + 2.0 * np.pi * sync_freq * dt; # last entry in phases for tausteps-1, since j+1!
                    # print('j+1:',j+1)
        elif init_hist == 'synced_highR':
            K[0:tausteps] = Kend;   # if synchronized system history with large K, coupling strength = Kend
            for i in xrange(N):
                for j in xrange(tausteps-1):                                    # xrange(tausteps-1) --> 0, ..., tausteps-2
                    phases[i,j+1] = phases[i,j] + 2.0 * np.pi * sync_freq * dt; # last entry in phases for tausteps-1, since j+1!
    ''' the initial phase history has been written up to jmax+1, where jmax=tausteps-2, i.e. it is set up to tausteps-1 '''
    # time.sleep(1);
    # print('LAST ENTRY HISTORY: phases[0,len(phases[0,:])-1]', phases[0,len(phases[0,:])-1], ' where len(phases[0,:])-1:', len(phases[0,:])-1);
    # print('last entry: phases[0,tausteps-1]',phases[0,tausteps-1],' next entry: phases[0,tausteps]',phases[0,tausteps],' and phases[0,tausteps+1]',phases[0,tausteps+1]);

    return {'phases': phases, 'K': K}

# function that iterates the all-to-all coupled system using the order parameter #################################
@cython.cdivision(True)    # enables c-division without 1/0 checking
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def evolvePhases(double Kstart, double Krate, int tausteps, int extrsteps, int trelsteps, int trevsteps, int adiabatic, int down,
            int N, double dt, double fc, double Tsim, int itersteps, np.ndarray[np.float64_t, ndim=1] int_freqs,
            np.ndarray[np.float64_t, ndim=2] phases, np.ndarray[np.float64_t, ndim=1] K):
    cdef int i;
    # define temporary variables for control signal, initial frequency and helpers
    cdef np.ndarray[np.float64_t, ndim=1] Xc = np.zeros(N, dtype=np.float64);
    # cdef np.ndarray[np.float64_t, ndim=1] Omeg0 = np.zeros(N, dtype=np.float64);
    cdef double Rx   = 0;
    cdef double Ry   = 0;
    cdef double beta = 2.0*np.pi*fc*dt;   # from numerical integration of the LF dynamics
    # prepare time-dependent coupling strength for adiabatic dynamics
    if adiabatic == 1:                                                          # prepare time-dependent coupling strength for adiabatic case
        K[tausteps:(trelsteps+tausteps+1)] = Kstart;                            # set the coupling strength constant for the relaxation time
        for i in xrange(trevsteps):                                             # set the coupling strength for the times of adiabatic change
            K[i+1+tausteps+trelsteps] = K[i+tausteps+trelsteps] + dt * Krate * np.sign( 0.5 * trevsteps - i );
        K[(tausteps+trelsteps+trevsteps):] = Kstart;                            # set the coupling strength constant after adiabatic change
    else:                                                                       # simulate with constant coupling strength
        K[:] = Kstart;
    # set initial state of the loop-filter
    if tausteps > 0:
        for i in xrange(N):
            # Omeg0[i] = ( phases[i,tausteps-1] - phases[i,tausteps-2] ) / dt;    # calculate initial frequencies
            # Xc[i]    = ( Omeg0[i]/(2.0*np.pi) - int_freqs[i] ) / Kstart;        # calculate initial control signal
            Xc[i]    = ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart;
            # print('initial control signal', ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart)
    else:
        Xc[:] = 0;
    # evolve the system
    cdef double OneOverN = (1.0/N);
    for i in xrange(tausteps-1,(tausteps+itersteps-1)):
        # calculate Rx and Ry
        Rx = OneOverN * np.sum( np.cos(phases[:,i-tausteps]) );
        Ry = OneOverN * np.sum( np.sin(phases[:,i-tausteps]) );
        # calculate control signal
        Xc[:] = (1.0 - beta) * Xc[:] + beta * ( Ry * np.cos(phases[:,i]) - Rx * np.sin(phases[:,i]) );
        # iterate phases
        phases[:,i+1] = phases[:,i] + 2.0 * np.pi * dt * ( int_freqs[:] + K[i] * Xc[:] );
        # print('K[',i,']', K[i])

    print('\nFinal mean (over PLLs) frequency: ', np.mean(phases[:,(tausteps+itersteps-1)] - phases[:,(tausteps+itersteps-1)-1])/dt )
    return {'phases': phases[:,::down], 'K': K[::down]}

# function that iterates the all-to-all coupled system using the order parameter #################################
@cython.cdivision(True)    # enables c-division without 1/0 checking
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def evolveOrder(double Kstart, double Krate, int tausteps, int extrsteps, int trelsteps, int trevsteps, int adiabatic, int down,
            int N, double dt, double fc, double Tsim, int itersteps, np.ndarray[np.float64_t, ndim=1] int_freqs,
            np.ndarray[np.float64_t, ndim=2] phases, np.ndarray[np.float64_t, ndim=1] K):
    cdef long int i;
    # define temporary variables for control signal, initial frequency and helpers
    cdef np.ndarray[np.float64_t, ndim=1] Xc = np.zeros(N, dtype=np.float64);
    cdef np.ndarray[np.float64_t, ndim=1] Rt = np.zeros(tausteps+itersteps, dtype=np.float64);
    # cdef np.ndarray[np.float64_t, ndim=1] Omeg0 = np.zeros(N, dtype=np.float64);
    cdef double Rx   = 0;
    cdef double Ry   = 0;
    cdef double beta = 2.0*np.pi*fc*dt;   # from numerical integration of the LF dynamics
    # prepare time-dependent coupling strength for adiabatic dynamics
    if adiabatic == 1:                                                          # prepare time-dependent coupling strength for adiabatic case
        K[tausteps:(trelsteps+tausteps+1)] = Kstart;                            # set the coupling strength constant for the relaxation time
        for i in xrange(trevsteps):                                             # set the coupling strength for the times of adiabatic change
            K[i+1+tausteps+trelsteps] = K[i+tausteps+trelsteps] + dt * Krate * np.sign( 0.5 * trevsteps - i );
        K[(tausteps+trelsteps+trevsteps):] = Kstart;                            # set the coupling strength constant after adiabatic change
    else:                                                                       # simulate with constant coupling strength
        K[:] = Kstart;
    # set initial state of the loop-filter
    if tausteps > 0:
        for i in xrange(N):
            # Omeg0[i] = ( phases[i,tausteps-1] - phases[i,tausteps-2] ) / dt;    # calculate initial frequencies
            # Xc[i]    = ( Omeg0[i]/(2.0*np.pi) - int_freqs[i] ) / Kstart;        # calculate initial control signal
            Xc[i]    = ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart;
            # print('initial control signal', ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart)
    else:
        Xc[:] = 0;
    # evolve the system
    cdef double OneOverN = (1.0/N);
    for i in xrange(tausteps-1,(tausteps+itersteps-1)):
        # calculate Rx and Ry
        Rx = OneOverN * np.sum( np.cos(phases[:,i-tausteps]) );
        Ry = OneOverN * np.sum( np.sin(phases[:,i-tausteps]) );
        # calculate control signal
        # print('i:', i, 'phases[0,i]', phases[0,i], 'phases[0,i+1]',phases[0,i+1]);
        # time.sleep(1);
        Xc[:] = (1.0 - beta) * Xc[:] + beta * ( Ry * np.cos(phases[:,i]) - Rx * np.sin(phases[:,i]) );
        # iterate phases
        phases[:,i+1] = phases[:,i] + 2.0 * np.pi * dt * ( int_freqs[:] + K[i] * Xc[:] );
        # print('K[',i,']', K[i])
        Rt[i-tausteps] = np.sqrt( Rx**2.0 + Ry**2.0 );

    # print final mean frequency of PLL system
    print('\nFinal mean (over PLLs) frequency: ', np.mean(phases[:,(tausteps+itersteps-1)] - phases[:,(tausteps+itersteps-2)])/dt )
    # calculate Rt for the missing time tau at the end
    # print(len(Rt[-tausteps:]), len(calcKuramotoOrderParameter(phases[:,-tausteps:])))
    Rt[-tausteps-1:] = calcKuramotoOrderParameter(phases[:,-tausteps-1:]);
    del phases; gc.collect();

    # print('Rt[0]:',Rt[0],'Rt[tausteps]:',Rt[tausteps],'Rt[tausteps+trelsteps]:',
    # Rt[tausteps+trelsteps],'Rt[tausteps+trelsteps+int(0.5*trevsteps)]:',Rt[tausteps+trelsteps+int(0.5*trevsteps)],)
    # print('K[0]:',K[0],'K[tausteps]:',K[tausteps],'K[tausteps+trelsteps]:',
    # K[tausteps+trelsteps],'K[tausteps+trelsteps+int(0.5*trevsteps)]:',K[tausteps+trelsteps+int(0.5*trevsteps)],)

    return {'Rt': Rt[::down], 'K': K[::down]}

# function that iterates the all-to-all coupled system using the order parameter #################################
@cython.cdivision(True)    # enables c-division without 1/0 checking
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def evolveOrderReduced(double Kstart, double Krate, int tausteps, int extrsteps, int trelsteps, int trevsteps, int adiabatic,
            int down, int N, double dt, double fc, double Tsim, int itersteps, np.ndarray[np.float64_t, ndim=1] int_freqs,
            np.ndarray[np.float64_t, ndim=2] phases, np.ndarray[np.float64_t, ndim=1] K):

    cdef long int i;
    # define temporary variables for control signal, initial frequency and helpers
    cdef np.ndarray[np.float64_t, ndim=1] Xc = np.zeros(N, dtype=np.float64);
    #print('TEST0, tausteps:',tausteps,'itersteps:',itersteps)
    cdef np.ndarray[np.float64_t, ndim=1] Rt = np.zeros((tausteps+itersteps), dtype=np.float64);
    #print('TEST1')
    # cdef np.ndarray[np.float64_t, ndim=1] Omeg0 = np.zeros(N, dtype=np.float64);
    cdef double Rx   = 0;
    cdef double Ry   = 0;
    cdef double beta = 2.0*np.pi*fc*dt;   # from numerical integration of the LF dynamics
    # prepare time-dependent coupling strength for adiabatic dynamics
    if adiabatic == 1:                                                          # prepare time-dependent coupling strength for adiabatic case
        K[tausteps:(trelsteps+tausteps+1)] = Kstart;                            # set the coupling strength constant for the relaxation time
        for i in xrange(trevsteps):                                             # set the coupling strength for the times of adiabatic change
            K[i+1+tausteps+trelsteps] = K[i+tausteps+trelsteps] + dt * Krate * np.sign( 0.5 * trevsteps - i );
        K[(tausteps+trelsteps+trevsteps):] = Kstart;                            # set the coupling strength constant after adiabatic change
    elif adiabatic == 0:                                                        # simulate with constant coupling strength
        K[:] = Kstart;
    else:
        print('\nERROR!')
    # plt.figure('K(t)'); plt.clf(); plt.plot(K);
    # plt.ylabel(r'$K(t)$ [radHz]', fontdict = labelfont);
    # plt.xlabel(r'$steps$', fontdict = labelfont);
    # plt.savefig('Kt_Kstart{0}_Fc{1}.png'.format(Kstart, fc), dpi=300)
    # plt.draw();
    # plt.show();
    # set initial state of the loop-filter
    if tausteps > 0:
        for i in xrange(N):
            # Omeg0[i] = ( phases[i,tausteps-1] - phases[i,tausteps-2] ) / dt;    # calculate initial frequencies
            # Xc[i]    = ( Omeg0[i]/(2.0*np.pi) - int_freqs[i] ) / Kstart;        # calculate initial control signal
            Xc[i]    = ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart;
            # print('initial control signal', ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / Kstart)
    else:
        Xc[:] = 0;
    # print('len(phases[0,:]):', len(phases[0,:]),' for tausteps:', tausteps)
    #print('IN evolveOrderReduced: phases[0,tausteps-1]:', phases[0,tausteps-1],' phases[0,tausteps]:', phases[0,tausteps],
    #        ' phases[0,tausteps+1]:', phases[0,tausteps+1],'phases[0,tausteps-2]:', phases[0,tausteps-2],
    #        'phases[0,len(phases[0,:])-1]',phases[0,len(phases[0,:])-1],'phases[0,len(phases[0,:])]',phases[0,len(phases[0,:])])
    # evolve the system
    cdef double OneOverN = (1.0/N);
    for i in xrange(tausteps-1,(tausteps+itersteps-1)):
        # print('i%tausteps=',i%tausteps,'(i+1)%tausteps',(i+1)%tausteps)
        # print('phases[0,i%tausteps=',i%tausteps,']:', phases[0,i%tausteps]);#,'phases[0,i+1]:', phases[0,i+1],'phases[0,i-1]:', phases[0,i-1])
        # print('phases[0,(i+1)%tausteps=',(i+1)%tausteps,']:', phases[0,(i+1)%tausteps]);
        # print('i:',i,' (i+1-tausteps)%tausteps:',(i+1-tausteps)%tausteps,' (i+1)%tausteps:',(i+1)%tausteps,' i%tausteps:',i%tausteps)
        # time.sleep(1)
        # calculate Rx and Ry
        Rx = OneOverN * np.sum( np.cos(phases[:,(i+1-tausteps)%tausteps]) );    # (i+1-tausteps)%tausteps: "-tausteps" may no be necessary
        Ry = OneOverN * np.sum( np.sin(phases[:,(i+1-tausteps)%tausteps]) );
        # save order parameter (magnitude)
        Rt[i+1-tausteps] = np.sqrt( Rx**2.0 + Ry**2.0 );
        # calculate control signal
        Xc[:] = (1.0 - beta) * Xc[:] + beta * ( Ry * np.cos(phases[:,i%tausteps]) - Rx * np.sin(phases[:,i%tausteps]) );
        # iterate phases
        #print('phases[0,',(i)%tausteps,']',phases[0,(i)%tausteps],'before update phases[0,',(i+1)%tausteps,']',phases[0,(i+1)%tausteps])
        phases[:,(i+1)%tausteps] = phases[:,i%tausteps] + 2.0*np.pi*dt* ( int_freqs[:] + K[i] * Xc[:] );
        #print('phases[0,',(i)%tausteps,']',phases[0,(i)%tausteps],'after update phases[0,',(i+1)%tausteps,']',phases[0,(i+1)%tausteps])
        # print('\nRx: ', Rx, '  Ry: ', Ry, 'Xc[0]: ', Xc[0], ' 2.0*np.pi*dt* ( int_freqs[:] + K[i] * Xc[:]:', 2.0*np.pi*dt* ( int_freqs[:] + K[i] * Xc[:] ) )
        #time.sleep(1)
        # print('K[',i,']', K[i])
        # print('mean (over PLLs) instantaneous frequency: ', np.mean(phases[:,(i+1)%tausteps] - phases[:,i%tausteps])/dt)

    # calculate last i of above loop, necessary to know where to start below when calculating the remaining order parameter
    cdef int last_i = (tausteps+itersteps-1)%tausteps;
    #print('last i-value:', last_i)
    # print final mean frequency of PLL system
    print('\nFinal mean (over PLLs) frequency: ', np.mean(phases[:,last_i] - phases[:,(last_i-1)])/dt )
    # calculate Rt for the missing time tau at the end; in the loop Rt has been filled up to index itersteps
    for i in xrange(itersteps,itersteps+tausteps):
        Rt[i] = np.abs(np.mean(np.exp(1j * phases[:,(i+last_i)%tausteps]), axis=0));
    del phases; gc.collect();

    return {'Rt': Rt[::down], 'K': K[::down]}

# function that iterates the all-to-all coupled system using the order parameter #################################
@cython.cdivision(True)    # enables c-division without 1/0 checking
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def evolveOrderReduced_discreteChangeK(double K, int tausteps, int extrsteps, int trelsteps, int trevsteps, int adiabatic,
                                int down, int N, double dt, double fc, double Tsim, int itersteps, np.ndarray[np.float64_t, ndim=1] int_freqs,
                                np.ndarray[np.float64_t, ndim=2] phases):

    #print('last and first phases (in function): ', phases[0,0], phases[0,tausteps])
    cdef long int i;
    # define temporary variables for control signal, initial frequency and helpers
    cdef np.ndarray[np.float64_t, ndim=1] Xc = np.zeros(N, dtype=np.float64);
    #print('TEST0, tausteps:',tausteps,'itersteps:',itersteps)
    cdef np.ndarray[np.float64_t, ndim=1] Rt = np.zeros((tausteps+itersteps), dtype=np.float64);
    #print('TEST1')
    # cdef np.ndarray[np.float64_t, ndim=1] Omeg0 = np.zeros(N, dtype=np.float64);
    cdef double Rx   = 0;
    cdef double Ry   = 0;
    cdef double beta = 2.0*np.pi*fc*dt;   # from numerical integration of the LF dynamics
    # plt.figure('K(t)'); plt.clf(); plt.plot(K);
    # plt.ylabel(r'$K(t)$ [radHz]', fontdict = labelfont);
    # plt.xlabel(r'$steps$', fontdict = labelfont);
    # plt.savefig('Kt_K{0}_Fc{1}.png'.format(K, fc), dpi=300)
    # plt.draw();
    # plt.show();
    # set initial state of the loop-filter
    if tausteps > 0:
        for i in xrange(N):
            # Omeg0[i] = ( phases[i,tausteps-1] - phases[i,tausteps-2] ) / dt;    # calculate initial frequencies
            # Xc[i]    = ( Omeg0[i]/(2.0*np.pi) - int_freqs[i] ) / K;             # calculate initial control signal
            Xc[i]    = ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / K;
            # print('initial control signal', ( (phases[i,tausteps-1] - phases[i,tausteps-2]) / (dt*2.0*np.pi) - int_freqs[i] ) / K)
    else:
        Xc[:] = 0;
    # print('len(phases[0,:]):', len(phases[0,:]),' for tausteps:', tausteps)
    #print('IN evolveOrderReduced: phases[0,tausteps-1]:', phases[0,tausteps-1],' phases[0,tausteps]:', phases[0,tausteps],
    #        ' phases[0,tausteps+1]:', phases[0,tausteps+1],'phases[0,tausteps-2]:', phases[0,tausteps-2],
    #        'phases[0,len(phases[0,:])-1]',phases[0,len(phases[0,:])-1],'phases[0,len(phases[0,:])]',phases[0,len(phases[0,:])])
    # evolve the system
    cdef double OneOverN = (1.0/N);
    for i in xrange(tausteps-1,(tausteps+itersteps-1)):
        # print('i%tausteps=',i%tausteps,'(i+1)%tausteps',(i+1)%tausteps)
        # print('phases[0,i%tausteps=',i%tausteps,']:', phases[0,i%tausteps]);#,'phases[0,i+1]:', phases[0,i+1],'phases[0,i-1]:', phases[0,i-1])
        # print('phases[0,(i+1)%tausteps=',(i+1)%tausteps,']:', phases[0,(i+1)%tausteps]);
        # print('i:',i,' (i+1-tausteps)%tausteps:',(i+1-tausteps)%tausteps,' (i+1)%tausteps:',(i+1)%tausteps,' i%tausteps:',i%tausteps)
        # time.sleep(1)
        # calculate Rx and Ry
        Rx = OneOverN * np.sum( np.cos(phases[:,(i+1-tausteps)%tausteps]) );    # (i+1-tausteps)%tausteps: "-tausteps" may no be necessary
        Ry = OneOverN * np.sum( np.sin(phases[:,(i+1-tausteps)%tausteps]) );
        # save order parameter (magnitude)
        Rt[i+1-tausteps] = np.sqrt( Rx**2.0 + Ry**2.0 );
        # calculate control signal
        Xc[:] = (1.0 - beta) * Xc[:] + beta * ( Ry * np.cos(phases[:,i%tausteps]) - Rx * np.sin(phases[:,i%tausteps]) );
        # iterate phases
        #print('phases[0,',(i)%tausteps,']',phases[0,(i)%tausteps],'before update phases[0,',(i+1)%tausteps,']',phases[0,(i+1)%tausteps])
        phases[:,(i+1)%tausteps] = phases[:,i%tausteps] + 2.0*np.pi*dt* ( int_freqs[:] + K * Xc[:] );
        #print('phases[0,',(i)%tausteps,']',phases[0,(i)%tausteps],'after update phases[0,',(i+1)%tausteps,']',phases[0,(i+1)%tausteps])
        # print('\nRx: ', Rx, '  Ry: ', Ry, 'Xc[0]: ', Xc[0], ' 2.0*np.pi*dt* ( int_freqs[:] + K[i] * Xc[:]:', 2.0*np.pi*dt* ( int_freqs[:] + K[i] * Xc[:] ) )
        #time.sleep(1)
        # print('K[',i,']', K[i])
        # print('mean (over PLLs) instantaneous frequency: ', np.mean(phases[:,(i+1)%tausteps] - phases[:,i%tausteps])/dt)

    # calculate last i of above loop, necessary to know where to start below when calculating the remaining order parameter
    cdef int last_i = (tausteps+itersteps-1)%tausteps;
    #print('last i-value:', last_i)
    # print final mean frequency of PLL system
    print('\nFor K =', K, ' radHz final mean (over PLLs) frequency: ', np.mean(phases[:,last_i] - phases[:,(last_i-1)])/dt )
    # calculate Rt for the missing time tau at the end; in the loop Rt has been filled up to index itersteps
    for i in xrange(itersteps,itersteps+tausteps):
        Rt[i] = np.abs(np.mean(np.exp(1j * phases[:,(i+last_i)%tausteps]), axis=0));

    # rearrange phases such that last calculated values are first vector entries via numpy.roll()
    return {'Rt': Rt[::down], 'K': K, 'phases': np.roll(phases, -last_i, axis=1)}

#function to calculate the Kuramoto order parameter #############################################################
def calcKuramotoOrderParameter(np.ndarray[np.float64_t, ndim=2] phi):
    '''Computes the Kuramoto order parameter r for in-phase synchronized states

       Parameters
       ----------
       phi:  np.array
            real-valued 2d matrix or 1d vector of phases
            in the 2d case the columns of the matrix represent the individual oscillators

       Returns
       -------
       r  :  np.array
            real value or real-valued 1d vector of the Kuramotot order parameter

       Authors
       -------
       Lucas Wetzel, Daniel Platz'''
    # Complex phasor representation
    z = np.exp(1j * phi);
    # Kuramoto order parameter
    r = np.abs(np.mean(z, axis=0));

    return r

#################################################################################################################
# t1 = time.time()
# main()
# print('time needed for entire program: ', (time.time()-t1), ' seconds')

#################################################################################################################

    # global parameters, number of PLLs, coupling strengths, cut-off frequency, transmission deĺay, ##############
    # mean intrinsic frequency, and scale parameter of intrinsic frequency
#    cdef str    init_phases  = 'iid'       # specifies initial phases of the PLLs [zero, gaussian, iid]
#    cdef str    init_hist    = 'free'      # specifies initial history for simulation [free, synced]
#    cdef str    int_frequen  = 'lorentz'   # specifies dist. of intrinsic frequencies [const, gaussian, lorentz]
#    cdef int    N            = 100;        # number of PLLs in the system
#    cdef double Tsim         = 750;        # simulation time in seconds, or Trelax the relaxation time (adiabtic case)
#    cdef double Kstart       = 0.002;      # starting coupling strength in Hz(/V) THIS IS NOT Kvco, K=Kvco/2
#    cdef double Kend         = 0.8;        # maximum coupling strength in Hz(/V)
#    cdef double fc           = 0.159;      # cut-off frequency in Hz
#    cdef double tau          = 0.0;        # transmission delay in Hz
#    cdef double m_init_phase = 3.14;       # mean initial phases of PLLs
#    cdef double var_init_phi = 6.28;       # variance of gaussian or iid distribution of initial phases
#    cdef double m_int_freqs  = 0.4779;     # mean intrinsic frequency of the PLLs in Hz
#    cdef double v_int_freqs  = 0.0001;     # variance of gaussian or scale of lorentz distributed int. frequencies
#    cdef double sync_freq    = 0.4779;     # initial frequency of the system in case of coupled initial state
#    cdef double dt           = 0.01;       # time-step for iteration of Euler-sheme / Runge-Kutta 4th order
