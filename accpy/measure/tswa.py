# -*- coding: utf-8 -*-
''' accpy.measure.tswa
author:
    Felix Kramer
'''
from __future__ import division, print_function
from numpy import (linspace, arange, concatenate, pi, complex128, zeros,
                   empty, angle, unwrap, diff, abs as npabs, ones, vstack, exp,
                   sin, log, linalg, polyfit)
import scipy.optimize as scop
import pyfftw
try:
    import cPickle as pickle
except ImportError:
    import pickle
from ..math.specialfuns import hanning

n = pyfftw.simd_alignment
effort = ['FFTW_ESTIMATE', 'FFTW_MEASURE', 'FFTW_PATIENT', 'FFTW_EXHAUSTIVE']

def init_pyfftw(x, effort=effort[0], wis=False):
    N = len(x)
    a = pyfftw.n_byte_align_empty(int(N), n, 'complex128')
    a[:] = x
    if wis is not False:
        pyfftw.import_wisdom(wis)
        fft = pyfftw.builders.fft(a, threads=8)
        ifft = pyfftw.builders.ifft(a, threads=8)
    else:
        fft = pyfftw.builders.fft(a, planner_effort=effort, threads=8)
        ifft = pyfftw.builders.ifft(a, planner_effort=effort, threads=8)
    return fft, ifft



def evaltswa(counts, bunchcurrents, effort='FFTW_MEASURE'):
    fs = 1.2495*1e6  # sampling frequency in Hz
    dt = 1/fs
    t = arange(len(counts))*dt*1e3  # time in ms
    bbfbcntsnorm = (counts.T/bunchcurrents).T

    N = len(bbfbcntsnorm[0, :])
    n = len(bbfbcntsnorm[:, 0])
    with open('calib_InterpolatedUnivariateSpline.pkl', 'rb') as fh:
        calib = pickle.load(fh)
    bbfbpos = [calib(bbfbcntsnorm[i, :]) for i in range(n)]

    init_pyfftw(bbfbpos[0], effort=effort)
    wisdom = pyfftw.export_wisdom()

    # Turn on the cache for optimum pyfftw performance
    pyfftw.interfaces.cache.enable()

    # create frequency vector
    fd = linspace(0, fs/2/1e3, N/2)

    # prepare frequency filter
    fcent, fsigm = 190, 50
    fleft, fright = fcent - fsigm, fcent + fsigm
    pts_lft = sum(fd < fleft )
    pts_rgt = sum(fd > fright)
    pts_roi = len(fd) - pts_lft - pts_rgt
    frequencyfilter = concatenate((zeros(pts_lft), hanning(pts_roi), zeros(pts_rgt+N/2)))

    # predefine lists
    fftx = empty([n, N], dtype=complex128)
    fftx_clipped = empty([n, N], dtype=complex128)
    fftx_filtered = empty([n, N], dtype=complex128)
    analytic_signal = empty([n, N], dtype=complex128)
    amplitude_envelope = empty([n, N-1])
    instantaneous_phase = empty([n, N])
    instantaneous_frequency = empty([n, N-1])

    for i in range(n): 
        # initialise pyfftw for both signals
        myfftw, myifftw = init_pyfftw(bbfbpos[i], wis=wisdom)

        # calculate fft of signal
        fftx[i, :] = myfftw(bbfbpos[i])

        # clip negative frequencies
        fftx_clipped[i, :] = fftx[i, :]
        fftx_clipped[i, N/2+1:] = 0

        # restore lost energy of negative frequencies
        fftx_clipped[i, 1:N/2] *= 2

        # apply frequency filter
        fftx_filtered[i, :] = fftx_clipped[i, :]*frequencyfilter

        # calculate inverse fft (analytical signal) of filtered and positive frequency only fft
        analytic_signal[i, :] = myifftw(fftx_filtered[i, :])
        amplitude_envelope[i, :] = npabs(analytic_signal[i, :])[:-1]
        instantaneous_phase[i, :] = unwrap(angle(analytic_signal[i, :]))
        instantaneous_frequency[i, :] = diff(instantaneous_phase[i, :]) / (2*pi) * fs

    ''' Damping time
    * amplitude damping time only half of center of mass damping time
    * chromaticity dependant
    * in horizontal plane dispersion and energy dependant
    * dephasing
    * landau damping
    * head tail damping
        > analytic two particle modell shows directional interaction from head to tail and position interchange
    '''

    beg, end = 23, 6000
    t2 = linspace(0, t[-1], N-1)[beg:end]
    amplit = [amplitude_envelope[i, beg:end] for i in range(n)]
    signal = [instantaneous_frequency[i, :][beg:end] for i in range(n)]
    fdamp = []
    initialamp = empty(n)
    for i in range(n):
        '''
        ln[A*e^(dt)] = ln(A) + d*t
        from linear fit: y = m*t + c we gain:
                     A = e^c
                     d = m
        '''
        M = vstack([t2, ones(len(t2))]).T
        tau_inverse, const = linalg.lstsq(M, log(amplit[i]))[0]
        tau_coherent = -1/tau_inverse
        initialamp[i] = exp(const)
        fdamp.append(lambda t, Amplitude=initialamp[i], tau_coherent=tau_coherent: Amplitude*exp(-t/tau_coherent))

    ''' Instantaneous frequency
    * square increase over amplitude
    * frequency is overlayed with synchrotron frequency (~7kHz)
    * filter out synchrotron frequency with a bandstop filter -> tricky (bad snr in fft)
    * fit assumed square funtion -> wrong
    f(amp) = a*amp**2 + b
    amp(t) = fdamp -> tau
    f(t) = a*exp(2*-t/tau) + b
    '''
    dump = 16
    fitfun = lambda t, a, b, c, f, d, tau: a + b*sin(c + t*f) - d*exp(2*-t/tau)
    f_instfreq = []
    for i in range(n-dump):
        popt, pcov = scop.curve_fit(fitfun, t2, signal[i]/1e3, p0=[194.5, 0.6, 1, 44.5, 7, 1.3])
        f_instfreq.append(lambda t, a=popt[0], d= popt[4], tau=popt[5]: a - d*exp(2*-t/tau))
    
    ''' Amolitude dependant tune shift
    '''
    tswa, b = empty(n-dump), empty(n-dump)
    for i in range(n-dump):
        ampl = fdamp[i](t2)
        freq = f_instfreq[i](t2)
        tswa[i], b[i] = polyfit(ampl**2, freq, 1)
        fitfun = lambda t: tswa[i]*t + b[i]
    
    return t, t2, bbfbcntsnorm, amplit, fdamp, signal, f_instfreq