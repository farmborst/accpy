# -*- coding: utf-8 -*-
''' accpy.measure.tunes
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (linspace, zeros, nan, sqrt, shape, nanmean, nanstd, array,
                   nanargmax)
from struct import unpack
from time import sleep
from ..dataio.hdf5 import load, save


def get_fftdata(instr, t, f_span, f_cent):  # measurement
    instr.write('TIM:HOR:POS {}'.format(t/1e3))    # set time position to t
    data = instr.ask_raw('CALC:MATH:DATA?')
    try:
        if data[0] == '#':
            n = int(data[1])
            N = int(data[2:2+n])
            fftdata = array(unpack('{}f'.format(N/4), data[2+n:2+n+N]))
    except:
        print('scope reading error')
    return fftdata


def setdata_errorbary(lines, x, y, y_error, empty=False):
    ln, (erry_top, erry_bot), (barsy,) = lines
    ln.set_data(x, y)
    yerr_top = y + y_error
    yerr_bot = y - y_error
    erry_top.set_xdata(x)
    erry_bot.set_xdata(x)
    erry_top.set_ydata(yerr_top)
    erry_bot.set_ydata(yerr_bot)
    if empty:
        new_segments_y = [array([[x, yerr_top], [x, yerr_bot]])]
    else:
        new_segments_y = [array([[xv, yt], [xv,yb]]) for xv, yt, yb in zip(x, yerr_top, yerr_bot)]
    barsy.set_segments(new_segments_y)


def updatefig(fig, ax, lines):
    for line in lines:
        ax.draw_artist(ax.patch)
        if len(line) == 1:
            ax.draw_artist(line[0])
        else:
            ax.draw_artist(line[0])
            ax.draw_artist(line[1][0])
            ax.draw_artist(line[1][1])
            ax.draw_artist(line[2][0])
    #ax.relim()
    #ax.autoscale_view(tight=True, scalex=True, scaley=True)
    ax.grid()
    fig.canvas.blit()
    fig.canvas.flush_events()


def findpeaks(f, fft, f_cent, f_span, points, fig, ax, line1, line2):
    center = round(points/2.)
    region = round(points/8.)
    lc = center - region
    rc = center + region
    region1 = round(points/6.)
    l = lc - region1
    r = rc + region1
    mu1 = nanargmax(fft[l:lc]) + l
    mu2 = nanargmax(fft[lc:rc]) + lc
    mu3 = nanargmax(fft[rc:r]) + rc
    args = [mu1, mu2, mu3]
    line1[0].set_data(f[args], fft[args])
    return args


def measure_tunes(figs, tunestr, mode, filename, f_rf, h, bunch, steps):
    f_span = 200.
    statistic = 6
    ''' RF parameters of BESSY II booster
    Revolution time:        t_rev = 96 m/(3e8 m/s) = 320 ns
    Revolution frequency:   f_rev = 1/t_rev = 3125 kHz
    Harmonic number:        h = 160 = maximum number of bunches
    Cavity frequency:       f_rf = h*f_rev = 160*312.5 = 500 MHz
    Synchrotron frequency:  f_s
    Synchrotron Tune:       Q_s = f_s/f_rev  '''
    f_rev = f_rf/h
    f_cent = f_rev*bunch
    # prepare figures
    fig = figs[2]
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.set_xlabel('Frequency / (kHz)')
    ax1.set_ylabel('Amplitude / (dB)')
    ax1.ticklabel_format(useOffset=False)
    ax1.set_xlim([f_cent-f_span/2, f_cent+f_span/2])
    ax1.set_ylim([-90, -10])
    ax2.set_xlabel('Time / (ms)')
    ax2.set_ylabel('Synchrotron frequency / (kHz)')
    ax2.set_xlabel('Time / (ms)')
    ax2.set_xlim([0, 100])
    ax2.set_ylim([30, 60])
    # fft data
    ln1_fft = ax1.plot([], [], '.k')
    ln1_thresh = ax1.plot([], [], '-b')
    ln1_peak = ax1.plot([], [], 'or', mfc='none', mec='red', ms=20, linewidth=15)
    ln1_avg = ax1.errorbar(nan, nan,yerr=nan, marker='.', mfc='b', ecolor='b', mec='b', ls='None')
    # synchrotron frequency
    ln2_fs = ax2.errorbar(nan, nan, yerr=nan, marker='.', mfc='b', ecolor='b', mec='b', ls='None')
    fig.canvas.draw()
    t = linspace(5, 90, steps)  # time in ms
    if mode == 'From File':
        fftdata = load(filename, False, 'fftdata')[0]
        points = shape(fftdata)[2]
    elif mode == 'Measurement':
        # test scope communication
        try:
            from vxi11 import Instrument
            instr = Instrument('scopez6g.ctl.bessy.de')
            print(instr.ask('*IDN?'))       # SCOPEZ6G/Diagnose (ROHDE&SCHWARZ)
            instr.write('TIM:RANGE 0.001')  # set time range to 1ms
            instr.write('TIM:HOR:POS {}'.format(t[0]/1e3))
            instr.write('CALC:MATH:FFT:CFR {} kHz'.format(f_cent))
            instr.write('CALC:MATH:FFT:SPAN {} kHz'.format(f_span))
            instr.write('FORMAT:DATA REAL,32')
        except:
            print('Communication with scope failed!')
            exit(1)
        points = len(get_fftdata(instr, t[0], f_span, f_cent))
        fftdata = zeros([steps, statistic, points])
    f = linspace(f_cent-f_span/2, f_cent+f_span/2, points)  # f in kHz
    fs = zeros([2, steps])
    mu, sigma = [zeros([3, steps]) for i in range(2)]
    for i in range(steps):
        ln1_peak[0].set_data([], [])
        setdata_errorbary(ln1_avg, nan, nan, nan, empty=True)
        init_mu = zeros((statistic, 3))
        for j in range(statistic):
            if mode == 'Measurement':
                fftdata[i, j, :] = get_fftdata(instr, t[i], f_span, f_cent)
                sleep(1.2)
            ln1_fft[0].set_data(f, fftdata[i, j, :])
            init_mu[j, :] = f[findpeaks(f, fftdata[i, j, :], f_cent, f_span, points, fig, ax1, ln1_peak, ln1_thresh)]
            fig.canvas.draw()
            #updatefig(fig, ax1, [ln1_fft])
            #updatefig(fig, ax1, [ln1_peak])
            sleep(.2)
        ln1_fft[0].set_data([], [])
        fft_mean = nanmean(fftdata[i, :, :], 0)
        fft_error = nanstd(fftdata[i, :, :], 0)
        setdata_errorbary(ln1_avg, f, fft_mean, fft_error)
        findpeaks(f, fft_mean, f_cent, f_span, points, fig, ax1, ln1_peak, ln1_thresh)
        fig.canvas.draw()
        #updatefig(fig, ax1, [ln1_fft, ln1_avg, ln1_peak])

        mu[:, i] = nanmean(init_mu, 0)
        sigma[:, i] = nanstd(init_mu, 0)

        fs[0, i] = (mu[2, i] - mu[0, i])/2
        tunestr[2].set('f_longitudinal={:g} kHz'.format(fs[0, i]))
        fs[1, i] = sqrt((sigma[0, i]/2)**2 + (sigma[2, i]/2)**2)
        setdata_errorbary(ln2_fs, t[:i+1], fs[0, :i+1], fs[1, :i+1])
        fig.canvas.draw()
        #updatefig(fig, ax2, [ln2_fs])
        sleep(.2)
    if mode == 'Measurement':
        instr.write('TIM:HOR:POS {}'.format(20/1e3))
    return figs