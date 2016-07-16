# -*- coding: utf-8 -*-
''' accpy.measure.tunes
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division


def measure_tunes(f_HF):
    # prepare figures
    fig = figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_xlabel('Frequency / (kHz)')
    ax1.set_ylabel('Amplitude / (dB)')
    ax1.ticklabel_format(useOffset=False)
    ax1.set_xlim([f_cent-f_span/2, f_cent+f_span/2])
    ax1.set_ylim([-90, -10])
    ax2 = fig.add_subplot(1, 2, 2)
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
    canvas = FigureCanvasTkAgg(fig, master=tabs[3])
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    fig.canvas.draw()
    t = linspace(5, 90, steps)  # time in ms
    if mode == 1:
        filename = filestr.get()
        fftdata = load(filename, fileinfo, 'fftdata')[0]
        points = shape(fftdata)[2]
    elif mode == 2:
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
            if mode == 2:
                fftdata[i, j, :] = get_fftdata(instr, t[i], f_span, f_cent)
                sleep(1.2)
            ln1_fft[0].set_data(f, fftdata[i, j, :])
            init_mu[j, :] = f[findpeaks(f, fftdata[i, j, :], f_cent, f_span, points, fig, ax1, ln1_peak, ln1_thresh)]
            # fig.canvas.draw()
            updatefig(fig, ax1, [ln1_fft])
            updatefig(fig, ax1, [ln1_peak])
            sleep(.2)
        ln1_fft[0].set_data([], [])
        fft_mean = nanmean(fftdata[i, :, :], 0)
        fft_error = nanstd(fftdata[i, :, :], 0)
        setdata_errorbary(ln1_avg, f, fft_mean, fft_error)
        findpeaks(f, fft_mean, f_cent, f_span, points, fig, ax1, ln1_peak, ln1_thresh)
        # fig.canvas.draw()
        updatefig(fig, ax1, [ln1_fft, ln1_avg, ln1_peak])

        mu[:, i] = nanmean(init_mu, 0)
        sigma[:, i] = nanstd(init_mu, 0)

        fs[0, i] = (mu[2, i] - mu[0, i])/2
        tunestr[2].set('f_longitudinal={:g} kHz'.format(fs[0, i]))
        fs[1, i] = sqrt((sigma[0, i]/2)**2 + (sigma[2, i]/2)**2)
        setdata_errorbary(ln2_fs, t[:i+1], fs[0, :i+1], fs[1, :i+1])
        # fig.canvas.draw()
        updatefig(fig, ax2, [ln2_fs])
        sleep(.2)
    if savedata:
        save('SynchroTune', fileinfo, fftdata=fftdata)
    if savefigs:
        figsave(fig, 'SynchroTune', filetype=filetype)
    if mode == 2:
        instr.write('TIM:HOR:POS {}'.format(20/1e3))
    return