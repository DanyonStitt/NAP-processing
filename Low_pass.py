# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 13:27:16 2022

@author: dst74
"""

from scipy.signal import butter, lfilter

def butter_lowpass(cutoff, fs, order=8):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=8):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y