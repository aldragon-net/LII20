import numpy as np
import scipy as sp
import os

def subtract_offset(signal):
    """determine and subtract offset from given signal"""
    offset_t = (signal.shape[-1]*5)//100   #5% of axis considered as offsetzone
    offset = signal[1,:offset_t].mean()    #mean level set as offset 
    noise = signal[1,:offset_t].std()      #STD set as noise level 
    mod_signal = np.copy(signal)
    mod_signal[1] = signal[1] - offset      #offset subtracted
    return mod_signal, offset
    
def flip_signal(signal):
    """flip if necessary tom make the peak positive"""
    s_min = np.amin(signal[1])             #if peak is positive or negative
    s_max = np.amax(signal[1])
    is_signal_flipped = False
    if abs(s_min) > abs(s_max):
        signal[1] = -signal[1]             #flip signal if negative
        is_signal_flipped = True
    return signal, is_signal_flipped    

def SG_smooth_signal(signal):
    """Savitsky-Golay smoothing of signal"""
    signal_points = signal[1]
    smoothed_signal = np.copy(signal)
    smoothed_signal[1] = sp.signal.savgol_filter(signal_points, 51, 3)
    return smoothed_signal
    
def cut_tails(signal):
    """cut zero tails of signal"""
    s_max = np.amax(signal[1])
    threshold = 0.3*s_max
    for i in range(signal.shape[-1]):
        if signal[1,i] > threshold:
            rise_index = i
            break
    left_cut_index = (rise_index*3) // 4
        
    noise = signal[1,:left_cut_index].std()
    
    for i in range(rise_index, signal.shape[-1]):
        if signal[1,i] < noise:
            fall_index = i
            break
    right_cut_index = fall_index
    
    trunc_signal = signal[:, left_cut_index:right_cut_index]
    
    return trunc_signal
    
def process_signals(signal_1, signal_2):
    """process and visualize input signals"""
    shifted_signal_1, offset_1 = subtract_offset(signal_1)
    flipped_signal_1, is_signal_1_flipped = flip_signal(shifted_signal_1)
    smoothed_signal_1 = SG_smooth_signal(flipped_signal_1)
    trunc_signal_1 = cut_tails(smoothed_signal_1)
    offset_line_1 =[[signal_1[0,0]*1e9, signal_1[0,-1]*1e9], [offset_1, offset_1]] 

    signal_2 = read_LIIfile(signal_path_2)
    shifted_signal_2, offset_2 = subtract_offset(signal_2)
    flipped_signal_2, is_signal_2_flipped = flip_signal(shifted_signal_2)
    smoothed_signal_2 = SG_smooth_signal(flipped_signal_2)
    trunc_signal_2 = cut_tails(smoothed_signal_2)
    offset_line_2 =[[signal_2[0,0]*1e9, signal_2[0,-1]*1e9], [offset_2, offset_2]] 

    plt.subplot(231)
    plt.plot(offset_line_1[0], offset_line_1[1], 'k--',
             signal_1[0]*1e9, signal_1[1], 'r-')       
    plt.ylabel('I, V')

    plt.subplot(234)
    plt.plot(offset_line_2[0], offset_line_2[1], 'k--',
             signal_2[0]*1e9, signal_2[1], 'b-')
    plt.ylabel('I, V')
    plt.xlabel('t')

    plt.subplot(232)
    plt.plot([flipped_signal_1[0,0]*1e9, flipped_signal_1[0,-1]*1e9], [0,0], 'k--',
             flipped_signal_1[0]*1e9, flipped_signal_1[1], 'r-',
             smoothed_signal_1[0]*1e9, smoothed_signal_1[1], 'y-')

    plt.subplot(235)
    plt.plot([flipped_signal_2[0,0]*1e9, flipped_signal_2[0,-1]*1e9], [0,0], 'k--',
             flipped_signal_2[0]*1e9, flipped_signal_2[1], 'b-',
             smoothed_signal_2[0]*1e9, smoothed_signal_2[1], 'y-')       
    plt.xlabel('t')

    plt.subplot(233)
    plt.plot([trunc_signal_1[0,0]*1e9, trunc_signal_1[0,-1]*1e9], [0,0], 'k--',
             trunc_signal_1[0]*1e9, trunc_signal_1[1], 'r-')

    plt.subplot(236)
    plt.plot([trunc_signal_2[0,0]*1e9, trunc_signal_2[0,-1]*1e9], [0,0], 'k--',
             trunc_signal_2[0]*1e9, trunc_signal_2[1], 'b-')      
    plt.xlabel('t')

    plt.suptitle('Input LII signals')
    plt.show()
    
    


       
    

   
