import numpy as np
import os

from tkinter import *
from tkinter import filedialog
from tkinter import font

def ask_for_signals():
    root = Tk()
    root.withdraw()
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file 1",filetypes = (("csv files","*.csv"),("all files","*.*")))
    print ('Signal 1 file: {0}'.format(root.filename))
    signal_path_1 = root.filename
    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file 2",filetypes = (("csv files","*.csv"),("all files","*.*")))
    print ('Signal 2 file: {0}'.format(root.filename))
    signal_path_2 = root.filename     
    root.destroy()
    return signal_path_1, signal_path_2

def read_LIIfile(filepath):
    """read signle file"""
    data = []
    with open(filepath, 'r') as inpfile:
        while True:
            line = inpfile.readline()
            if line.strip() == '' :
                break   #end of file  
            try:    
                data.append([float(line.split(',')[0]), float(line.split(',')[1])])
            except Exception: continue
    inpfile.close()
    signal = np.array(data)
    signal = signal.transpose()
    return signal
    
def read_LIIdir(folderpath):
    """read all files from given folder"""
    data = []
    files = os.listdir(folderpath)
    for file in files:
        filepath = folderpath+'\\'+file
        signal = read_LIIfile(filepath)
        data.append(signal)
    try:
        signals = np.stack(data)
    except Exception:
        return 'files of different length'
    for n in range(signals.shape[0]):
        if (signals[n,0,1]-signals[n,0,0]) - (signals[0,0,1]-signals[0,0,0]) > 1e-15:
            return 'files of different time steps'
    return signals
  
def average_LIIsignals(signals):
    """returns averaged signal"""
    rise_indices = []
    for n in range(signals.shape[0]):
        signal = signals[n,...]
        offset_t = (signal.shape[-1]*5)//100   #5% of axis considered as offsetzone
        offset = signal[1,:offset_t].mean()    #mean level set as offset 
        noise = signal[1,:offset_t].std()      #STD set as noise level 
        
        signal[1] = signal[1] - offset         #offset subtracted
        
        s_min = np.amin(signal[1])             #if peak is positive or negative
        s_max = np.amax(signal[1])
          
        if abs(s_min) > abs(s_max):
            signal[1] = -signal[1]             #flip signal if negative
            s_max = abs(s_min)                  
            
        if s_max < noise*3:
            return 'too noisy signal'          #ERROR: signal too noisy
        
        threshold = max(0.3*s_max, noise*3)    # 

        for i in range(signal.shape[-1]):
            if signal[1,i] > threshold:
                rise_indices.append(i)
                break
           
    indent = max(rise_indices)-min(rise_indices) 

    for n in range(len(rise_indices)-1):
        shift = rise_indices[n]-rise_indices[-1]
        np.roll(signals[n], -shift)
        
    signals = signals[:,:,indent:-indent-1]         #cutting non-overlapped edges of signals
    
    aver_signal = signals[0].copy()                 #time_axis
    aver_signal[1] = np.mean(signals[:,1], axis=0)  #averaged_signal_axis
    
    return aver_signal
        
#=============================================================================#

if __name__ == "__main__":      #for autonomous work        
    dirs = os.listdir(".")
    for folderpath in dirs:
        if os.path.isdir(folderpath):           #folders
            signals = read_LIIdir(folderpath)   #read files from folder
            if  isinstance(signals, str):
                print('Error processing folder ' + folderpath + ': '+signals)
                continue
            result_signal = average_LIIsignals(signals)     #averaging signal
            if  isinstance(result_signal, str):
                print('Error processing folder ' + folderpath + ': '+result_signal)
                continue
            resultpath = folderpath + '.csv'
            try:
                with open(resultpath, 'w') as resultfile:
                    for i in range(result_signal.shape[-1]):
                        resultfile.write(str(result_signal[0,i]) +', '+
                                         str(result_signal[1,i])+'\n')       
                print('{} signals from folder {} were averaged. File {} was written'.format(signals.shape[0],
                                                   folderpath, resultpath))
                resultfile.close()
            except PermissionError:
                print('ERROR: Access denied for file' + resultpath)
    input()
#=============================================================================#
    
    

   
