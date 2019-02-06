import numpy as np
import os

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import font

from groundf import Cp2Cp_int

def LIndI_greetings():
    print('*'*72)
    print(' '*33, 'LIndI')
    print('*'*72)
    print('\nLIndI is starting…')

def confrim_settings(partfilepath, mixfilepath, lasfilepath, detfilepath):
    """confrim reading from settings.inp"""
    if all((partfilepath, mixfilepath, lasfilepath, detfilepath)):
        print('Successfully read settings from settings.inp')
    else:
        print('ERROR reading settings from settings.inp\n\nPress any key to exit')
        input()
        sys.exit()
        
def confirm_data(part_data, mix_data, la_data, det_data):
    """confrim reading from data files"""
    if all((part_data, mix_data, la_data, det_data)):
        print('Successfully read data')
    else:
        print('ERROR')
        input()
        sys.exit()

def ask_for_signals(det_data):
    """GUI for signals input"""
    (_, band_1, band_2, _) = det_data
    signal_path_1 = signal_path_2 = ''
    while not all((signal_path_1, signal_path_2)):
        print(
              '\nPlease input paths to signal files\n'
         + '(signal 1 in {:.0f}-{:.0f} nm,'.format(band_1[0]*1e9, band_1[1]*1e9)
         + ' signal 2 in {:.0f}-{:.0f} nm)'.format(band_2[0]*1e9, band_2[1]*1e9)   
              )
        root = Tk()
        root.withdraw()                                      #hiding root window
        root.filename = filedialog.askopenfilename(
              initialdir = "/",
              title = "Select file 1 ({:.0f}-{:.0f} nm)".format(band_1[0]*1e9,
                                                                band_1[1]*1e9),
              filetypes = (("csv files","*.csv"),("all files","*.*"))
              )
        print ('Signal 1 file: {0}'.format(root.filename))
        signal_path_1 = root.filename                           
        fst_dir = signal_path_1.rpartition('/')[0]             #remember 1st dir
        root.filename = filedialog.askopenfilename(
              initialdir = fst_dir,                            #staying same dir
              title = "Select file 1 ({:.0f}-{:.0f} nm)".format(band_2[0]*1e9,
                                                                band_2[1]*1e9),
              filetypes = (("csv files","*.csv"),("all files","*.*"))
              )
        print ('Signal 2 file: {0}'.format(root.filename))
        signal_path_2 = root.filename     
        root.destroy()
    return signal_path_1, signal_path_2

def ask_for_la_energy(la_data):
    """GUI ask for laser energy in given experiment"""
    (la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data) = la_data
    
    def save_and_quit(event, *args):
        """closing the window"""
        energy_float.set(float(energy_value.get()) * 1e-3) #milliJoules to Joules
        root.destroy()
    
    print('\nSpecify laser energy:')
    
    root = Tk()
    energy_float = DoubleVar(root)
    energy_float.set(la_energy)
    root.title("Input laser energy")
    mainframe = ttk.Frame(root, padding="10 10 10 10")
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    
    energy_value = StringVar()
    
    energy_value.set('{:.1f}'.format(la_energy*1e3)) #Joules to milliJoules
       
    ttk.Label(mainframe, text="Energy = ").grid(column=0, row=0, sticky=E)
    E_entry = ttk.Entry(mainframe, width=7,
                        textvariable=energy_value).grid(column=1,row=0,sticky=E)
    ttk.Label(mainframe, text="mJ").grid(column=2, row=0, sticky=W)
    ok_button = ttk.Button(mainframe, text = 'OK')
    ok_button.grid(row = 0, column = 4)
    ok_button.bind("<Button-1>", save_and_quit)
    
    root.mainloop()
    
    la_energy = energy_float.get()
    print('I = {:.1f} mJ'.format(la_energy*1e3)) #Joules to milliJoules
    la_data = (la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data)
    
    return la_data

def ask_for_T0_P0(mix_data):
    """GUI ask for T0 and P0 in given experiment"""
    
    (composition, gas_weight, gas_Cp_data,
     gas_Cpint_data, alpha_data, T0, P0) = mix_data #unpacking data
    
    def save_and_quit(event, *args):
        """closing the window"""
        T0_float.set(float(T0_value.get())) 
        P0_float.set(float(P0_value.get())*1e5)     #bar to Pa 
        root.destroy()
    print('\nSpecify ambient gas temperature and pressure:')
    root = Tk()
    T0_float = DoubleVar(root)
    T0_float.set(T0)
    P0_float = DoubleVar(root)
    P0_float.set(P0)
    root.title("Input T0 and P0")
    mainframe = ttk.Frame(root, padding="10 10 10 10")
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    
    T0_value = StringVar()
    P0_value = StringVar()
    
    T0_value.set('{:.0f}'.format(T0))
    P0_value.set('{:.3f}'.format(P0*1e-5))   #Pa to bar
       
    ttk.Label(mainframe, text="T0 = ").grid(column=0, row=0, sticky=E)
    E_entry = ttk.Entry(mainframe, width=7, textvariable=T0_value).grid(column=1, row=0, sticky=E)
    ttk.Label(mainframe, text="K").grid(column=2, row=0, sticky=W)
    
    ttk.Label(mainframe, text="P0 = ").grid(column=0, row=1, sticky=E)
    E_entry = ttk.Entry(mainframe, width=7, textvariable=P0_value).grid(column=1, row=1, sticky=E)
    ttk.Label(mainframe, text="bar").grid(column=2, row=1, sticky=W)
    
    
    ok_button = ttk.Button(mainframe, text = 'OK')
    ok_button.grid(row = 1, column = 4)
    ok_button.bind("<Button-1>", save_and_quit)
    
    root.mainloop()
    
    T0 = T0_float.get()
    P0 = P0_float.get()
    
    print('T = {:.0f} K  |  P = {:.3} bar'.format(T0, P0*1e-5))  #Pa to bar
    
    
    gas_Cpint_data = Cp2Cp_int(gas_Cp_data, T0)

    mix_data = (composition, gas_weight, gas_Cp_data,
                gas_Cpint_data, alpha_data, T0, P0)     #packing data
    
    return mix_data

   
