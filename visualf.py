import numpy as np
import os

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import font

def ask_for_signals():
    """GUI for signals"""
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

def ask_for_la_energy(la_data):
    """GUI ask for laser energy in given experiment"""
    (la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data) = la_data
    
    def save_and_quit(event, *args):
        energy_float.set(float(energy_value.get()) * 1e-3) #milliJoules to Joules
        print('Energy from function = ', energy_float.get())
        root.destroy()
    
    root = Tk()
    energy_float = DoubleVar(root)
    root.title("Input laser energy")
    mainframe = ttk.Frame(root, padding="10 10 10 10")
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    
    energy_value = StringVar()
    
    energy_value.set('{:.0f}'.format(la_energy*1e3)) #Joules to milliJoules
       
    ttk.Label(mainframe, text="Energy = ").grid(column=0, row=0, sticky=E)
    E_entry = ttk.Entry(mainframe, width=7, textvariable=energy_value).grid(column=1, row=0, sticky=E)
    ttk.Label(mainframe, text="mJ").grid(column=2, row=0, sticky=W)
    ok_button = ttk.Button(mainframe, text = 'OK')
    ok_button.grid(row = 0, column = 4)
    ok_button.bind("<Button-1>", save_and_quit)
    
    root.mainloop()
    
    print(energy_float)
    la_energy = energy_float.get()
    print('Energy from function before return = ', la_energy)

    la_data = (la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data)
    
    return la_data

   
