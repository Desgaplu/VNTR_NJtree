# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:52:43 2020

@author: pldesgagne
"""


import tkinter as tk
from tkinter import filedialog
import pandas as pd

root= tk.Tk()

canvas1 = tk.Canvas(root, width = 300, height = 300, bg = 'lightsteelblue')
canvas1.pack()

def getExcel ():
    global df
    
    import_file_path = filedialog.askopenfilename(title="Select an Excel File", filetypes=(("Excel files","*.xlsx"),("All files","*.*")))
    df = pd.read_excel (import_file_path)
    print (df)
    root.destroy()
    
browseButton_Excel = tk.Button(text='Import Excel File', command=getExcel, bg='green', fg='white', font=('helvetica', 12, 'bold'))
canvas1.create_window(150, 150, window=browseButton_Excel)

root.mainloop()