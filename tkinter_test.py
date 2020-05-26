# -*- coding: utf-8 -*-
"""
Created on Tue May  5 17:18:52 2020

@author: pldesgagne
"""

from tkinter import *

root = Tk()
root.title("This is my tkinter")
root.iconbitmap('phylotree.ico')

hello_count = 0

def clickMe(hello=None): # when no parameter is passed, take content of entree
    if hello is None:
        hello = 'Hello ' + entree.get() + '!' # take the content of entree
    entree.delete(0, END) # remove content until pos 0 to END
    entree.insert(0, "Who?") # insert at pos 0 the following string
    myLabel = Label(root, text=hello)
    myLabel.pack() # pack the content of hello
    global hello_count # link to the global scope hello_count
    hello_count += 1
    global myLabel1 # use the global myLabel1 instead of a new scope instance
    myLabel1.destroy()
    myLabel1 = Label(frame, text='Hello World! x'+str(hello_count))
    myLabel1.pack()
    
# Creating a Label Widget
frame = LabelFrame(root, text="My Frame..", padx=50, pady=50) # inner padding of the frame
myLabel1 = Label(frame, text="Hello World!") # inside of frame instead of root
myLabel2 = Label(root, text="Enter your name:")
myButton1 = Button(root, text="Click Me!", command=clickMe) # pass the function itself, not a call
myButton2 = Button(root, text="Can't click me!", state=DISABLED)
largeButton = Button(root, text="So big!", padx=100, pady=30, bg='light blue', fg='red', command=lambda:clickMe('So big!'))
entree = Entry(root, borderwidth=5)
button_quit = Button(root, text='Exit', command=root.destroy)


# Shoving into screen
# =============================================================================
# myLabel1.grid(row=0, column=0, columnspan=2) # place lock on a grid
# myLabel2.grid(row=1, column=1)
# myButton1.grid(row=3, column=1)
# myButton2.grid(row=4, column=1)
# largeButton.grid(row=5, column=1)
# =============================================================================

# Packing on screen one after the other
frame.pack(padx=100) # outer padding of the frame
myLabel1.pack() # grid/pack can be used inside of frame
myLabel2.pack()
entree.pack()
myButton1.pack()
myButton2.pack()
largeButton.pack()
button_quit.pack()




root.mainloop()
