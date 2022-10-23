import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin
from tkinter import *
from Graph import *

root = Tk()
root.title("Stress Element Visualizer")

Fx = Entry(root)
Fy = Entry(root)
Fz = Entry(root)
Vxy = Entry(root)
Vxz = Entry(root)
Vyz = Entry(root)

def draw_graph():
    try:
        stressGraph(int(Fx.get()), int(Fy.get()), int(Fz.get()), int(Vxy.get()), int(Vxz.get()), int(Vyz.get()))
    except:
        print("nah")

l1 = Label(root, text="Normal Stresses")
l2 = Label(root, text="Shear Stresses")

lx = Label(root, text="Fx: ")
ly = Label(root, text="Fy: ")
lz = Label(root, text="Fz: ")

lxy = Label(root, text="Txy: ")
lxz = Label(root, text="Txz: ")
lyz = Label(root, text="Tyz: ")

bt1 = Button(root, text="Solve", command = draw_graph)


l1.grid(row=0, column=1), l2.grid(row=0, column=3)

lx.grid(row=1, column=0), ly.grid(row=2, column=0), lz.grid(row=3, column=0)
lxy.grid(row=1, column=2), lxz.grid(row=2, column=2), lyz.grid(row=3, column=2)

Fx.grid(row=1,column=1), Fy.grid(row=2,column=1), Fz.grid(row=3,column=1)
Vxy.grid(row=1,column=3), Vxz.grid(row=2,column=3), Vyz.grid(row=3,column=3)

bt1.grid(columnspan=4)


root.mainloop()