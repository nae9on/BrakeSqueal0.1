import Tkinter
import threading
import matplotlib.backends.backend_tkagg

root = Tkinter.Tk()

class Plotter():
    def __init__(self,fig):
        t = threading.Thread(target=self.PlottingThread,args=(fig,))
        t.start()       

    def PlottingThread(self,fig):        
        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

        Tkinter.mainloop()

if __name__ == "__main__":
    import time

    fig1 = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
    fig1.gca().plot([1,2,3])     
    fig2 = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
    fig2.gca().plot([3,2,1])

    #Shows fig1 and not fig2, just like it's supposed to
    Plotter(fig1)

    time.sleep(1)