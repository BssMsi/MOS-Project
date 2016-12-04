from Tkinter import Tk, Frame, Button, Label, Scrollbar, N, S, E, W, HORIZONTAL, StringVar, Canvas, Entry, SUNKEN, \
    IntVar, Radiobutton, RAISED
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plot
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np


# import serial


class AutoScrollbar(Scrollbar):
    # A scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)


class PlotGraph:
    breadth = 1.0
    MAX_STRESS = 44

    def __init__(self):
        def onresize(event):
            data_canvas.configure(scrollregion=data_canvas.bbox("all"))
            self.root.update_idletasks()

        def on_mousewheel(event):
            data_canvas.yview_scroll(-1 * (event.delta / 120), "units")
            self.root.update_idletasks()

        self.root = Tk()
        self.root.wm_title("Micro-Tensile Testing Machine")
        self.root.state("zoomed")

########################################################################################################################
        # Title
        header = Frame(self.root)
        Label(header, text="Micro-Tensile Testing Machine", font=("Helvetica", 40, "bold"), bg="green", relief=RAISED).grid(row=0)
        header.columnconfigure(0, weight=1)
        header.grid(row=0, sticky=E + W)
########################################################################################################################
        # Parameters
        par_frame = Frame(self.root)
        choice_frame = Frame(par_frame, bd=2, relief=RAISED)
        self.choice = IntVar()
        self.choice.set(1)
        Radiobutton(choice_frame, text="ELONGATION", font=("Courier", 15),
                    variable=self.choice, value=1).grid(row=0, column=0, padx=50)
        Radiobutton(choice_frame, text="BENDING", font=("Courier", 15),
                    variable=self.choice, value=2).grid(row=1, column=0, padx=50)
        choice_frame.grid(row=0, column=0, sticky=W)

        spec_frame = Frame(par_frame, bd=2, relief=RAISED)
        Label(spec_frame, text="Specimen\nDimensions", font=("arial", 15)).grid(row=0, column=0, rowspan=2)
        Label(spec_frame, text="LENGTH(mm)", font=("Courier", 15),
              relief=SUNKEN, bg="yellow", fg="brown").grid(row=0, column=1, ipadx=5, ipady=1)
        self.L = Entry(spec_frame, fg="blue", font=("Courier", 15), width=7)
        self.L.insert(0, 25.4)
        self.L.grid(row=0, column=2)
        Label(spec_frame, text="Area(mm^2)", font=("Courier", 15),
              relief=SUNKEN, bg="yellow", fg="brown").grid(row=1, column=1, ipadx=5, ipady=1)
        self.area_box = Entry(spec_frame, fg="blue", font=("Courier", 15), width=7)
        self.area_box.insert(0, 40.3765)
        self.area_box.grid(row=1, column=2)
        spec_frame.grid(row=0, column=1, rowspan=2, sticky=E)
        par_frame.grid_columnconfigure(1, weight=1)
        par_frame.grid(row=1, sticky=E + W)
########################################################################################################################
        # Main Frame containing Curve and data display
        curve = Frame(self.root)

        self.no_reading = 0
        # Stress, Strain value arrays
        self.stress = np.array([0.0])
        self.strain = np.array([0.0])

        # Figure

        main_curve = plot.figure(figsize=(18.5, 10))
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 2])

        ###################################################################
        # # Live diagram
        live_plot = main_curve.add_subplot(gs[0, 0])
        live_plot.set_title("Live Diagram")
        live_plot.grid(True)
        live_plot.set_xlabel("in mm")
        live_plot.set_xlim([0, float(self.L.get())])
        live_plot.set_ylim([0, 50.0])  # Set y axis limits for live diagram

        self.lines = live_plot.plot([0.0], [0.0], 'b-', [0.0], [0.0], 'r--')
        self.live_plot = live_plot
        ###################################################################
        # # Stress-Strain Curve
        curve_plot = main_curve.add_subplot(gs[0, 1])
        curve_plot.set_title("Curve")
        curve_plot.grid(True)
        curve_plot.set_xlim([0.0, 0.01])
        curve_plot.set_ylim([0.0, 100.0])
        curve_plot.set_xlabel("")
        curve_plot.set_ylabel("")

        self.points = curve_plot.plot(self.strain, self.stress, 'k-')
        self.curve_plot = curve_plot
        ###################################################################
        # # Save and show Figure
        main_curve.tight_layout()
        self.main_curve = main_curve
        self.canvas_main = FigureCanvasTkAgg(self.main_curve, master=curve)
        self.canvas_main.show()
        self.canvas_main.get_tk_widget().grid(row=0, column=0, sticky=N + S + E + W)

        # Displaying data
        disp_data = Frame(curve, relief=RAISED, bd=2)
        yscrollbar = AutoScrollbar(disp_data)
        yscrollbar.grid(row=0, column=1, rowspan=4, sticky=N + S)
        xscrollbar = AutoScrollbar(disp_data, orient=HORIZONTAL)
        xscrollbar.grid(row=2, column=0, sticky=E + W)
        data_canvas = Canvas(disp_data, yscrollcommand=yscrollbar.set, xscrollcommand=xscrollbar.set)
        data_canvas.grid(row=0, column=0, rowspan=4, sticky=N + S + W + E)
        yscrollbar.config(command=data_canvas.yview)
        xscrollbar.config(command=data_canvas.xview)
        label_frame = Frame(disp_data)
        Label(label_frame, text="DATA\nS.No  Load(N)  Elongation(mm)\n",
              font=("COMIC SANS MS", 17, "bold")).grid(row=0, column=0, sticky=N + S + W + E)
        self.data = StringVar()
        self.data_str = ""
        self.data.set(self.data_str)
        Label(label_frame, textvariable=self.data,
              font=("Helvetica", 16, "italic")).grid(row=1, column=0, sticky=N + S + W + E)
        data_canvas.create_window(0, 0, window=label_frame)
        label_frame.bind("<Configure>", onresize)
        label_frame.bind_all("<MouseWheel>", on_mousewheel)

        disp_data.rowconfigure(1, weight=1)
        disp_data.columnconfigure(0, weight=1)
        disp_data.columnconfigure(1, weight=1)

        disp_data.grid(row=0, column=1, sticky=N + S + W + E)

        curve.columnconfigure(0, weight=1)
        curve.grid(row=2, sticky=N + S + E + W)
########################################################################################################################
        footer = Frame(self.root)
        Button(footer, text="Run", font=('Comic Sans MS', 25, "bold italic"),
               command=self.run_program).grid(row=0, column=0)
        footer.columnconfigure(0, weight=1)
        footer.grid(row=3, sticky=N + S + E + W)

        self.root.mainloop()

    # Start Program
    def run_program(self):
        # Initialize
        # # Reset all values to initial value
        try:
            self.arrow1.remove()
            self.arrow2.remove()
            self.label1.remove()
            self.label2.remove()

            self.material_type.remove()
            self.zoom.remove()
            self.load.remove()
            self.load_arrow.remove()
            self.delta.remove()
            self.del_arrow.remove()

        except AttributeError:
            # If this is first run, above values have not been initialized, hence AttributeError will be raised
            pass

        self.data_str = ""
        self.no_reading = 0
        L = float(self.L.get())

        self.stress = np.array([0.0])
        self.strain = np.array([0.0])
        self.points[0].set_data(self.strain, self.stress)

        x_shift = 3.0

        self.live_plot.set_xlim([0, L + x_shift + 5.0])
        if self.choice.get() == 1:
            y_shift = 10.0
            self.curve_plot.set_title("Stress-Strain Curve")
            self.curve_plot.set_xlabel("STRAIN ->")
            self.curve_plot.set_ylabel("STRESS (MPa) ->")
            self.live_plot.set_ylim([0, 50.0])
            self.lines[0].set_data([x_shift, L + x_shift, L + x_shift, x_shift, x_shift],
                                   [y_shift, y_shift, y_shift + self.breadth, y_shift + self.breadth, y_shift])
            xlim1, xlim2, ylim1, ylim2 = L + x_shift - 0.5, L + x_shift + 2, y_shift - 0.2, y_shift + self.breadth + 0.2

            self.zoom = zoomed_inset_axes(self.live_plot, 6)
            self.zoom_lines = self.zoom.plot([x_shift, L + x_shift, L + x_shift, x_shift, x_shift],
                                             [y_shift, y_shift, y_shift + self.breadth, y_shift + self.breadth,
                                              y_shift], 'b-',
                                             [0], [0], 'r--')
            self.zoom.grid(True)
            self.zoom.set_xlim([xlim1, xlim2])
            self.zoom.set_ylim([ylim1, ylim2])

            self.load_arrow = self.live_plot.annotate('', xy=(x_shift + L, y_shift + self.breadth / 2),
                                                      xytext=(x_shift + L, y_shift + self.breadth / 2),
                                                      arrowprops=dict(arrowstyle="->"))
            self.load = self.live_plot.text(x_shift + L + 0.25, y_shift + self.breadth / 2, "")

            mark_inset(self.live_plot, self.zoom, loc1=3, loc2=4, fc="none", ec="0.5")
        else:
            y_shift = 1.0
            self.curve_plot.set_title("Force-Displacement Curve")
            self.curve_plot.set_xlabel("DISPLACEMENT (mm) ->")
            self.curve_plot.set_ylabel("FORCE (N) ->")
            self.lines[0].set_data([x_shift, L + x_shift, L + x_shift, x_shift, x_shift],
                                   [y_shift, y_shift, y_shift + self.breadth / 10, y_shift + self.breadth / 10,
                                    y_shift])
            self.live_plot.set_ylim([0, 5.0])
            self.load_arrow = self.live_plot.annotate('', xy=(x_shift + L / 2, y_shift + self.breadth / 10),
                                                      xytext=(x_shift + L / 2, y_shift + self.breadth / 10),
                                                      arrowprops=dict(arrowstyle="<-"))
            self.load = self.live_plot.text(x_shift + L / 2 + 0.2, y_shift * 3 / 2, "")
            self.del_arrow = self.live_plot.annotate('', xy=(x_shift + L / 2, y_shift),
                                                     xytext=(x_shift + L / 2, y_shift),
                                                     arrowprops=dict(arrowstyle="-", linestyle=':'))
            self.delta = self.live_plot.text(x_shift + L / 2 + 0.05, y_shift, "")

        self.main_curve.tight_layout()
        # Dynamic Plot
        with open("raw.csv") as f:
            for line in f:
                load, delta_l = map(float, line.split(","))
                self.update_plot(load, delta_l / 10, L, x_shift, y_shift)

        f.close()
        ind = np.argmax(self.stress)

        if self.choice.get() == 1:
            self.arrow1 = self.curve_plot.annotate('', xy=(self.strain[ind], self.stress[ind]),
                                                   xytext=(self.strain[ind], 0),
                                                   arrowprops=dict(arrowstyle="<->"))
            self.label1 = self.curve_plot.text(self.strain[ind] + 0.001, self.stress[ind] / 2, 'Ultimate Strength')
            self.arrow2 = self.curve_plot.annotate('', xy=(self.strain[self.no_reading], self.stress[self.no_reading]),
                                                   xytext=(self.strain[self.no_reading], 0),
                                                   arrowprops=dict(arrowstyle="<->"))
            self.label2 = self.curve_plot.text(self.strain[self.no_reading] + 0.001,
                                               self.stress[self.no_reading] / 2, 'Fracture')

            self.material_type = self.curve_plot.text(self.curve_plot.get_xlim()[1] / 2, self.stress[ind] + 100
                                                      , 'Ductile Material', fontsize=15,
                                                      color='black', bbox=dict(facecolor='red', edgecolor='black'
                                                                               , boxstyle='round'))

        self.main_curve.canvas.draw()
        self.root.update_idletasks()

    def update_plot(self, load, delta_l, L, x_shift, y_shift):
        # Plot graph
        
        self.no_reading += 1
        self.data_str += '{:3}.\t{:.1f} \t {:=1.6f}\n'.format(self.no_reading, load, delta_l)
        self.data.set(self.data_str)

        # Update live diagram
        # # Elongation
        if self.choice.get() == 1:
            # # Update Live Diagram
            self.lines[1].set_data([L + x_shift, L + x_shift + delta_l, L + x_shift +
                                    delta_l, L + x_shift],
                                   [y_shift, y_shift, y_shift + self.breadth, y_shift + self.breadth])
            self.zoom_lines[1].set_data([L + x_shift, L + x_shift + delta_l, L + x_shift + delta_l, L + x_shift],
                                        [y_shift, y_shift, y_shift + self.breadth, y_shift + self.breadth])

            self.load.set_text("P = {} N".format(load))
            self.load.set_x(L + x_shift + delta_l + 0.25)
            self.load_arrow.xy = (L + x_shift + delta_l + 0.25, y_shift + self.breadth / 2)
            self.load_arrow.xytext = (L + x_shift + delta_l, y_shift + self.breadth / 2)

            # # Update Stress-Strain Curve
            stress = load / float(self.area_box.get())
            strain = delta_l / L

            self.stress = np.append(self.stress, [stress])
            self.strain = np.append(self.strain, [strain])

            self.points[0].set_data(self.strain, self.stress)
            if stress >= np.max(self.stress):
                self.curve_plot.set_ylim([0, stress + 0.2 * stress])
            self.curve_plot.set_xlim([0, strain + 0.2 * strain])

        # # Bending
        # ################ IMP: Note equation used for bending is an approximation for visualization purposes ##########
        elif self.choice.get() == 2:
            # # Update Live Diagram
            # Assume equation f(x) = px2 + qx + r
            # f(x_shift) = y_shift
            # f(x_shift + L) = y_shift
            # (4pr - q2) / 4p = delta_l - y_shift

            p = 4 * delta_l / L ** 2
            q = -p * (L + 2 * x_shift)
            r = y_shift - delta_l + q ** 2 / (4 * p)

            x_points = np.linspace(x_shift, x_shift + L, 10)
            y_points = (x_points ** 2) * p + x_points * q + r

            self.load_arrow.xy = (x_shift + L / 2, y_shift + self.breadth + load / 100000)
            self.del_arrow.xy = (x_shift + L / 2, y_shift - delta_l)
            self.delta.set_y(y_shift - delta_l / 2)
            self.delta.set_text(r'$\delta = {}$'.format(delta_l))
            self.load.set_y(y_shift + (self.breadth + delta_l) / 2)
            self.load.set_text("P = {} N".format(load))
            self.lines[1].set_data(np.append(x_points, np.flipud(x_points) + self.breadth / 10),
                                   np.append(y_points, np.flipud(y_points) + self.breadth / 10))

            # # Update Force-Displacement Curve
            self.stress = np.append(self.stress, [load])
            self.strain = np.append(self.strain, [delta_l])

            self.points[0].set_data(self.strain, self.stress)
            if load >= np.max(self.stress):
                self.curve_plot.set_ylim([0, load + 0.2 * load])
            self.curve_plot.set_xlim([0, delta_l + 0.2 * delta_l])

        self.main_curve.canvas.draw()
        self.root.update()

'''
arduino = serial.Serial('COM1', 115200, timeout=.1)
data = arduino.readline()[:-2] #the last bit gets rid of the new-line chars
arduino.write("Hello from Python!")
while True:
	data = arduino.readline()
	if data:
		print data.rstrip('\n') #strip out the new lines for now
'''
c = PlotGraph()

'''
For all the normal stress-strain curves you see (I guess you could call these "engineering stress-strain curves"), the stress for all points on the curve is computed using the original cross section area.
The curve eventually reaches a maximum peak stress (ultimate stress) and for larger values of strain to the right of the peak, stress starts declining.
This apparent decrease in strength is due to "necking" where the cross section of the test specimen gets smaller due to excessive stretching.
 A true stress-strain curve uses a recomputed cross section area at each point of strain, so the stress plotted is the applied force divided by the instantaneous area instead the original area.
 A true stress-strain curve will show that the actual metal keeps getting stronger all the way to failure (the stress strain curve always has a positive slope),
 but this information has no real practical value outside of the classroom and it is a much harder curve to make. You have to keep measuring the diameter and computing the area for each strain data point.
  In the real world, engineers want to predict how strong things are and when they are going to break. Since you can't change the cross section of things while they are being loaded,
  the engineering stress-strain curve is the one that predicts how part of a real structure will act.
'''