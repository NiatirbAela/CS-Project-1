## CS Project 1.py
# App using GUI and PIL
    # to compute strain, stress, and axial force in a single 2D truss element

import math
import tkinter as tk
import numpy as np
from tkinter import ttk, messagebox
from PIL import Image, ImageDraw, ImageTk

#defining functs. to compute math
def compute_stiffness_matrix(E, A, L, theta_deg):
    if L <= 0:
        raise ValueError("Length L must be positive.")  # sanity check

    th = math.radians(theta_deg)        # convert degrees to radians
    c, s = math.cos(th), math.sin(th)   # direction cosines
    k = (E * A) / L                     # axial stiffness multiplier

    # Global stiffness matrix in standard format
    k_global = k * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s],
    ], dtype=float)

    return k_global

def compute_displacement(theta_deg, u1, v1, u2, v2):
    th = math.radians(theta_deg)
    c, s = math.cos(th), math.sin(th)

    u_prime_1 = c * u1 + s * v1
    u_prime_2 = c * u2 + s * v2
    return u_prime_1, u_prime_2

def compute_axial(E, A, L, u_prime_1, u_prime_2):
    """Return (strain, stress, axial_force)."""
    if L <= 0:
        raise ValueError("Length L must be positive.")
    
    eps = (u_prime_2 - u_prime_1) / L
    sigma = E * eps
    N = sigma * A
    return eps, sigma, N

# defining a funct to create image (display later)
def draw_bar_image(theta_deg: float, w=520, h=420) -> Image.Image:
    # initialize the image to use later
    img = Image.new("RGB", (w, h), "white")              #blank image
    d = ImageDraw.Draw(img)                              #Drawing on image
    ox, oy = w // 2, h // 2                              #origin in center
    #set up grid bg
    step = 40                                            #grid spacing
    for x in range(0, w, step):                          #vertical lines
        d.line([(x, 0), (x, h)], fill="#e9eef7")
    for y in range(0, h, step):                          #horizontal lines
        d.line([(0, y), (w, y)], fill="#e9eef7")
    #lying out x and y axis
    d.line([(0, oy), (w, oy)], fill="black", width=2)       # x axis
    d.line([(ox, 0), (ox, h)], fill="black", width=2)       # y axis
    d.text((w - 20, oy - 18), "x", fill="black")            # x label
    d.text((ox + 8, 5), "y", fill="black")                  # y label
    # calculating bar element endpoint - dynamic (display length only)
    th = math.radians(theta_deg)
    Lp = 160                        # arbitrary pixel length (not porportional)
    x2 = ox + Lp * math.cos(th)     # bar element ending x value
    y2 = oy - Lp * math.sin(th)     # bar element ending y value
    # drawing Bar element - dynamic (origin to endpoint)
    d.line([(ox, oy), (x2, y2)], fill="#1f77b4", width=6) 
    # orgin dot from (-5,+5) pixels in x and y - static
    d.ellipse([(ox - 5, oy - 5), (ox + 5, oy + 5)], fill="#1f77b4")
    # Angle arc 
    r = 50                                      
    bbox = [ox - r, oy - r, ox + r, oy + r]     #boundary - static
    ang = -(theta_deg + 180) % 360 - 180        #angle of arc - dynamic
    color = "red" if ang >= 0 else "blue"       # choose color by sign (red ≥ 0, blue < 0)
    start, end = (ang, 0) if ang < 0 else (0, ang) # change quadrants for neg or pos ang
    d.arc(bbox, start=start, end=end, fill=color, width=2) # drawing arc
    d.text((ox + r + 10, oy - 10), f"θ = {theta_deg:+.1f}°", fill=color) #angle label
    return img

# create child class from ttk.frame for our Truss App
class TrussApp(ttk.Frame):                         
    def __init__(self, master): # initalizing truss App class
        super().__init__(master, padding=10)        # initalize parent class
        master.title("2D Truss Element — Axial Stress Calculator")  #window title
        self.grid(sticky="nsew")                    # init geo as grid (dashboard)
        master.columnconfigure(0, weight=1)         # sizing the cols
        master.rowconfigure(0, weight=1)            # sizing the rows
        self._build_ui()                            # init UI setup funct. (use later)
        self._build_display()                       # init display funct. (use later)
        self.load_a()                               # init default a (US Units)
        self.load_b()                               # init default b (Metric Units)

    # defining UI setup as funct.
    def _build_ui(self):
        f = ttk.LabelFrame(self, text="Inputs")                # label input area
        f.grid(row=0, column=0, sticky="nsew", padx=6, pady=6) # place input area

        def add_row(r, label, var):      #creating a funct. to add rows for inputs
            ttk.Label(f, text=label).grid(row=r, column=0, sticky="e", padx=4, pady=2) # label
            e = ttk.Entry(f, textvariable=var, width=16)                          # input area
            e.grid(row=r, column=1, sticky="w", padx=4, pady=2)                # placing input
            return e
        
        # create variables to be used later
        self.sE = tk.StringVar()
        self.sA = tk.StringVar()
        self.sL = tk.StringVar()
        self.sth = tk.StringVar()
        self.su1 = tk.StringVar()
        self.sv1 = tk.StringVar()
        self.su2 = tk.StringVar()
        self.sv2 = tk.StringVar()
        # add input rows: row num, text label, assign variable
        add_row(0, "E (psi or Pa):", self.sE)
        add_row(1, "A (in² or m²):", self.sA)
        add_row(2, "L (in or m):", self.sL)
        add_row(3, "θ (deg):", self.sth)
        ttk.Separator(f, orient="horizontal").grid(row=4, columnspan=2, sticky="ew", pady=6)
        add_row(5, "u₁ (in or m):", self.su1)
        add_row(6, "v₁ (in or m):", self.sv1)
        add_row(7, "u₂ (in or m):", self.su2)
        add_row(8, "v₂ (in or m):", self.sv2)

        # Command buttons - commands defined later
        btns = ttk.Frame(self)
        btns.grid(row=1, column=0, sticky="w", pady=6)
        ttk.Button(btns, text="Compute", command=self.compute).grid(row=0, column=0, padx=4)
        ttk.Button(btns, text="Preset (a)", command=self.load_a).grid(row=0, column=1, padx=4)
        ttk.Button(btns, text="Preset (b)", command=self.load_b).grid(row=0, column=2, padx=4)
        
        # Results label - placeholder
        self.out = tk.StringVar(value="Results will appear here.") # creates a var for use later
        ttk.Label(self, textvariable=self.out, justify="left").grid(row=2, column=0, sticky="w", padx=6, pady=4)

    # defining display setup as funct.
    def _build_display(self):
        # creating the area to display the image function defined before the class
        vis = ttk.LabelFrame(self, text="Bar Element")
        vis.grid(row=0, column=1, rowspan=3, sticky="nsew", padx=6, pady=6)
        self.image_label = ttk.Label(vis)
        self.image_label.grid(row=0, column=0, padx=4, pady=4)
        self.columnconfigure(1, weight=1)  # let the right column stretch
        self.update_display()  # draw initial image, theta = 0
    
    # defining a funct. to update display
    def update_display(self):
        # redrawing the Pillow image to theta = angle
        try:
            theta = float(self.sth.get())
        except ValueError:
            theta = 0.0
        img = draw_bar_image(theta, w=520, h=420)
        tkimg = ImageTk.PhotoImage(img)          #display img on tkinter
        self.image_label.configure(image=tkimg)  #updating the image
        self.image_label.image = tkimg           #saving reference to display

    # setting defaults
    def load_a(self):           # US units
        #using problem 3.18a as a default
        self.sE.set(str(30e6))  # psi
        self.sA.set("2.0")      # in²
        self.sL.set("60.0")     # in
        self.sth.set("45.0")
        self.su1.set("0.0"); self.sv1.set("0.0")
        self.su2.set("0.02"); self.sv2.set("0.04")
        self.out.set("Loaded preset (a). Click Compute.")  #replaces result placeholder

    def load_b(self):           # Metric Units
        #using problem 3.18b as a default
        self.sE.set(str(210e9))  # Pa
        self.sA.set(str(3e-4))   # m²
        self.sL.set("3.0")       # m
        self.sth.set("30.0")
        # displacements in millimeters — auto-converted
        self.su1.set("0.25")
        self.sv1.set("0.0")
        self.su2.set("1.00")
        self.sv2.set("0.0")
        self.out.set("Loaded preset (b). Click Compute.")   #also replace result placeholder

    # Solving
    def compute(self):          #and print
        try:                                #assigning input values to variables
            E = float(self.sE.get())
            A = float(self.sA.get())
            L = float(self.sL.get())
            th = float(self.sth.get())
            u1 = float(self.su1.get())
            v1 = float(self.sv1.get())
            u2 = float(self.su2.get())
            v2 = float(self.sv2.get())
        except ValueError as e:
            messagebox.showerror("Input error", f"Please enter valid numbers.\n{e}")
            return

        # auto-detect units (psi vs Pa)
        if E < 1e9:
            unit = "US"
        else:
            unit = "SI"

        # convert mm → m for SI preset
        if unit == "SI" and max(abs(u1), abs(v1), abs(u2), abs(v2)) > 1e-4:
            u1 /= 1000; v1 /= 1000; u2 /= 1000; v2 /= 1000

        #solve for stiffness matrix
        try:
            K = compute_stiffness_matrix(E, A, L, th)
        except Exception as e: 
            messagebox.showerror("I. Stiffness Error", str(e))
            return
        #solve for displacements
        try: u_prime_1,u_prime_2 = compute_displacement(th, u1, v1, u2, v2)
        except Exception as e:
            messagebox.showerror("II. Displacement Error")
            return
        #solve for strain/stress/force
        try:
            eps, sigma, N = compute_axial(E, A, L, u_prime_1, u_prime_2)
        except Exception as e:
            messagebox.showerror("III. Stress or Force error", str(e))
            return
                            #errors used to pinpoint problems/debug

        def fmt4x4(M):      #formatting the stiffness matrix
            return ("\n".join("  [" + ", ".join(f"{v: .4e}" for v in row) + "]" for row in M))
        
        # correct formatting for results depending on uniits
        if unit == "US":  
            glob_k = f"Global Stiffness K (lb/in):\n{fmt4x4(K)}"
            disp =f"  u'₁ = {u_prime_1:.6f} in\n  u'₂ = {u_prime_2:.6f} in"
            sig_s = f"{sigma:,.4f} psi"
            N_s = f"{N:,.4f} lbf"
        else:
            glob_k = f"Global Stiffness K (N/m):\n{fmt4x4(K)}"
            disp =f"  u'₁ = {u_prime_1:.6f} m\n  u'₂ = {u_prime_2:.6f} m"
            sig_s = f"{sigma/1e6:,.4f} MPa"
            N_s = f"{N/1e3:,.4f} kN"

        # Output Results
        self.out.set(f"{glob_k}\nDisplacement: \n{disp}\n Strain:\n  ε = {eps:.6e}\nStress:\n  σ = {sig_s}\nForce:\n  N = {N_s}  (+: tension)")
        # update display based on results
        self.update_display()


if __name__ == "__main__":      #main funt. to run app
    root = tk.Tk()                  # Creates window for the rest of the code

    style = ttk.Style()             # not necessary 
    try:
        style.theme_use("clam")
    except Exception:
        pass                        #if chosen theme not avail, use sys default

    TrussApp(root)                  # uses window from line 247 to run the app on
    root.mainloop()                 # starts the app, and keeps it running til close
