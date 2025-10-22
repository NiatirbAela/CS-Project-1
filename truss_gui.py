# truss_gui.py
# GUI to compute strain, stress, and axial force in a single 2D truss element
# Works on Python 3.9–3.13+

import math
import tkinter as tk
from tkinter import ttk, messagebox

def compute_axial(E, A, L, theta_deg, u1, v1, u2, v2):
    """Return (strain, stress, axial_force)."""
    if L <= 0:
        raise ValueError("Length L must be positive.")
    theta = math.radians(theta_deg)
    c, s = math.cos(theta), math.sin(theta)

    # Project global displacements onto element axis
    ubar1 = c * u1 + s * v1
    ubar2 = c * u2 + s * v2

    eps = (ubar2 - ubar1) / L
    sigma = E * eps
    N = sigma * A
    return eps, sigma, N


class TrussApp(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=10)
        master.title("2D Truss Element — Axial Stress Calculator")
        self.grid(sticky="nsew")
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self._build_ui()
        self.load_a()

    def _build_ui(self):
        f = ttk.LabelFrame(self, text="Inputs")
        f.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)

        def add_row(r, label, var):
            ttk.Label(f, text=label).grid(row=r, column=0, sticky="e", padx=4, pady=2)
            e = ttk.Entry(f, textvariable=var, width=16)
            e.grid(row=r, column=1, sticky="w", padx=4, pady=2)
            return e

        self.sE = tk.StringVar()
        self.sA = tk.StringVar()
        self.sL = tk.StringVar()
        self.sth = tk.StringVar()
        self.su1 = tk.StringVar()
        self.sv1 = tk.StringVar()
        self.su2 = tk.StringVar()
        self.sv2 = tk.StringVar()

        add_row(0, "E (psi or Pa):", self.sE)
        add_row(1, "A (in² or m²):", self.sA)
        add_row(2, "L (in or m):", self.sL)
        add_row(3, "θ (deg):", self.sth)
        ttk.Separator(f, orient="horizontal").grid(row=4, columnspan=2, sticky="ew", pady=6)
        add_row(5, "u₁ (in or m):", self.su1)
        add_row(6, "v₁ (in or m):", self.sv1)
        add_row(7, "u₂ (in or m):", self.su2)
        add_row(8, "v₂ (in or m):", self.sv2)

        # Buttons
        btns = ttk.Frame(self)
        btns.grid(row=1, column=0, sticky="w", pady=6)
        ttk.Button(btns, text="Compute", command=self.compute).grid(row=0, column=0, padx=4)
        ttk.Button(btns, text="Preset (a)", command=self.load_a).grid(row=0, column=1, padx=4)
        ttk.Button(btns, text="Preset (b)", command=self.load_b).grid(row=0, column=2, padx=4)

        # Results label
        self.out = tk.StringVar(value="Results will appear here.")
        ttk.Label(self, textvariable=self.out, justify="left").grid(row=2, column=0, sticky="w", padx=6, pady=4)

    # ---------- Presets ----------
    def load_a(self):
        """Problem (a): US units"""
        self.sE.set(str(30e6))  # psi
        self.sA.set("2.0")      # in²
        self.sL.set("60.0")     # in
        self.sth.set("45.0")
        self.su1.set("0.0"); self.sv1.set("0.0")
        self.su2.set("0.02"); self.sv2.set("0.04")
        self.out.set("Loaded preset (a). Click Compute.")

    def load_b(self):
        """Problem (b): SI units"""
        self.sE.set(str(210e9))  # Pa
        self.sA.set(str(3e-4))   # m²
        self.sL.set("3.0")       # m
        self.sth.set("30.0")
        # displacements in millimeters — auto-converted
        self.su1.set("0.25")
        self.sv1.set("0.0")
        self.su2.set("1.00")
        self.sv2.set("0.0")
        self.out.set("Loaded preset (b). Click Compute.")

    # ---------- Compute ----------
    def compute(self):
        try:
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

        try:
            eps, sigma, N = compute_axial(E, A, L, th, u1, v1, u2, v2)
        except Exception as e:
            messagebox.showerror("Computation error", str(e))
            return

        if unit == "US":
            sig_s = f"{sigma:,.4f} psi"
            N_s = f"{N:,.4f} lbf"
        else:
            sig_s = f"{sigma/1e6:,.4f} MPa"
            N_s = f"{N/1e3:,.4f} kN"

        self.out.set(f"ε = {eps:.6e}\nσ = {sig_s}\nN = {N_s}  (tension +)")


if __name__ == "__main__":
    root = tk.Tk()
    style = ttk.Style()
    try:
        style.theme_use("clam")
    except Exception:
        pass
    TrussApp(root)
    root.mainloop()
