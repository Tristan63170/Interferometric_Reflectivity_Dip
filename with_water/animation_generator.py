# this script generate an animation of the beam narrowing effect, representing the magnetic field module as function of the incidence angle in water
# the animation purpose is to illustrate the effect, the optimal gold layer thickness is set, the optimal waist is considered, a range of incidences around the optimal one is scanned
# only the reflection is simulated, the animation is made from the field right after it

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the working directory to the script's location
os.chdir(script_dir)

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftshift, ifft2
import PyMoosh as PM
import my_module as my_mod
import imageio
import pickle

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)

# function to display a progress bar in the for loops
def progress_bar(iteration, total, bar_length=40):
    progress = iteration / total
    arrow = '=' * int(round(progress * bar_length))
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f'\r[{arrow}{spaces}] {int(progress * 100)}%')
    sys.stdout.flush()

filenames = [] # initialisation of a list used to store the file names
max_H = [] # initialisation of a list used to build the common colormap

# initialisation of the needed parameters, other than the incidence,
n1 = 1 # air refractive index
n2 = 1.515  # glass refractive index
wavelength = 632.8
k0 = 2 * np.pi / wavelength
N = 100
alpha0 = 0
beta0 = 0
x0 = 0
y0 = 0
w = 46089
wx = w
wy = w
d = 400000

# the x and y-axes considered
x = np.linspace(-N, N, 2 * N + 1) * (d / 2) / N
y = np.linspace(-N, N, 2 * N + 1) * (d / 2) / N
# grids of the axes
X, Y = np.meshgrid(x, y)

# definition of the spectral grid
pask = 2 * np.pi / d  # definition of the spectral step
a = np.arange(-N, N + 1) * pask  # alpha (kx) range in nm-1
b = a  # beta (ky) range in nm-1
A, B = np.meshgrid(a, b)  # alpha and beta grid

# definition of the structure on which the reflection is done
# get the material from refractiveindex database (RII) when possible
Cr = PM.Material(["main", "Cr", "Rakic-BB"], specialType="RII")
Au = PM.Material(["main", "Au", "Rakic-BB"], specialType="RII")
mat = ["BK7",Cr,Au,1.33**2]  # list with all the materials needed, BK7 taken from PM .json because the one in RII don't match the wavelength, and 1.33 is for water
stack = [0, 1, 2, 3]  # stacking order of the materials in the  multilayered structure
thickness = [4000, 2, 47, 1000]  # thickness of each layer
structure = PM.Structure(mat, stack, thickness, verbose=False)

inc_range = np.linspace(72.5,74.5, int((74.5-72.5)/0.01)+1)
kz = ((n2 * k0) ** 2 - A**2 - B**2) ** (1 / 2)  # determination of the wave vector z component (kz or gamma) for all incident plane wave
print("Simulation of the reflection for all the incidence angles:")
it_loop = 0
total_it = len(inc_range)
for inc in inc_range:
    # reflection
    theta = np.arcsin(np.sqrt((A * np.cos(np.deg2rad(inc)) + kz * np.sin(np.deg2rad(inc))) ** 2 + B**2)/ (n2 * k0))  # array with the incidence angle for all incident plane wave
    cr = np.zeros((2 * N + 1, 2 * N + 1), dtype=complex)
    # determination of the reflection coefficients for each plane wave
    for i in range(len(a)):
        for j in range(len(b)):
            r, t, R, T = PM.coefficient(structure, wavelength, theta[i,j], 1)
            cr[i, j] = r
    U1 = my_mod.gauss2D(X, Y, wx, wy, x0, y0, alpha0, beta0) # spatial beam before reflection
    U1k = fftshift(fftshift(fft2(U1), axes=0), axes=1) # spectral beam before reflection
    U1ki=cr*U1k # spectral beam after reflection
    U1i=ifft2(U1ki) # spatial beam after reflection
    # Save data in a temporary file to use a common colormap at the end
    with open(f"Plots/for_the_animation/data_{inc}.pkl", "wb") as f:
        pickle.dump(U1i, f)
    
    # in order to be able to use a common colormap and still being able to see the effect we display a root of the module instead of the module
    max_H.append(np.max((np.abs(U1i))))
    # display progress bar
    it_loop += 1
    progress_bar(it_loop, total_it)
print("")
# determination of the upper limit of the global colorbar
max_H_value = np.max(max_H)
print(f"The maximum H field module is {max_H_value}")

print("Save of all the png files:")
it_loop = 0
# load of the pickle files, set of the common colorbar and save of the png files for the animation
for inc in inc_range:
    with open(f"Plots/for_the_animation/data_{inc}.pkl", "rb") as f:
        U1i = pickle.load(f)
    
    # save of the figure in pickle files (temporary) to be able to set a common colorbar afterward
    fig = plt.figure(figsize=(18, 12))
    ax = fig.add_subplot(111)
    absU1i = np.abs(U1i)
    absU1i[absU1i<1e-8]=1e-8
    contour = ax.contourf(X/1000, Y/1000, absU1i, levels=np.logspace(np.log10(1e-8), np.log10(max_H_value), 500), cmap="jet", norm="log")
    contour.set_edgecolor("face")
    cbar = fig.colorbar(contour)
    cbar.set_label("|H|")
    cbar.set_ticks([1e-8, 1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    cbar.set_ticklabels([r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_xlabel(r"x ($\mathrm{\mu m}$)")
    ax.set_ylabel(r"y ($\mathrm{\mu m}$)")
    ax.set_title(f"Magnetic field module after reflection at an incidence of {inc}°")
    plt.savefig(f"Plots/for_the_animation/H_after_reflection_{inc}.png", bbox_inches="tight")
    filenames.append(f"Plots/for_the_animation/H_after_reflection_{inc}.png")
    plt.close(fig)
    # display progress bar
    it_loop += 1
    progress_bar(it_loop, total_it)
print("")
# Gif creation
# with imageio.get_writer("Plots/H_module.gif", mode="I", duration=0.2) as writer:
#     for filename in filenames:
#         image = imageio.imread(filename)
#         writer.append_data(image)

# it is better to use the mp4 format to keep the best color resolution
with imageio.get_writer("Plots/H_module.mp4", fps=10) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

print("Animation created : Plots/H_module.mp4")

# temporary files deletion
for inc in inc_range:
    os.remove(f"Plots/for_the_animation/data_{inc}.pkl")
print("Temporary files deleted")
