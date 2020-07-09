import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

# Specifying simulation directory and the directory to save results in
wdir = str(input('Enter simulation directory path: '))
sdir = wdir + str(input('Enter path of storage directory relative to simulation directory: '))

# Specifying snapshot index
sim_index = int(input('Enter snapshot index: '))

# Importing data from the snapshot
part = gizmo.io.Read.read_snapshots(['star', 'gas', 'dark'], 'index', sim_index, assign_hosts_rotation = True, 
                                    simulation_directory = wdir)

# Getting halo properties
halo_properties = ut.particle.get_halo_properties(part, 'all')


# Determining some key properties of the galaxy

# Finding radial distance, temperature, number density, and mass of grid cells
r = part['gas'].prop('host.distance.principal.spherical')[:,0]
T = part['gas'].prop('temperature')
n = part['gas'].prop('number.density')
m = part['gas'].prop('mass')

# Phase diagram of the galaxy

select_halo = (r < halo_properties['radius'])

fig, ax = plt.subplots()
fig.set_size_inches(8,6)

im = ax.hexbin(np.log10(n[select_halo]), np.log10(T[select_halo]), C = m[select_halo], reduce_C_function = np.sum, 
              vmin = 1.0E4, vmax = 2.0E7, norm = colors.LogNorm(), gridsize=40, cmap='viridis')

cb = fig.colorbar(im)
cb.set_label(r'Mass in Bin (M$_{\odot}$)')

ax.set_xlabel(r'log(n)')
ax.set_ylabel(r'log(T)')
plt.title('Phase diagram for the galaxy')

plt.savefig(sdir + 'phase_diag.png')

print('Completed rendering phase diagram')