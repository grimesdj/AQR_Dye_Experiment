#patch Dolfyn compatibility for NumPy 2.0 and SciPy
import numpy as np

np.NaN = np.nan  #restore np.NaN
from numpy.polynomial.polyutils import RankWarning as _RankWarning
np.RankWarning = _RankWarning
import scipy.integrate as _si
def _cumtrapz(y, x=None, initial=0.0, axis=-1):
    y = np.asarray(y)
    if x is None:
        dx = 1.0
        trap = (y[:-1] + y[1:]) * dx * 0.5
    else:
        x = np.asarray(x)
        dx = np.diff(x)
        trap = (y[:-1] + y[1:]) * dx * 0.5
    return np.concatenate(([initial], np.cumsum(trap, axis=axis)), axis=axis)
_si.cumtrapz = _cumtrapz

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import welch
from datetime import datetime, timedelta
import contextlib


import io

# +
#load ADCPDolfynDataset
print("Setting up imports...")
#sys.path.append(r"/Users/eberreyes/Documents/adcpy5")
#sys.path.append(r"/Users/derekgrimes/git/kelp/Python_code/adcpy5")
sys.path.append(r"./adcpy5")


from ADCPDolfynDataset import ADCPDolfynDataset
# -

#load M3
print("Changing directory...")
#os.chdir(r"/Users/eberreyes/Downloads/")
#os.chdir(r"/Users/derekgrimes/OneDriveUNCW/KELP-vingilote/data/FullExperiment/raw/RDI")
os.chdir(r"../../../Kelp_data/data/FullExperiment/raw/RDI")
f = 'RDI_SN_1071_8m_data.000'
nens = 999  #load all ensembles
a = ADCPDolfynDataset()

print("Reading M3 binary file...")
with contextlib.redirect_stdout(io.StringIO()):
    a.read_raw(raw_file=f, nens=nens)
print("Successfully read M3!")

# save as matfile
from scipy.io import savemat

mdict = {
    name: a.Dataset[name].values
    for name in a.Dataset.data_vars
}

mdict.update({
    name: a.Dataset.coords[name].values
    for name in a.Dataset.coords
})
savestr = f"ADCP_M3_{nens}_ens.mat"
savemat(savestr, mdict)


return

#extract velocity and range
print("Extracting velocity and range...")
vel = a.velocity
rng = -a.bin_center_elevation  #flip to positive Height Above Bottom (HAB)
u = vel[:, :, 0]  # East
v = vel[:, :, 1]  # North

#remove gray area (>8m above bottom)
print("Removing gray area (above 8 meters from bottom)...")
mask = rng <= 8
rng_filtered = rng[mask]
u_filtered = u[:, mask]  #apply mask
v_filtered = v[:, mask]

print(f"u_filtered shape: {u_filtered.shape}")
print(f"rng_filtered shape: {rng_filtered.shape}")

#rotate velocities along principal axis
print("Computing principal axis rotation...")
theta = 0.5 * np.arctan2(2 * np.nanmean(u_filtered * v_filtered),
                         np.nanmean(u_filtered**2 - v_filtered**2))
u_rot = u_filtered * np.cos(theta) + v_filtered * np.sin(theta)

#build correct timestamps manually

print("Building manual timestamps...")
start_date = datetime(2024, 7, 1, 0, 0, 0)  #start of deployment
dt_seconds = 1  #assuming 1 second per ping (can be adjusted later)
times = [start_date + timedelta(seconds=i*dt_seconds) for i in range(u_rot.shape[0])]
times = pd.to_datetime(times)
print(f"Time range: {times.min()} to {times.max()}")

#build DataFrame
df_vel = pd.DataFrame({'Date': times})
for i in range(u_rot.shape[1]):
    df_vel[f'bin_{i+1}'] = u_rot[:, i]

#resample to 5-minute averages
print("Resampling to 5-minute averages...")
df_vel = df_vel.set_index('Date').resample('5min').mean().reset_index()

print(f"Number of points after resampling: {len(df_vel)}")
print(df_vel.head())

print("Plotting principal velocity (M3)...")
fig, ax = plt.subplots(figsize=(12, 6))
X, Y = np.meshgrid(df_vel['Date'], rng_filtered)

pc = ax.pcolormesh(X, Y, df_vel.iloc[:,1:].T, shading='auto', cmap='coolwarm', vmin=-0.05, vmax=0.05)
plt.colorbar(pc, label='Velocity (m/s)')
ax.set_xlabel('Time')
ax.set_ylabel('Height Above Bottom (m)')
ax.set_title('M3 Principal Flow Velocity (5-min average)')
ax.invert_yaxis()
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()

#compute and plot spectrum
print("Computing and plotting spectrum...")
u_ave = np.nanmean(df_vel.iloc[:,1:].values, axis=1)

fs = 1/300  # 5-min interval = 300 seconds
nperseg = min(len(u_ave)//5, len(u_ave))
f, Pxx = welch(u_ave, fs=fs, nperseg=nperseg)

plt.figure(figsize=(8,6))
plt.loglog(f, Pxx)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density (m²/s²/Hz)')
plt.title('Spectra of Principal Velocity (M3)')
plt.grid(True)
plt.tight_layout()
plt.show()

print("M3 principal velocity plotting complete!")

