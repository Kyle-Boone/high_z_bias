import sys
sys.path.append('/home/kboone/software/Galaxy_Bias/Bispectrum/Convert_Bispectrum/')

import gc
import os
import h5py
import numpy as np
import pandas as pd
import scipy.fft as fft
from Visual_Helpers import *
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d as inter
from scipy.interpolate import CubicSpline as spline

L = 2000
N = 2048
N_comp = 256
kF = 2*np.pi/L
kN = N*np.pi/L
kS = kN/2

model_power_file = '../Two_Points_Data/CLASS_power_z1' # Location of power spectrum file
k_model_orig, Pk_model_orig = read_power(model_power_file)
max_k_model = np.max(k_model_orig)
Pk_model_orig = Pk_model_orig * (0.2322)**2 # Normalization to z=8

P_inter = spline(k_model_orig, Pk_model_orig, bc_type='natural')

log_k = np.log(k_model_orig); log_Pk = np.log(Pk_model_orig)
log_P = inter(log_k, log_Pk, bounds_error=False, fill_value='extrapolate')

def P(k, kS=kS, P_inter=P_inter, log_P=log_P, max_k=max_k_model):
    k = np.array(k)
    interpolate = k < max_k_model
    smooth = k > kS
    retvals = np.empty_like(k, dtype=np.float32)

    retvals[interpolate] = P_inter(k[interpolate])

    retvals[~interpolate] = np.exp(log_P(np.log(k[~interpolate])))

    retvals[smooth] = 0

    return retvals

def generate_density_field(N, L, Pk_func, workers=32, only_phase=True):
    """
    Generates a realization of the matter density field δ(x) from a given power spectrum P(k).

    Parameters:
    - N (int): Number of grid points in each dimension.
    - L (float): Physical size of the simulation box (in units consistent with k).
    - Pk_func (callable): Function that returns P(k) given k.

    Returns:
    - delta_x (ndarray): Real-space matter density field δ(x) on a 3D grid.
    """
    with fft.set_workers(workers):
        dx = L / N
        dk = 2*np.pi / L
        k  = np.astype(fft.fftfreq(N, d=dx), np.float32) * 2*np.pi           # angular frequencies
        kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
        k_mag = np.sqrt(kx**2 + ky**2 + kz**2)

        uniq_k_mag, inv = np.unique(k_mag, return_inverse=True)
        inv = inv.astype(np.int32)
    
        # A ──────────────────────────────────────────────────────────────
        vol_prefac = (L)**1.5              # √(Δk³)
        Pk_uniq = Pk_func(uniq_k_mag).astype(np.float32)
        Pk = Pk_uniq[inv]
        if only_phase:
            rnd = np.exp(2*np.pi*1j*np.random.uniform(size=(N,N,N)).astype(np.float32)) * np.sqrt(2)
        else:
            rnd = (np.random.normal(size=(N,N,N)) + 1j*np.random.normal(size=(N,N,N))).astype(np.complex64)
        delta_k = rnd * np.sqrt(Pk/2) * vol_prefac  # correct variance
        delta_k = enforce_hermitian_symmetry_fast(delta_k)
        # ────────────────────────────────────────────────────────────────
    
        # inverse FFT (includes its own 1/N³)
        delta_x = np.real(fft.ifftn(delta_k))
    
        # B ──────────────────────────────────────────────────────────────
        delta_x *= (N / L)**3                          
        # ────────────────────────────────────────────────────────────────
        return delta_x, delta_k, np.array([kx, ky, kz])


def enforce_hermitian_symmetry_fast(a):
    """
    In-place Hermitian symmetry for a full (N,N,N) complex cube `a` in NumPy FFT ordering.
    For each pair (i,j,k) and (−i,−j,−k) mod N, copy only one side (preserves variance).
    Works for even/odd N. Returns `a`.
    """
    N = a.shape[0]
    assert a.shape == (N, N, N)
    h = N // 2

    # partner index along one axis: p[0]=0; p[q]=N-q for q>0
    p = np.arange(N)
    p[1:] = np.arange(N-1, 0, -1)  # [0, N-1, N-2, ..., 1]

    # ---- 1) Pair i with N-i for i=1..h-1 (even N) or 1..h (odd N)
    if N > 2:
        # indices to copy FROM (the "positive" side)
        isrc = np.arange(1, h + (N % 2))      # 1..h-1 (even) or 1..h (odd)
        idst = (N - isrc) % N                  # their partners

        # vectorized assign over j,k using partner mapping on axes 1 and 2
        # a[idst, :, :] = conj( a[isrc, partner(j), partner(k)] )
        a[idst, :, :] = np.conj(a[isrc][:, :, :][:, p][:, :, p])

    # ---- 2) Inside i=0 (and i=h if even), pair j with N-j for j=1..h-1 (even) or 1..h (odd)
    special_is = [0] + ([h] if N % 2 == 0 else [])
    if N > 2:
        jsrc = np.arange(1, h + (N % 2))
        jdst = (N - jsrc) % N
        for i in special_is:
            # a[i, jdst, :] = conj( a[i, jsrc, partner(k)] )
            a[i, jdst, :] = np.conj(a[i, jsrc, :][:, p])

    # ---- 3) On lines where i and j are special, pair k with N-k for k=1..h-1 (even) or 1..h (odd)
    if N > 2:
        ksrc = np.arange(1, h + (N % 2))
        kdst = (N - ksrc) % N
        for i in special_is:
            for j in special_is:
                # a[i, j, kdst] = conj( a[i, j, ksrc] )
                a[i, j, kdst] = np.conj(a[i, j, ksrc])

    # ---- 4) Make only the truly self-conjugate points real
    # (i,j,k) where each coord is 0 (odd N) or in {0, h} (even N)
    for i in special_is:
        for j in special_is:
            for k in special_is:
                a[i, j, k] = a[i, j, k].real + 0j

    return a

d1x, d1k, k = generate_density_field(N,L,P,only_phase=True)
sig2 = np.std(d1x)**2
k_mag = np.linalg.norm(k,axis=0)
k_mag_div = np.copy(k_mag)
k_mag_div[np.where(k_mag_div == 0)] = np.max(k_mag) * 1_000_000

crop = np.full(N, False, dtype=bool)
crop[:int(N_comp/2)] = True
crop[-int(N_comp/2):] = True
crop1, crop2, crop3 = np.meshgrid(crop, crop, crop, indexing='ij')
full_crop = crop1&crop2&crop3

d1k_crop = np.reshape(d1k[full_crop], (N_comp, N_comp, N_comp))
with fft.set_workers(32):
    d1x_crop = np.real(fft.ifftn(d1k_crop)) * ((N_comp/L)**3)
print('For d, average/std: ' + str(np.average(d1x_crop)/np.std(d1x_crop)))
d1x_crop = d1x_crop - np.average(d1x_crop)
np.save('/mnt/nvme1/kboone/Sims/d.npy', d1x_crop)

d2x_NS = d1x**2
with fft.set_workers(32):
    d2k_NS = fft.fftn(d2x_NS) * (L/N)**3
    
    del d2x_NS
    gc.collect()
    
    d2k_NS[k_mag > kS] = 0
    d2x = np.real(fft.ifftn(d2k_NS)) * (N/L)**3

    del d2k_NS
    gc.collect()
    
    d2x = d2x - sig2
    d2k = fft.fftn(d2x) * (L/N)**3

    del d2x
    gc.collect()
    
    d2k_crop = np.reshape(d2k[full_crop], (N_comp, N_comp, N_comp))

    del d2k
    gc.collect()
    
    d2x_crop = np.real(fft.ifftn(d2k_crop)) * ((N_comp/L)**3)
print('For d^2, average/std: ' + str(np.average(d2x_crop)/np.std(d2x_crop)))
d2x_crop = d2x_crop - np.average(d2x_crop)
np.save('/mnt/nvme1/kboone/Sims/d2.npy', d2x_crop)

d3x_NS = d1x**3
with fft.set_workers(32):
    d3k_NS = fft.fftn(d3x_NS) * (L/N)**3
    
    del d3x_NS
    gc.collect()

    d3k_NS[k_mag > kS] = 0
    d3x = np.real(fft.ifftn(d3k_NS)) * (N/L)**3

    del d3k_NS
    gc.collect()
    
    d3x = d3x - 3*sig2*d1x
    d3k = fft.fftn(d3x) * (L/N)**3

    del d3x
    gc.collect()
    
    d3k_crop = np.reshape(d3k[full_crop], (N_comp, N_comp, N_comp))

    del d3k
    gc.collect()
    
    d3x_crop = np.real(fft.ifftn(d3k_crop)) * ((N_comp/L)**3)
print('For d^3, average/std: ' + str(np.average(d3x_crop)/np.std(d3x_crop)))
d3x_crop = d3x_crop - np.average(d3x_crop)
np.save('/mnt/nvme1/kboone/Sims/d3.npy', d3x_crop)

# Generation of d^(2)
with fft.set_workers(32):
    d_2x_NS = (5/7 * d1x**2).astype(np.complex64)
    
    for i in np.arange(len(k)):
        d_2x_NS += (fft.ifftn(k[i] * d1k) * (N / L)**3) * (fft.ifftn(k[i] * d1k / k_mag_div**2) * (N / L)**3)
    
    for i in np.arange(len(k)):
        for j in np.arange(len(k)):
            d_2x_NS += 2/7 * (fft.ifftn(k[i] * k[j] * d1k / k_mag_div**2) * (N / L)**3)**2
    
    d_2k = fft.fftn(d_2x_NS) * (L / N)**3

    del d_2x_NS
    gc.collect()
    
    # d_2k[k_mag > kS] = 0 Unnecessary, already zerod out by the following crop

    d_2k_crop = np.reshape(d_2k[full_crop], (N_comp, N_comp, N_comp))

    del d_2k
    gc.collect()
    
    d_2x_crop = np.real(fft.ifftn(d_2k_crop)) * (N / L)**3
print('For d^(2), average/std: ' + str(np.average(d_2x_crop)/np.std(d_2x_crop)))
d_2x_crop = d_2x_crop - np.average(d_2x_crop)
np.save('/mnt/nvme1/kboone/Sims/d_2.npy', d_2x_crop)

# Generation of K^2 Field
with fft.set_workers(32):
    for i in np.arange(len(k)):
        for j in np.arange(len(k)):
            if i == j:
                if i == 0:
                    K2x_NS = (fft.ifftn(k[i] * k[j] * d1k / k_mag_div**2 - d1k/3) * (N / L)**3)
                K2x_NS += (fft.ifftn(k[i] * k[j] * d1k / k_mag_div**2 - d1k/3) * (N / L)**3)
            K2x_NS += (fft.ifftn(k[i] * k[j] * d1k / k_mag_div**2) * (N / L)**3)

    K2k_NS = fft.fftn(K2x_NS) * (L/N)**3
    
    del K2x_NS
    gc.collect()
    
    K2k_NS[k_mag > kS] = 0
    K2x = np.real(fft.ifftn(K2k_NS)) * (N/L)**3

    del K2k_NS
    gc.collect()
    
    K2x = K2x - 2*sig2/3
    K2k = fft.fftn(K2x) * (L/N)**3

    del K2x
    gc.collect()
    
    K2k_crop = np.reshape(K2k[full_crop], (N_comp, N_comp, N_comp))

    del K2k
    gc.collect()
    
    K2x_crop = np.real(fft.ifftn(K2k_crop)) * ((N_comp/L)**3)
print('For K^2, average/std: ' + str(np.average(K2x_crop)/np.std(K2x_crop)))
K2x_crop = K2x_crop - np.average(K2x_crop)
np.save('/mnt/nvme1/kboone/Sims/K2.npy', K2x_crop)