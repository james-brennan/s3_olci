"""
brdf.py
This performs a BRDF correction seperately from the regularisation
of NBAR

Maybe consider a c-factor correction too?
"""
import numpy as np
import scipy.optimize

global MODIS_C_obs
# This is now true variance...
MODIS_C_obs = np.array([0.004, 0.015, 0.003, 0.004, 0.013, 0.010, 0.006])**2

OLCI_C_obs = np.array([0.004, 0.003, 0.004, 0.015])**2

def do_brdf_corr_multi_sns(doys, qa, refl, K):
    """
    try a simple fast solver...
    """
    iso = np.zeros_like(refl)
    nBands = refl.shape[1]
    band_rmse = np.zeros(nBands)
    for band in range(nBands):
        fix_qa = np.logical_and(qa, refl[:, band] > 0)
        y = refl[fix_qa, band]
        KK = K.T[fix_qa]
        #x = scipy.optimize.lsq_linear(KK, y, bounds=(0, 1), method='bvls').x
        x = np.linalg.lstsq(KK, y, rcond=None)[0]
        if np.logical_or(np.any(x<0), np.any(x>0.6)):
            # soemething wrong with this do a constrained one
            x = scipy.optimize.lsq_linear(KK, y, bounds=(0, .6), method='bvls').x
        res = (KK.dot(x) - y)**2
        rmse = np.sqrt(res.sum() / res.shape[0])
        band_rmse[band] = rmse
        #import pdb; pdb.set_trace()
        brdf = KK[:, 1:].dot(x[1:])
        _iso = y - brdf
        iso[fix_qa, band] = _iso
        #Z[fix_qa, band]=z
        #import pdb; pdb.set_trace()
    return iso, band_rmse


def do_brdf_corr_multi_sns_old(doys, qa, refl, K):
    """
    try a simple fast solver...
    """
    iso = np.zeros_like(refl)
    nBands = refl.shape[1]
    band_rmse = np.zeros(nBands)
    for band in range(nBands):
        fix_qa = np.logical_and(qa, refl[:, band] > 0)
        y = refl[fix_qa, band]
        KK = K.T[fix_qa]
        #x = scipy.optimize.lsq_linear(KK, y, bounds=(0, 1), method='bvls').x
        #x = scipy.optimize.least_squares(KK, y, bounds=(0,1), loss='soft_l1')
        #import pdb; pdb.set_trace()
        x = np.linalg.solve(KK.T.dot(KK), KK.T.dot(y))
        #x = np.linalg.lstsq(KK, y, rcond=None)[0]
        # do BRDF correction
        res = (KK.dot(x) - y)**2
        rmse = np.sqrt(res.sum() / res.shape[0])
        band_rmse[band] = rmse
        # z score -- remove bad obs
        z = np.sqrt(res) / rmse
        z_T = 2
        zTfrac = 100 * (z > z_T).sum() / z.shape[0]
        if zTfrac < 10:
            y[z > z_T] = 0
        else:
            z[:]=0
        # run it again to fix estimate
        KK = K.T[fix_qa][y > 0]
        y = y[y > 0]
        #x = np.linalg.lstsq(KK, y, rcond=None)[0]
        x = scipy.optimize.lsq_linear(KK, y, bounds=(0, 1), method='bvls').x
        #import pdb; pdb.set_trace()
        brdf = KK[:, 1:].dot(x[1:])
        _iso = y - brdf
        tmp = iso[fix_qa, band]
        tmp[z <= z_T] = _iso
        iso[fix_qa, band] = tmp
        #Z[fix_qa, band]=z
        #import pdb; pdb.set_trace()
    return iso, band_rmse

