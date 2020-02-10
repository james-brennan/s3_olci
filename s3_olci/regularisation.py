import numpy as np
from kernels import *
import scipy.linalg
import glob
import datetime
import numba as nb
import matplotlib.pyplot as plt
import time
np.seterr(all='ignore')
import scipy.sparse as sp
import scipy.sparse.linalg as sl
import bandmat as bm  # nice package!
import scipy.stats
import scipy.ndimage
from scipy.fftpack import fft, fftshift
import scipy.signal


MODIS_C_obs = np.array([0.004, 0.015, 0.003, 0.004, 0.013, 0.010, 0.006])**2
OLCI_C_obs = np.array([0.004, 0.003, 0.004, 0.015])**2



@nb.njit
def prepare_iso(doy, iso, nT):
    # Get a optimal 1 ob per day record...
    # not ideal of course...
    uni_doys = np.unique(doy)
    iso_ret = 0*np.zeros((nT, iso.shape[1]))
    for idoy, _doy in enumerate(uni_doys):
        idx = np.where(doy == _doy)[0]
        # take the mean? better for SNR
        ii = iso[idx]
        #vv = vza[idx]
        #_z = Z[idx]
        for band in range(iso.shape[1]):
            r = np.nanmean(iso[:, band][idx])
            iso_ret[idoy, band] = r
    return iso_ret

def prepare_mats(refl, alpha=1e3):
    """
    Does some prep to produce matrices for the TDMA
    """
    y = refl
    Iobs = y > 0
    nT, ys, xs = y.shape
    # Create mats
    AC, BC, CC, DC = [  np.zeros([nT-1, ys, xs], dtype=np.float32),
                        np.zeros([nT, ys, xs], dtype=np.float32),
                        np.zeros([nT-1, ys, xs],  dtype=np.float32),
                        np.zeros([nT, ys, xs],  dtype=np.float32)]
    # Create Diagonals matrices
    # for ease just use bandmat for this for now
    # Form diagonals
    I = np.eye(nT)  #np.diag (np.ones(nx))
    D = (I - np.roll(I, -1)).T
    #D = D.T.dot(D)
    D1A = np.zeros((nT, 2))
    D1A[1:, 0] = np.diag(D, 1)
    D1A[:, 1] = np.diag(D, 0)
    # convert to banded matrices
    D = bm.BandMat(0, 1, D1A.T, transposed=False)
    D2 = bm.dot_mm(D.T, D)
    a = alpha * D2.data[2][:-1]
    b = alpha * D2.data[1]
    c = alpha * D2.data[0][1:]
    # add these
    AC[:]=a[:, None, None]
    BC[:]=b[:, None, None]+ Iobs
    CC[:]=c[:, None, None]
    DC[:]=y
    return AC, BC, CC, DC


def TDMA_MAT(A_, B_, C_, D_):
    """
    Thomas algorithm for a matrix
    Arguments need to be supplied correctly
    to keep this code generic

    AC -- matrix for the lower diagonals
    BC -- matrix for the middle diagonals
    CC -- matrix of the upper diagonals
    DD -- matrix of the RHS
    """
    # these two get over-written
    DC = np.copy(D_)
    BC = np.copy(B_)
    AC = np.copy(A_)
    CC = np.copy(C_)
    nT, ys, xs = DC.shape
    for it in range(1, nT):
        MC = AC[it-1]/(BC[it-1])
        BC[it] = (BC[it]) - MC*CC[it-1]
        DC[it] = DC[it] - MC*DC[it-1]
    XC = BC
    XC[-1] = DC[-1]/(BC[-1])
    for il in range(nT-2, -1, -1):
        XC[il] = (DC[il]-CC[il]*XC[il+1])/(BC[il])
    return XC


def retrieve_along_axis(ys, xs, refl, IDX, PRES, POSTS):
    for y in range(ys):
        for x in range(xs):
            # this is wrong still
            idx = IDX[y, x]
            # clip it
            idx = np.clip(idx, 30, 60)
            pre = refl[:, y, x][idx-16:idx]
            post = refl[:, y, x][idx+1:idx+17]
            # step check
            PRES[:, y,x]=pre
            POSTS[:, y,x]=post


#@numba.njit
def refine_edges(W, X,refl):
    """
    This updates the edges W to stengthen drops etc
    """
    nT, ys, xs= X.shape
    PRES = np.zeros((16, ys, xs))
    POSTS = np.zeros((16, ys, xs))
    IDX = np.nanargmin(W[30:61], axis=0)
    IDX+=30
    # fill these matrices
    retrieve_along_axis(ys, xs, refl, IDX, PRES, POSTS)
    PRES[PRES==0]=np.nan
    POSTS[POSTS==0]=np.nan
    """
    do the MAD test
    """
    dch = np.nanmedian(POSTS, axis=0) - np.nanmedian(PRES, axis=0)
    preMAD = 1.426 * np.nanmedian((np.abs(np.nanmedian(PRES, axis=0)-PRES)), axis=0)
    postMAD = 1.426 * np.nanmedian((np.abs(np.nanmedian(POSTS, axis=0)-POSTS)), axis=0)
    # make robust Z SCORE
    Sch =  np.abs(dch/(preMAD + postMAD))
    # check we actually have some obs for this
    COND = np.logical_and((PRES>0).sum(axis=0)>3, (POSTS>0).sum(axis=0)>3)
    KeepEdges=COND
    for y in range(ys):
        for x in range(xs):
            idx = IDX[y,x]
            if Sch[y,x]>1 and COND[y,x]:
                """
                So this is probably a step
                fix it up a bit
                """
                W[:, y, x] = 1
                W[:, y, x][idx] = 0.01
            # and if the COND is not met remove the step
            if not COND[y,x]:
                pass
                #W[:, y,x]=1
                #Sch[y,x]=0
    return W, KeepEdges, Sch


def solve_band(refl, solve_edge=True, W=None, alpha=1000, unc=False, drop=True):
    """
    test run function
    """
    Iobs = (refl>0)
    nObs = (Iobs).sum(axis=0)
    # scale alpha
    alpha = 1
    alpha = alpha * np.nanmean(nObs)
    alpha = np.clip(alpha, 20, 250)
    if solve_edge:
        T = .6 * alpha
        nT, ys, xs = refl.shape
        # make D matrix
        I = np.eye(nT)  #np.diag (np.ones(nx))
        D = (I - np.roll(I, -1)).T
        #D = D.T.dot(D)
        D1A = np.zeros((nT, 2))
        D1A[1:, 0] = np.diag(D, 1)
        D1A[:, 1] = np.diag(D, 0)
        # convert to banded matrices
        _D = bm.BandMat(0, 1, D1A.T, transposed=False)
        # create mats
        A, B, C, D = prepare_mats(refl, alpha=alpha)
        # Initial run with no smoothing
        X = TDMA_MAT(A, B, C, D)
        W0 = np.ones_like(X)*100
        WW = np.ones_like(X)
        X0 = np.zeros_like(D)
        CONV = np.zeros((ys, xs)).astype(np.bool)
        Nits =  np.zeros((ys, xs), dtype=np.int)
        for i in range(5):
            # calculate weights
            order_n = 1
            #dx = np.diff(X, axis=0, n=order_n)* alpha**2
            dx = np.gradient(X, axis=0) * alpha**2
            # enforce direction constraint
            if drop != None:
                if drop:
                    dx[dx>0]=0
                else:
                    dx[dx<=0]=0
            # eg frst 30 days and last 30 days
            _w = np.exp(-dx**2 / T)
            """
            Scaling from CERC
            """
            y,x=30, 57
            wScale = nObs / np.sum(_w, axis=0)
            _w = _w * wScale
            _w = np.clip(_w, 0., 1)
            np.place(_w, np.isnan(_w), 1)
            WW[:]=_w
            """
            Update the regularisation matrix
            with the weights

            -- This involves a dot product
            with a diag eg D.T W D
                so i don't how to do this atm

            This is the slowest so if this can be done fast much better
            """
            for y in range(ys):
                for x in range(xs):
                    if not CONV[y,x]:
                        # no convergence keep going
                        n = alpha * bm.dot_mm(_D.T, _D, diag=WW[:, y, x].astype(np.float))
                        a = n.data[2][:-1]
                        b = n.data[1]  + Iobs[:, y, x]
                        c = n.data[0][1:]
                        # update matrices
                        A[:, y, x]=a
                        B[:, y, x]=b
                        C[:, y, x]=c
                        Nits[y,x]+=1
            """
            Run again with new weights
            """
            X = TDMA_MAT(A, B, C, D)
            """
            assess convergence
            -eg are the weights still changing
            """
            cow = np.sum(np.abs((WW - W0) / W0), axis=0)+ 1e-6
            CONV[cow<1e-2]=True
            MINITER=3
            if i< MINITER:
                CONV[:]=False
            CONV[WW.min(axis=0)<0.1]=True
            W0 = WW
            X0 = X
        """
        Want to check steps and re-inforce real steps
        """
        W, keptEdges, Sch = refine_edges(WW, X, refl)
        """
        Resolve with refined edges
        """
        for y in range(ys):
            for x in range(xs):
                    n = alpha * bm.dot_mm(_D.T, _D, diag=W[:, y, x].astype(np.float))
                    a = n.data[2][:-1]
                    b = n.data[1]  + Iobs[:, y, x]
                    c = n.data[0][1:]
                    # update matrices
                    A[:, y, x]=a
                    B[:, y, x]=b
                    C[:, y, x]=c
                    Nits[y,x]+=1
        """
        Run again with the refined weights
        """
        #y,x=43, 62
        if unc:
            X = TDMA_MAT(A, B, C, D)
            DD = np.copy(D).astype(np.float32)
            Inv = np.zeros_like(X)
            # only need the middle month tbh
            for t in range(20, nT-20):
                DD[:]=0
                DD[t]=1
                Inv[t]=TDMA_MAT(A, B, C, DD)[t]
            return X, Inv, W, CONV, Nits, Sch
        else:
            X = TDMA_MAT(A, B, C, D)
            return X, W, CONV, Nits, Sch
    else:
        """
        Edge has been provided so use that
        """
        Iobs = (refl>0)
        nObs = (Iobs).sum(axis=0)
        nT, ys, xs = refl.shape
        # make D matrix
        I = np.eye(nT)  #np.diag (np.ones(nx))
        D = (I - np.roll(I, -1)).T
        #D = D.T.dot(D)
        D1A = np.zeros((nT, 2))
        D1A[1:, 0] = np.diag(D, 1)
        D1A[:, 1] = np.diag(D, 0)
        # convert to banded matrices
        _D = bm.BandMat(0, 1, D1A.T, transposed=False)
        # create mats
        A, B, C, D = prepare_mats(refl, alpha=alpha)
        # add the weights matrix
        for y in range(ys):
            for x in range(xs):
                n = alpha * bm.dot_mm(_D.T, _D, diag=W[:, y, x].astype(np.float))
                a = n.data[2][:-1]
                b = n.data[1]  + Iobs[:, y, x]
                c = n.data[0][1:]
                # update matrices
                A[:, y, x]=a
                B[:, y, x]=b
                C[:, y, x]=c
        """
        Run again with the refined weights
        """
        X = TDMA_MAT(A, B, C, D)
        return X, W


def edge_preserving(sensor, refl, band_rmse=None):
    """

    """
    nT, nBands, ys, xs = refl.shape
    alpha=10
    if sensor == "MODIS" or sensor=='VIIRS':
        solutions = np.zeros((nT, 7, ys, xs))
        uncs = np.zeros((nT, 7, ys, xs))
        """
        Do edge preserving on both bands
        """
        iso = np.copy(refl)
        # do 1 and 4 first and get w
        X1, W1, C, N, Z1 = solve_band(iso[:, 1],
                                   alpha=alpha,solve_edge=True, drop=True)
        X4, Inv, W4, C, N, Z4  = solve_band(iso[:, 4],
                                  alpha=alpha, solve_edge=True, unc=True, drop=True)
        pick = np.argmax([Z1, Z4], axis=0)
        idx = np.where(pick==0)
        W4[:, idx[0], idx[1]]=1
        idx = np.where(pick==1)
        W1[:, idx[0], idx[1]]=1
        W = np.minimum(W1, W4)
        for band in range(7):
            X, WW =  solve_band(iso[:, band],
                                  alpha=alpha,solve_edge=False, W=W, )
            # save them
            solutions[:, band] = X
            uncs[:, band] = Inv * MODIS_C_obs[band]
    elif sensor =="OLCI":
        solutions = np.zeros((nT, 4, ys, xs))
        uncs = np.zeros((nT, 4, ys, xs))
        """
        Do edge preserving on both bands
        """
        iso = np.copy(refl)
        X4, Inv, W4, C, N, Z4  = solve_band(iso[:, 3],
                                  alpha=alpha, solve_edge=True, unc=True, drop=True)
        for band in range(4):
            X, WW =  solve_band(iso[:, band],
                                  alpha=alpha,solve_edge=False, W=W4)
            # save them
            solutions[:, band] = X
            uncs[:, band] = Inv * OLCI_C_obs[band]
    return solutions, uncs


