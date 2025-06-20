import numpy as np
from matplotlib import pyplot as plt
import scipy.io

def xfer(t, u, y, dt_tol=2e-3):
    n = len(t)
    dt = (t[-1] - t[0])/(n-1)
    dt_err  = (np.diff(t.flatten())/dt)-1
    dt_emax = np.max(np.abs(dt_err))
    if dt_emax > dt_tol:
       plt.plot(t[1:], dt_err)
       plt.show()

       raise Exception('Signals must be sampled at uniform rate!')
    
    uh = u*np.hanning(n)
    yh = y*np.hanning(n)

    power = np.ceil(np.log2(n))
    nfft = int(2**(power))

    xft = np.fft.fft(uh,nfft)
    yft = np.fft.fft(yh,nfft)
    xfun = yft[1:int(nfft/2)+1] / xft[1:int(nfft/2)+1]

    omsamp = 2*np.pi/dt
    omega = (omsamp/nfft) * np.arange(1, int(nfft/2)+1)
    # print(omega)
    # omega = 2*np.pi*np.fft.fftfreq(2*nfft, dt)[1:2048+1]


    return omega, xfun

def fitxf(omega, xfun, dencof, omglim, order):
    print(f"ncr = {np.floor(order / 2)}")
    print(f"nci = {np.ceil(order / 2)}")

    n = len(omega)
    delomg = omega[n - 1] / (n - 1)
    nfit = int(np.floor(omglim / delomg)) + 1
    omgfit = omega[:nfit]

    denom = np.polyval(dencof, 1j * omega)
    numer = xfun * denom

    omgmat = np.zeros((2 * nfit, order + 1), dtype=complex)
    for n in range(order + 1):
        k = order - n
        coln = (1j * omgfit) ** k / denom[:nfit]
        omgmat[:nfit, n] = np.real(coln)
        omgmat[nfit:2 * nfit, n] = np.imag(coln)

    xfvec = np.zeros(2 * nfit, dtype=complex)
    xfvec[:nfit] = np.real(xfun[:nfit])
    xfvec[nfit:2 * nfit] = np.imag(xfun[:nfit])

    numcof = np.linalg.lstsq(omgmat, xfvec, rcond=None)[0]

    numer = np.polyval(numcof, 1j * omega)
    xffit = numer / denom

    return numcof, xffit


def T1():
    data = scipy.io.loadmat('./SPPO_deltae_to_theta.mat')
    q = data["Ptchrt"]
    t = data["Time"]
    e = data["Elevator"]

    q -= q[0]
    e -= e[0]

    plt.plot(t, q)
    plt.plot(t, e)
    plt.show()

    print(xfer(t, e, q))

def T2():
    data = scipy.io.loadmat('./Phugoid.mat')
    q = data["Nz"]
    t = data["Time"]
    e = data["Elevator"]

    tp = t[t>20]
    q = q[t>20]
    e = e[t>20]


    q -= q[0]
    e -= e[0]

    # plt.plot(tp, q)
    # plt.plot(tp, e)
    # plt.show()

    print(xfer(tp, e, q))



if __name__ == "__main__":
    T2()