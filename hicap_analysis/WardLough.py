import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('./tests/data')




def WardLoughDepletion(T1, S1, K, lambd, x, y, t, NSteh1):

    # Inverse Fourier transform
    DeltaQ = StehfestCoeff(1, NSteh1) * if1_dQ(T1, S1, K, lambd, np.log(2) / t, x, y)
    for jj in range(2, NSteh1 + 1):
        DeltaQ += StehfestCoeff(jj, NSteh1) * if1_dQ(T1, S1, K, lambd, jj * np.log(2) / t, x, y)
    DeltaQ = 2 * np.pi * lambd * DeltaQ * np.log(2) / t

    return DeltaQ

def WardLoughDrawdown(T1, S1, K, lambd, x, y, t, NSteh1, NSteh2):
    
    # Initialize output arrays
    s1 = np.zeros_like(t)
    s2 = np.zeros_like(t)

    # Inverse Fourier transform
    for ii in range(len(t)):
        try:
            s1[ii] = StehfestCoeff(1, NSteh1) * if1(T1, S1, K, lambd, x, y, np.log(2) / t[ii])
            for jj in range(2, NSteh1 + 1):
                s1[ii] += StehfestCoeff(jj, NSteh1) * if1(T1, S1, K, lambd, x, y, jj * np.log(2) / t[ii])
            s1[ii] *= np.log(2) / t[ii]
        except OverflowError as e:
            print(f"Overflow error in s1 calculation at index {ii}: {e}")
            s1[ii] = np.nan  # Assign NaN if there's an overflow

        try:
            s2[ii] = StehfestCoeff(1, NSteh2) * if2(T1, S1, K, lambd, x, y, np.log(2) / t[ii])
            for jj in range(2, NSteh2 + 1):
                s2[ii] += StehfestCoeff(jj, NSteh2) * if2(T1, S1, K, lambd, x, y, jj * np.log(2) / t[ii])
            s2[ii] *= np.log(2) / t[ii]
        except OverflowError as e:
            print(f"Overflow error in s2 calculation at index {ii}: {e}")
            s2[ii] = np.nan  # Assign NaN if there's an overflow

    return s1, s2

def if1_dQ(T1, S1, K, lambda_, p, x, y):
    return kernel1(T1,S1,K,lambda_,0,0,p)+kernel2(T1,S1,K,lambda_,0,0,p)

def if1(T1, S1, K, lambd, x, y, p):
    G = lambda phi: 2 * (kernel1(T1, S1, K, lambd, x, np.tan(phi), p) +
                          kernel2(T1, S1, K, lambd, x, np.tan(phi), p)) * \
                          np.cos(np.tan(phi) * y) / np.cos(phi) ** 2

    s1InvFour, _ = quad(G, 0, np.pi / 2, epsrel=1e-1, epsabs=1e-1, limit=10000)
    return s1InvFour

def if2(T1, S1, K, lambd, x, y, p):
    H = lambda phi: 2 * (coeff_s1_1(T1, S1, K, lambd, np.tan(phi), p) * 
                          kernel1(T1, S1, K, lambd, x, np.tan(phi), p) +
                          coeff_s1_2(T1, S1, K, lambd, np.tan(phi), p) *
                          kernel2(T1, S1, K, lambd, x, np.tan(phi), p)) * \
                          np.cos(np.tan(phi) * y) / np.cos(phi) ** 2

    s2InvFour, errbnd = quad(H, 0, np.pi / 2, epsrel=1e-1, epsabs=1e-1, limit=10000)
    return s2InvFour

def coeff_s1_1(T1, S1, K, lambd, theta, p):
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = coeffs(T1, S1, K, lambd, theta, p)
    B1 = (mu1 * T1 - b11) / b12
    return B1

def coeff_s1_2(T1, S1, K, lambd, theta, p):
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = coeffs(T1, S1, K, lambd, theta, p)
    B2 = (mu2 * T1 - b11) / b12
    return B2

def kernel1(T1, S1, K, lambd, x, theta_or_y, p):
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = coeffs(T1, S1, K, lambd, theta_or_y, p)

    if x < 0:
        F1 = A1 * np.exp(x * np.sqrt(mu1))
    elif 0 <= x <= 1:
        F1 = A1 * np.exp(-x * np.sqrt(mu1)) + beta1 / (2 * np.sqrt(mu1) * l1) * \
             (np.exp((x - 1) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    else:
        F1 = A1 * np.exp(-x * np.sqrt(mu1)) + beta1 / (2 * np.sqrt(mu1) * l1) * \
             (np.exp((1 - x) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    return F1

def kernel2(T1, S1, K, lambd, x, theta_or_y, p):
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = coeffs(T1, S1, K, lambd, theta_or_y, p)

    if x < 0:
        F2 = A2 * np.exp(x * np.sqrt(mu2))
    elif 0 <= x <= 1:
        F2 = A2 * np.exp(-x * np.sqrt(mu2)) + beta2 / (2 * np.sqrt(mu2) * l2) * \
             (np.exp((x - 1) * np.sqrt(mu2)) - np.exp(-(x + 1) * np.sqrt(mu2)))
    else:
        F2 = A2 * np.exp(-x * np.sqrt(mu2)) + beta2 / (2 * np.sqrt(mu2) * l2) * \
             (np.exp((1 - x) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    return F2

def coeffs(T1, S1, K, lambd, theta_or_y, p):
    b11 = T1 * theta_or_y ** 2 + S1 * p + K
    b12 = -K
    b22 = theta_or_y ** 2 + p + K

    mu1 = (b11 / T1 + b22) / 2 + np.sqrt((b11 / T1 + b22) ** 2 / 4 + (b12 ** 2 - b11 * b22) / T1)
    mu2 = (b11 / T1 + b22) / 2 - np.sqrt((b11 / T1 + b22) ** 2 / 4 + (b12 ** 2 - b11 * b22) / T1)
    l1 = T1 + ((mu1 * T1 - b11) / b12) ** 2
    l2 = T1 + ((mu2 * T1 - b11) / b12) ** 2

    beta1 = (mu1 * T1 - b11) / (b12 * 2 * np.pi * p)
    beta2 = (mu2 * T1 - b11) / (b12 * 2 * np.pi * p)

    Delta = 4 * np.sqrt(mu1 * mu2) + 2 * lambd * (np.sqrt(mu1) / l2 + np.sqrt(mu2) / l1)

    A1 = ((lambd / l2 + 2 * np.sqrt(mu2)) * beta1 * np.exp(-np.sqrt(mu1)) - 
           lambd * beta2 / l2 * np.exp(-np.sqrt(mu2))) / Delta / l1

    A2 = (-lambd * beta1 / l1 * np.exp(-np.sqrt(mu1)) + 
          (lambd / l1 + 2 * np.sqrt(mu1)) * beta2 * np.exp(-np.sqrt(mu2))) / Delta / l2

    return b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2

def safe_factorial(n):
    """Calculate factorial using logarithmic method to avoid overflow."""
    if n < 0:
        return float('inf')
    elif n < 2:
        return 1
    else:
        return np.exp(gammaln(n + 1))

def StehfestCoeff(jj, N):
    LowerLimit = (jj + 1) // 2
    UpperLimit = min(jj, N // 2)

    V = 0
    for kk in range(LowerLimit, UpperLimit + 1):
        denominator = (safe_factorial(N // 2 - kk) *
                       safe_factorial(kk) *
                       safe_factorial(kk - 1) *
                       safe_factorial(jj - kk) *
                       safe_factorial(2 * kk - jj))
        if denominator != 0:  # Prevent division by zero
            V += (kk**(N // 2) * safe_factorial(2 * kk) / denominator)
    
    V *= (-1) ** (N // 2 + jj)
    return V

if __name__ == "__main__":
    T1 = 1.0  # Example parameter
    S1 = 1000  # Example parameter
    K = 0.1   # Example parameter
    lambd = 0.1  # Example parameter
    x = 0.5  # Example parameter
    y =1  # Example parameter
    t = np.logspace(-2, 8, 100)  # Example time array
    NSteh1 = 2  # Example number for Stehfest coefficients
    NSteh2 = 2  # Example number for Stehfest coefficients

    from ward_lough_data import *



    s1, s2 = WardLoughDrawdown(T1, S1, K, lambd, x, y, t, NSteh1, NSteh2)
    
    df = pd.DataFrame(index=t, data={'s1':s1,'s2':s2})
    df.to_csv('s1s2.csv')
    ax = df.plot()
    ax.plot(s1_t,s1_obs,'*')
    ax.plot(s2_t,s2_obs,'o')
    
    # ax.set_ylim(0,.4)
    ax.set_ylim(0,1)
    plt.grid()
    ax.set_xscale('log')
    ax.set_xlim([1e-2,1e8])
    plt.savefig('tmp.pdf')
    T1=0.0001
    dQ = WardLoughDepletion(0.0001, S1, K, lambd, x, y, t, NSteh1)
    
    df = pd.DataFrame(index=t, data={'dQ':dQ})
    df.to_csv('dQ.csv')
    ax = df.plot()
    ax.plot(q_t_0p0001,dQ_0p0001, 'o')
    # ax.set_ylim(0,.4)
    ax.set_ylim(0,1)
    plt.grid()
    ax.set_xscale('log')
    ax.set_xlim([1e-2,1e8])
    plt.savefig('tmp_Q.pdf')
    