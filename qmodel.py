from __future__ import division

from numpy import zeros, eye, linalg

def s_sum(i, j, N):
    return i * (N + 1) + j;

def complexity2throughput(N, q, a):
    M = N + 1
    L = zeros((M**2, M**2))
    for i in range(M):
        for j in range(M):
            Pc = q * i * j / N / N
            Qc = 1 - Pc
            #Pcc = Pc**2
            #Qcc = 1 - Pcc
            Pj = (N - j) / N
            Pi = (N - i) / N
            Qj = 1 - Pj
            Qi = 1 - Pi
            s = s_sum(i, j, N)
            if i > 0:
                L[s_sum(i - 1, j, N)][s] = a * Pc * Qi * Qj 
                if j < N:
                    L[s_sum(i - 1, j + 1, N)][s] = a * Pc * Qi * Pj
            
            if j > 0:
                L[s_sum(i, j - 1, N)][s] = (1 - a) * Pc * Qi * Qj 
                if i < N:
                    L[s_sum(i + 1, j - 1, N)][s] = (1 - a) * Pc * Pi * Qj

            if i < N:
                L[s_sum(i + 1, j, N)][s] = Qc * Pi * Qj + (1 - a) * Pc * Pi * Pj
                if j > 0:
                    L[s_sum(i + 1, j - 1, N)][s] = (1 - a) * Pc * Pi * Qj

            if j < N:
                L[s_sum(i, j + 1, N)][s] = Qc * Qi * Pj + a * Pc * Pi * Pj
                if i > 0:
                    L[s_sum(i - 1, j + 1, N)][s] = a * Pc * Qi * Pj

            if (i < N) and (j < N):
                L[s_sum(i + 1, j + 1, N)][s] = Qc * Pi * Pj

    for s in range(M*M):
        P = 0
        for l in range(M*M):
            P += L[l][s]
        L[s][s] = 1 - P
    B = zeros((M*M, 1))
    B[0] = 1.0
    I = eye(M*M)
    Z = I - L
    for i in range(M**2):
        Z[0][i] = 1.0
    P = linalg.solve(Z, B)
    Cp = 0.0
    Cm = 0.0
    Cr = 0.0
    for i in range(M):
        for j in range(M):
            Cp += i * P[s_sum(i, j, N)][0]
            Cm += j * P[s_sum(i, j, N)][0]
            Cr += (i + j) * P[s_sum(i, j, N)][0]

    return (Cr, Cp, Cm, P)

if __name__ == "__main__":
    N = 10
    Cr, Cp, Cm, P = complexity2throughput(N, 0.4, 0.5)
    print Cr/2/N


