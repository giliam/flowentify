#!/usr/bin/python
# coding:utf-8
from itertools import *
DEBUG = False

etatsPossibles = ["bourgogne","vin","bordeaux","cheval","hebergement","alsace"]
PI = [0.2,0.1,0.2,0.1,0.2,0.2]
A = [
        [0.01,  0.6,    0.01,   0.1,     0.27,    0.01],
        [0.4,   0.01,   0.4,    0.04,    0.05,    0.1],
        [0.01,  0.6,    0.01,   0.12,    0.25,    0.01],
        [0.3,   0.05,   0.4,    0.01,    0.14,    0.1],
        [0.3,   0.04,   0.3,    0.05,    0.01,    0.3],
        [0.01,  0.6,    0.01,   0.17,    0.2,     0.01],
    ]
visibles = ["degustation","hebergement","balade"]
E = [
        [0.5,   0.25,   0.25],
        [0.87,   0.12,   0.01],
        [0.5,   0.25,   0.25],
        [0.06,   0.14,   0.8],
        [0.03,   0.9,   0.07],
        [0.5,   0.25,   0.25],
    ]
N = len(A)
S = len(visibles)
etats = range(N)


def argmax(values):
    return max(enumerate(values), key=lambda x:x[1])

def forward (PI, A, E, sequence):
    T = len(sequence)
    temps = range(T)

    alpha = [[0. for t in temps] for i in etats]
    for i in etats:
        alpha[i][0] = PI[i] * E[i][sequence[0]]

    for t in range(1,T):
        for j in etats:
            for i in etats:
                alpha[j][t] += alpha[i][t-1]*A[i][j]
            alpha[j][t] *= E[j][sequence[t]]
    proba = 0.
    for j in etats:
        proba += alpha[j][T-1]
    if DEBUG:
        print "---------------"
        print "FORWARD"
        print "Pour la séquence : "
        print sequence
        print "Proba de : "
        print proba
        print "---------------\n\n"
    return proba,alpha


def backward (PI, A, E, sequence):
    T = len(sequence)
    temps = range(T)

    beta = [[0. for t in temps] for i in etats]
    for i in etats:
        beta[i][T-1] = 1

    for t in range(T-2, -1, -1):
        for j in etats:
            r = 0
            for i in etats:
                r += beta[i][t+1]*A[j][i]*E[i][sequence[t+1]]
            beta[j][t]
    proba = 0.
    for i in etats:
        proba += PI[i]*E[i][sequence[0]]*beta[i][0]
    if DEBUG:
        print "---------------"
        print "BACKWARD"
        print "Pour la séquence : "
        print sequence
        print "Proba de : "
        print proba
        print "---------------\n\n"
    return beta

def viterbi(PI,A,E,sequence):
    T = len(sequence)
    temps = range(T)

    delta = [[0. for t in temps] for i in etats]
    for i in etats:
        delta[i][0] = PI[i] * E[i][sequence[0]]
    psi = [[0. for t in temps] for i in etats]
    for t in range(0,T-1):
        for j in etats:
            for i in etats:
                psi[j][t],delta[j][t+1] = argmax([delta[i][t]*A[i][j] for i in etats])
                delta[j][t+1] *= E[j][sequence[t+1]]
    path = [0. for j in temps]
    proba = 1
    for t in range(T-1, -1, -1):
        path[t],p = argmax([delta[j][t] for j in etats])
        proba *= p
    return path,proba

def baum_welch(PI,A,E,sequence):
    T = len(sequence)
    temps = range(T)

    xi = [[[0. for t in temps] for i in etats] for j in etats]
    gamma = [[0. for i in etats] for t in temps]
    proba,alpha = forward(PI,A,E,sequence)
    beta = backward(PI,A,E,sequence)
    for t in range(T-1):
        for i in etats:
            s = 0.
            for j in etats:
                xi[i][j][t] = alpha[i][t]*A[i][j]*E[j][sequence[t+1]]*beta[j][t+1]/proba
                s += xi[i][j][t]
            gamma[t][i] = s

    

    nPI = gamma[0]
    nA = [[0. for i in etats] for j in etats]
    nE = [[E[j][i] for i in range(S)] for j in etats]

    sGamma = [0. for i in etats]
    for i in etats:
        for t in range(T-1):
            sGamma[i] += gamma[t][i]
        for j in etats:
            sXi = 0.
            for t in range(T-1):
                sXi += xi[i][j][t]
            nA[i][j] = sXi/sGamma[i]
    print gamma
    for j in etats:
        for o in range(S):
            s = 0.
            print o
            for t in range(T):
                print t
                if o == sequence[t]:
                    s += gamma[t][j]
            nE[j][o] = s/(sGamma[j]+gamma[T-1][j])
            print "--"
    print E
    print nE
    nProba,nAlpha = forward(nPI,nA,nE,sequence)

    print "---------------"
    print "BAUM WELCH"
    print "Pour la séquence : "
    print sequence
    print "AVANT"
    print "Proba de : "
    print proba
    print "APRES"
    print "Proba de : "
    print nProba
    print "---------------\n\n"
if __name__ == "__main__":
    etatsPossibles = ["bourgogne","vin","bordeaux","cheval","hebergement","alsace"]
    PI = [0.2,0.1,0.2,0.1,0.2,0.2]
    A = [
            [0.01,  0.6,    0.01,   0.1,     0.27,    0.01],
            [0.4,   0.01,   0.4,    0.04,    0.05,    0.1],
            [0.01,  0.6,    0.01,   0.12,    0.25,    0.01],
            [0.3,   0.05,   0.4,    0.01,    0.14,    0.1],
            [0.3,   0.04,   0.3,    0.05,    0.01,    0.3],
            [0.01,  0.6,    0.01,   0.17,    0.2,     0.01],
        ]
    visibles = ["degustation","hebergement","balade"]
    E = [
            [0.4,   0.35,   0.25],
            [0.87,   0.12,   0.01],
            [0.5,   0.25,   0.25],
            [0.06,   0.14,   0.8],
            [0.03,   0.9,   0.07],
            [0.5,   0.25,   0.25],
        ]

    bf = False
    vi = False
    bw = True
    if bf:
        backward(PI,A,E,[0])
        backward(PI,A,E,[0,1])
        backward(PI,A,E,[1])
        backward(PI,A,E,[2])
        forward(PI,A,E,[0])
        forward(PI,A,E,[0,1])
        forward(PI,A,E,[1])
        forward(PI,A,E,[2])
    if vi:
        viterbi(PI,A,E,[0,1])
    if bw:
        baum_welch(PI,A,E,[0,1])