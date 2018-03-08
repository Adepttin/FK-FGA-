import numpy as np

#Functions for calculating the second order contribution to dual self energy

def setUpExponent(nk):
#creates the coefficients exp(ikx) in the matrix exponent(nk,nk)
#or exp(-ikx) in the matrix exponent_neg(nk,nk)
    
    #exp(ikx)
    exponent=np.fromfunction(lambda i, k: 1j*(k-nk/2)*np.pi/(nk/2)*(i-nk/2), (nk,nk), dtype=complex)
    exponent=np.exp(exponent)
    
    #exp(-ikx)
    exponent_neg=np.fromfunction(lambda i, k: (-1)*1j*(k-nk/2)*np.pi/(nk/2)*(i-nk/2), (nk,nk), dtype=complex)
    exponent_neg=np.exp(exponent_neg)
    
    return (exponent,exponent_neg)

def readGdual(nk, nv):
#reads in "G0dual" from DMFT calculations and returns it as (nk,nk,2*nv) matrix
    
    Gdual = np.fromfile("G0dual", dtype = complex)
    Gdual = Gdual.reshape((nk,nk,2*nv))
    return Gdual

def readFloc(nv):
#reads in "LambdaUp" from DMFT calculations, slices it to the size of the nv-grid
#and returns it as (2*nv,2*nv) matrix

    F = np.fromfile("LambdaUp", dtype = complex)
    s = int(np.sqrt(len(F)))
    F=F.reshape((s,s))
    #slicing of center part of (2*nv,2*nv)
    F = F[(s/2-nv):(s/2+nv),(s/2-nv):(s/2+nv)]
    return F

def fourierTransform(space, momentum, nk, nv, exponent):
#Fourier Transformation of the matrix momentum(nk,nk,nv) to the matrix space
#using the exponential values of exponent(nk,nk)

    for v in range(0,2*nv):
        for x in range(0,nk):
            for y in range(0,nk):
                
                dsum = 0
                
                for kx in range(0,nk):
                    for ky in range(0,nk):
                        dsum += momentum[kx][ky][v]*exponent[x][kx]*exponent[y][ky]
                
                #normalisation
                dsum = dsum/float(nk*nk)
                space[x][y][v]=dsum
				
def calculateSigma(Sigma, F, Gx, Gx_neg, nk, nv, beta):
#calculates 2nd order diagram of Sigma(nk,nk,2*nv) (in real space) by using local vertex function F,
#dual Green's function G(r) Gx and G(-r) Gx_neg

    for v in range(0,2*nv):
        for x in range(0,nk):
            for y in range(0,nk):
                
                Sigmasum = 0
                
                #term for w = 0, sum over v'
                for vv in range(0,2*nv):
                    Sigmasum -= F[v,vv]*F[v,vv]*Gx[x][y][v]*Gx[x][y][vv]*Gx_neg[x][y][vv]
                    
                #term for v=v', sum over w
                for w in range(0,2*nv):
                    Sigmasum -= F[v,w]*F[v,w]*Gx[x][y][w]*Gx[x][y][v]*Gx_neg[x][y][w]
                
                #normalisation
                Sigmasum = Sigmasum/beta
                #Sigmasum = Sigmasum/float(nk*nk)
                Sigma[x][y][v] = Sigmasum
				
