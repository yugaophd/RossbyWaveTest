def  VERT_FSFB3(N2_uniform, Pmid_uniform):
    
    import numpy as np
    from numpy import linalg as LA
    import scipy
    import gsw
    
    '''
    VERT_FSFB3.py - Yu Gao, April 18, 2022 
    This function solves the discretized wave projection problem, given the vertical profiles of Temperature, Salinity, Pressure and the depth inteval length. Note that it wants REGULARLY spaced depths. The seawater function sw_bfrq can be used to calculate N2.
    
    # gsw-python package: https://github.com/TEOS-10/GSW-Python.
    # gsw can be installed from a clone of the repo using
    # $ pip install gsw
    
    # Arguments:
    # T = temperature vector at same depths as salinity and pressure.
    # S = salinity vector at same depths as temperature and pressure.
    # P = pressure vector at same depths as temperature and salinity.
    # Dz = length of depth interval in meters.
    # Returns:
    # c2 = vector of square of the wavespeed.
    # Psi = matrix of eigenvectors (horizontal velocity and pressure structure functions).
    # G  =  matrix of integral of eigenvectors (vertical velocity structure functions).
    # N2 = Brunt-Vaisla frequency squared calculated at the midpoint pressures.
    # Pmid = midpoint pressures.
    '''

    
    import numpy as np
    import gsw
    from numpy import linalg as LA
    
    for i in range(len(N2_uniform)):
        if N2_uniform[i] < 0:
            N2_uniform[i] = np.min(np.absolute(N2_uniform))
   
    Dz = np.median(np.diff(Pmid_uniform))
    Dzrng = np.max(np.diff(Pmid_uniform)) - np.min(np.diff(Pmid_uniform))
    
    if(Dzrng > 1e-3 * Dz):
        print('Dz is not constant!!')
    # else: print('Dz is constant.')
    
    # add a point for the surface
    # note that this means that Psi has one more depth than N2!!
    M = len(N2_uniform) + 1
    #  Fill in D - the differential operator matrix.
    #  Surface (repeat N2 from midpoint depth)
    D = np.zeros([M, M]) # Make a M x M zero array
    
    D[0, 0] = -2/N2_uniform[0]
    D[0, 1] = 2/N2_uniform[0]

    # Interior
    for i in range(1, M-1):
        D[i, i-1] = 1/N2_uniform[i-1]
        D[i, i] = -1/N2_uniform[i-1] - 1/N2_uniform[i]
        D[i, i+1] = 1/N2_uniform[i]

    # Bottom
    #D[M-1, M-2] = 1/N2_uniform[M-2]
    D[M-1, M-2] = 2/N2_uniform[M-2]
    D[M-1, M-1]  = -2/N2_uniform[M-2]
    D = -1 * D / (Dz * Dz)
    
    # BDC: D should be a tridiagonal matrix
    # Calculate generalized eigenvalue problem
    Lambda, Psi = LA.eig(D)
    # LA.eig() Compute the eigenvalues and right eigenvectors of a square array.
    # Lambda:    The eigenvalues, each repeated according to its multiplicity.
    # Psi: The normalized (unit "length") eigenvectors, such that the 
    # column ``Psi[:,i]`` is the eigenvector corresponding to the eigenvalue Lambda[i].
    
    # Remove Barotropic Mode 
    # Lambda = 1/(radius of deformation * f)^2 (m/s)^{-2}
    # BT mode, rd = infinity. hence lambda = 0.
    ind = np.where(Lambda >= 1e-11)
    Lambda = Lambda[ind]
    Psi = Psi[:, ind].reshape([D.shape[0], len(ind[0])])
    
    # Sort eigenvalues and eigenvectors
    idx = np.argsort(Lambda)
    Lambda_sorted = Lambda[idx]
    Psi_sorted = Psi[:,idx]
    c2 = 1 / Lambda_sorted # Modal Speed
    
    # Reference eigencevtor matrix to top value
    for i in range(Psi_sorted.shape[1]):
        Psi_sorted[:, i] = Psi_sorted[:, i] / Psi_sorted[0, i]
    
    # Generates an integration matrix
    INT = INTEGRATOR(M, Dz)
    
    # 
    G = np.matmul(INT, Psi_sorted)

    return c2, Psi, G, N2_uniform, Pmid_uniform


def INTEGRATOR(M, Dz): 
    
    """
    Translated from Gabriel A. Vecchi - June 7, 1998 by Yu Gao, - April 18, 2022 
    This function Generates an integration matrix, 
    and Integrates from first point to each point.
    """
    
    import numpy as np
    
    INT = np.tril(np.ones(M))
    INT = INT - 0.5*(np.eye(M))
    INT[:,0] = INT[:,0] - 0.5;
    INT = INT * Dz
    
    return INT