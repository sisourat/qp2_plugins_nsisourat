   def tql2(n, d, e, V):
        #  This is derived from the Algol procedures tql2, by
        #  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        #  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        #  Fortran subroutine in EISPACK.
    
        num_opt = False  # using vectors from numpy make it faster

        if not num_opt:
            for i in range(1, n):  # (int i = 1; i < n; i++):
                e[i-1] = e[i]
        else:
            e[0:n-1] = e[1:n]
        e[n-1] = 0.0
    
        f = 0.0
        tst1 = 0.0
        eps = 2.0**-52.0
        for l in range(n):  # (int l = 0; l < n; l++) {

            # Find small subdiagonal element
     
            tst1 = max(tst1, abs(d[l]) + abs(e[l]))
            m = l
            while m < n: 
                if abs(e[m]) <= eps*tst1:
                    break
                m += 1
 
            # If m == l, d[l] is an eigenvalue,
            # otherwise, iterate.
     
            if m > l:
                iiter = 0
                while 1:  # do {
                    iiter += 1  # (Could check iteration count here.)
     
                    # Compute implicit shift
     
                    g = d[l]
                    p = (d[l+1] - g) / (2.0 * e[l])
                    r = (p**2 + 1)**0.5  # hypot(p, 1.0)
                    if p < 0:
                        r = -r
 
                    d[l] = e[l] / (p + r)
                    d[l+1] = e[l] * (p + r)
                    dl1 = d[l+1]
                    h = g - d[l]
                    if not num_opt: 
                        for i in range(l+2, n):
                            d[i] -= h
                    else:
                        d[l+2:n] -= h
 
                    f = f + h
     
                    # Implicit QL transformation.
     
                    p = d[m]
                    c = 1.0
                    c2 = c
                    c3 = c
                    el1 = e[l+1]
                    s = 0.0
                    s2 = 0.0
 
                    # hh = V.T[0].copy()  # only with num_opt
                    for i in range(m-1, l-1, -1):
                        # (int i = m-1; i >= l; i--) {
                        c3 = c2
                        c2 = c
                        s2 = s
                        g = c * e[i]
                        h = c * p
                        r = (p**2 + e[i]**2)**0.5  # hypot(p,e[i])
                        e[i+1] = s * r
                        s = e[i] / r
                        c = p / r
                        p = c * d[i] - s * g
                        d[i+1] = h + s * (c * g + s * d[i])
     
                        # Accumulate transformation.
     
                        if not num_opt:  # overall factor 3 in 30-D
                            for k in range(n):  # (int k = 0; k < n; k++){
                                h = V[k][i+1]
                                V[k][i+1] = s * V[k][i] + c * h
                                V[k][i] = c * V[k][i] - s * h
                        else:  # about 20% faster in 10-D
                            hh = V.T[i+1].copy()
                            # hh[:] = V.T[i+1][:]
                            V.T[i+1] = s * V.T[i] + c * hh
                            V.T[i] = c * V.T[i] - s * hh
                            # V.T[i] *= c
                            # V.T[i] -= s * hh
 
                    p = -s * s2 * c3 * el1 * e[l] / dl1
                    e[l] = s * p
                    d[l] = c * p
     
                    # Check for convergence.
                    if abs(e[l]) <= eps*tst1:
                        break
                # } while (Math.abs(e[l]) > eps*tst1);
 
            d[l] += f
            e[l] = 0.0

        # Sort eigenvalues and corresponding vectors.
        if 11 < 3:
            for i in range(n-1):  # (int i = 0; i < n-1; i++) {
                k = i
                p = d[i]
                for j in range(i+1, n):  # (int j = i+1; j < n; j++) {
                    if d[j] < p:  # NH find smallest k>i
                        k = j
                        p = d[j]

                if k != i: 
                    d[k] = d[i]  # swap k and i
                    d[i] = p   
                    for j in range(n):  # (int j = 0; j < n; j++) {
                        p = V[j][i]
                        V[j][i] = V[j][k]
                        V[j][k] = p
