import numpy as np
import math

def results_mtx(function, interval, resmat):
    """Print the Romberg results in the form Steps, Stepsize and Results"""
    i = j = 0
    print('Romberg integration from ' + str(interval[0]) + ' to ' + str(interval[1]))
    print
    for i in xrange(len(resmat)):
        print('Steps: ' + str((2**i)))
	print('StepSize: ' + str((interval[1]-interval[0])/(2.**i)))
	newresults = []
        for j in xrange(i+1):
	    newresults.append(resmat[i][j])
	print('Results: ' + str(newresults))
	print('')
    print('')
    print('The final result is ' + str(resmat[i][j]))
    print('after ' + str(2**(len(resmat)-1)+1) +' function evaluations.')


def romberg_diff(b, c, k):
    """
    Compute the differences for the Romberg quadrature corrections.
    """
    tmp = 4.0**k
    return (tmp * c - b)/(tmp - 1.0)


def difftrap(function, interval, numtraps):
    """
    Perform part of the trapezoidal rule to integrate a function.
    """
    if numtraps <= 0:
        raise ValueError("numtraps must be > 0 in difftrap().")
    elif numtraps == 1:
        return 0.5*(function(interval[0])+function(interval[1]))
    else:
        numtosum = numtraps/2
        h = float(interval[1]-interval[0])/numtosum
        lox = interval[0] + 0.5 * h
        points = lox + h * np.arange(numtosum)
        s = np.sum(function(points), axis=0)
        return s


def vectorize(func, args=(), vec_func=False):
    """Vectorize the call to a function.
    This is an internal utility function used by `romberg` 
    to create a vectorized version of a function.
    If `vec_func` is True, the function `func` is assumed to take vector
    arguments.
    """
    if vec_func:
        def vfunc(x):
            return func(x, *args)
    else:
        def vfunc(x):
            if np.isscalar(x):
                return func(x, *args)
            x = np.asarray(x)
            # call with first point to get output type
            y0 = func(x[0], *args)
            n = len(x)
            dtype = getattr(y0, 'dtype', type(y0))
            output = np.empty((n,), dtype=dtype)
            output[0] = y0
            for i in xrange(1, n):
                output[i] = func(x[i], *args)
            return output
    return vfunc


def romberg(f, a, b, tol=1.0e-6):
    """
	I, n_panels = romberg(f, a, b, tol=1.0e-6)
	Return the integral from a to b of f(x), by Romberg integration,
	as well as the number of panels used.
    """
    if np.isinf(a) or np.isinf(b):
        raise ValueError("Romberg integration only available "
                         "for finite limits.")
    show = True
    vec_func = False
    vfunc = vectorize(f, (), vec_func=vec_func)
    n = 1
    interval = [a, b]
    intrange = b - a
    ordsum = difftrap(vfunc, interval, n)
    result = intrange * ordsum
    resmat = [[result]]
    err = np.inf
    last_row = resmat[0]
    for i in xrange(1, 10+1):
        n *= 2
        ordsum += difftrap(vfunc, interval, n)
        row = [intrange * ordsum / n]
        for k in xrange(i):
            row.append(romberg_diff(last_row[k], row[k], k+1))
        result = row[i]
        lastresult = last_row[i-1]
        if show:
            resmat.append(row)
        err = abs(result - lastresult)
        if err < tol:
            break
        last_row = row

    if show:
        results_mtx(vfunc, interval, resmat)
    return result, str(2**(len(resmat)-1)+1)


if __name__ == "__main__":

	def myfunc(x):
		return (x**2)*math.pi
	print(romberg(myfunc, 0, 1))
