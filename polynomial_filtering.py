#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
This file contains the code that filters polynomials.

The main method is ``find_polynomial_candidates``.
"""


import numpy as np
from sage.all import *
import pyximport
pyximport.install()
from recursive import find_polynomials_recursive

MAXINT = 1000000 # just some big integer

# Choosing too small TOLERANCE has the danger incorrectly eliminating
# good polynomials, for instance, because their largest root, which is
# real, is found to be a complex number with small imaginary part
# On the flip side, a big TOLERANCE might allow bad polynomials to
# creep in.
TOLERANCE = 1e-3

# Choosing too big EPSILON has the danger of incorrectly eliminating
# polynomials, for instance, because if the difference between the two
# largest root is found to be less than EPSILON, the polynomial is
# eliminated. With a very small ESPILON, bad polynomials might also accepted.
EPSILON = 1e-10

UPPER = 0
LOWER = 1


def passes_newton_test(polynomial, largest_root_bound):
    """
    Test if a polynomial passes the Newton test.

    The Newton test is the following. A monic polynomial with a positive Perron root R is positive, increasing and convex on `[R,\infty]`. So Newton's method is an efficient way to find the R if we start from a number in the interval `[R,\infty]`. The polynomial fails the test if at any point during the iteration
    - the value of the function is negative
    - the derivative is not positive
    - our point is very close to 1 or smaller than 1 (in the cases that we can reasonably test, the pA stretch factor cannot get much smaller than 1.1).

    The say that the polynomial passes the test if a point and its iterate are very close to each other. If this does not happen in ``MAX_ITERATIONS`` steps, then the test also fails.

    INPUT:

    - ``polynomial`` -- a polynomial
    - ``largest_root_bound`` -- an upper bound for the Perron root

    EXAMPLES:

    sage: from polynomial_filtering import passes_newton_test
    sage: passes_newton_test(x^3-x^2-x-1, 2)
    True
    sage: passes_newton_test(x^3-1, 2)
    False

    """
    # by Penner's lower bound for the dilatation for
    # pseudo-Anosovs, dilatations smaller than 1.001 may occur
    # only if the genus is larger than 50, so it is a safe choice.
    LOWER_BOUND_FOR_PERRON_ROOT = 1.001

    # the number of iterations before raising an error. The test
    # should always finish in just a few iterations.
    MAX_ITERATIONS = 100

    derivative = polynomial.derivative()
    X = largest_root_bound
    for i in xrange(MAX_ITERATIONS):
        prev_X = X

        # the function value must be strictly positive
        y = polynomial(x=X)
        if y < 0:
            return False

        # the slope must be strictly positive (even at the largest root)
        m = derivative(x=X)
        if m <= EPSILON:
            return False

        X = X - y/m
        if X < LOWER_BOUND_FOR_PERRON_ROOT:
            return False
        # print x, y, m
        # when two subsequent x-values get too close to each other
        # without failing before, the polynomial passes
        if abs(X - prev_X) < 0.0001:
            return True

    # generally the test finishes in just a few iterations, so
    # this should not happen.
    assert(False)


def is_good_poly(poly, largest_root_bound, is_orientable, ring=CDF):
    """
    Decides if a polynomial satisfies some necessary conditions. These tests are:
    - the Newton test (see ``passes_newton_test``)
    - there should be Perron root (see ``get_PF_root``)
    - the Perron root should be positive
    - the Perron root should be less than ``largest_root_bound``

    INPUT:

    - ``poly`` -- the polynomial to be tested
    - ``largest_root_bound`` -- the roots of the polynomial should be less than this
    - ``is_orientable`` -- whether the surface is orientable or not (the effect of this is whether the reciprocal of the Perron root is allowed to also be a root of not, see ``get_PF_root``)
    - ``ring`` (default: CDF)-- the ring over which the root is computed. CDF is fast but inaccurate, QQbar is exact but slow.

    EXAMPLES:

        sage: from polynomial_filtering import is_good_poly
        sage: is_good_poly(x^3-x^2-x-1, 1.84, False, CDF)
        True
        sage: is_good_poly(x^3-x^2-x, 1.84, False, CDF)
        False

    """
    if not passes_newton_test(poly, largest_root_bound):
        return False
    pf_root = get_PF_root(poly, ring, is_orientable)
    if pf_root == None or pf_root < 0:
        return False
    if not abs(pf_root) <= largest_root_bound:
        return False
    return True



def get_trace_bound(degree, largest_root_bound, power, is_orientable,
                    is_orientation_rev=False,
                    upper_or_lower=UPPER,
                    optimized_bounds=True):
    """
    Return an upper or lower bound for the trace of the power of a matrix whose eigenvalues are bounded by the specified bound.

    INPUT:
    - ``degree`` -- the degree of the polynomial
    - ``largest_root_bound`` -- it is assumed that all roots of the polynomial
    are at most this bound in absolute value
    - ``power`` -- the power of the matrix whose trace is considered
    - ``is_orientable`` -- if True, then we consider orientable surfaces, if
    False, then nonorientable surfaces.
    - ``is_orientation_rev`` -- if True, the orientation of the orientable surface is reversed, otherwise it is preserved.
    - ``upper_or_lower`` -- whether an upper of lower bound is returned

    EXAMPLES:

        sage: from polynomial_filtering import get_trace_bound, UPPER, LOWER
        sage: get_trace_bound(3, 1.82, 1, False, False, UPPER)
        3
        sage: get_trace_bound(3, 1.82, 1, False, False, LOWER)
        0
        sage: get_trace_bound(3, 1.82, 4, False, False, UPPER)
        12
        sage: get_trace_bound(3, 1.82, 4, False, False, LOWER)
        0


    REFERENCES:

    For the bounds on orientable surfaces, see Lemma A.1 in

    ```
    Lanneau, Thiffeault: On the minimum dilatation of pseudo-Anosov
    homeromorphisms on surfaces of small genus (2011)
    ```

    For nonorientable surfaces, see Section 4.2 of

    ```
    Liechti, Strenner: Minimal pseudo-Anosov stretch factors for nonorientable surfaces and orientation-reversing maps
    ```

    REMARKS:

    In the nonorientable case, there isn't a huge difference between using carefully optimized bounds and trivial bounds. Here are some running times.

    Degree  Time (optimized bounds)     Time (trivial bounds)
    3       0.041198968887              0.0217459201813
    5       0.031574010849              0.0337769985199
    7       0.0774021148682             0.116263866425
    9       0.225617170334              0.426881074905
    10      1.67028188705               3.74016213417
    11      1.58136796951               4.64646601677
    13      21.6418228149               93.6907081604
    15      561.597373009               2443.97552204
    """
    p = largest_root_bound**power
    h = degree // 2
    if is_orientable:
        assert(degree % 2 == 0 and degree >= 4)
        if not is_orientation_rev:
            if upper_or_lower == UPPER:
                bound = h*p+h/p
            else:
                bound = -(h-2)*p-(h-2)/p
        else:
            # These bound can be improved, but for safety, we use the trivial bounds.
            if upper_or_lower == UPPER:
                bound = h*p+h/p
            else:
                bound = -h*p-h/p
    else:
        if optimized_bounds:
            if upper_or_lower == UPPER:
                bound = h*p+h/p + degree % 2
            else:
                bound = min(-(h-2)*p-h/p - degree % 2,
                    2-2*h - degree%2)
        else:
            if upper_or_lower == UPPER:
                bound = h*p+h/p + degree % 2
            else:
                bound = -h*p-h/p - degree % 2
    if upper_or_lower == UPPER:
        return floor(bound)
    else:
        return floor(bound + 1)


class MiddleCoeffError(Exception):
    pass


def construct_poly(middle_coeffs, is_orientable, is_orientation_rev=False):
    """
    Construct a monic polynomial (sometimes reciprocal) with sufficiently many middle coefficients.

    For orientation-preserving maps, a reciprocal polynomial is constructed. For orientation-reversing maps, a resulting polynomial has the property the the terms with odd powers for a reciprocal polynomial and the terms with even powers are reciprocal up to negation. This condition forces the middle coefficient of the resulting polynomial to be zero if it has even power. If this is not satisfied, a MiddleCoeffError is returned.

    INPUT:

    - ``middle_coeffs`` -- the middle coefficients of the polynomial
    - ``is_orientable`` -- whether the surface is orientable
    - ``is_orientation_rev`` (default: False)-- whether the map is orientation-reversing (only has an effect for orientable surfaces).

    EXAMPLES:

        sage: from polynomial_filtering import construct_poly
        sage: construct_poly([5,6,7],False)
        x^3 + 5*x^2 + 6*x + 7

        sage: construct_poly([5,6,7],True,False)
        x^6 + 5*x^5 + 6*x^4 + 7*x^3 + 6*x^2 + 5*x + 1

        sage: construct_poly([5,6,7],True,True)
        x^6 + 5*x^5 + 6*x^4 + 7*x^3 - 6*x^2 + 5*x - 1

        sage: construct_poly([5,6,7,8],True,True)
        Traceback (most recent call last)
        ...
        MiddleCoeffError: The middle coefficient has to be zero if it has even power and the map is orientation-reversing.

    """
    R = PolynomialRing(IntegerRing(), 'x')
    if is_orientable:
        coeffs = [1] + list(middle_coeffs) + list(middle_coeffs)[-2::-1] + [1]
        if is_orientation_rev:
            if len(middle_coeffs) % 2 == 0:
                # The power of the middle coefficient is even, so
                # the middle coefficient has to be zero.
                if middle_coeffs[-1] != 0:
                    raise MiddleCoeffError("The middle coefficient has to be zero if it has even power and the map is orientation-reversing.")
            n = len(coeffs)
            for i in range(n-1, n//2, -2):
                coeffs[i] *= -1
    else:
        coeffs = [1] + list(middle_coeffs)
    return R(coeffs[::-1])




def get_PF_root(poly, ring, is_orientable):
    r"""Return the Perron root of the polynomial if is exists.

    Let `r_1,\ldots,r_n` be the roots with multiplicity of ``poly``,
    sorted increasingly accoring to their absolute values, i.e. `r_n`
    is the largest root in absolute value. The polynomial is good if
    all of the following are true:

    - `r_n` is real
    - `r_n` is a singe root
    - `|r_n|>|r_{n-1}|`
    - `|r_1| >= 1/|r_n|` (with strict inequality in the nonorientable case)

    INPUT:

    - ``poly`` -- a polynomial with integer coefficients

    - ``ring`` -- the ring in which roots are calculated. ``CDF`` is
      for fast approximate complex roots, ``QQbar``, which is much
      slower, is for exact algebraic numbers.

    - ``is_orientable`` -- whether the polynomial is a candidate for
      the action on homology of a pseudo-Anosov on an orientable or
      non-orientable surface. The difference is that in the orientable
      case |r_1|=1/|r_n| is okay, while in the nonorientable case
      |r_1| > 1/|r_n| must hold.

    OUTPUT:

    the Perron root, if exists, in the specified ``ring``. If it does
    not exists, ``None`` is returned.

    EXAMPLES:

    This one passes the tests::

        sage: from polynomial_filtering import get_PF_root
        sage: get_PF_root(x^3-x^2-x-1, CDF, False)
        1.83928675521

        sage: get_PF_root(x^3-x^2-x-1, QQbar, False)
        1.839286755214161?

    The Perron root might be negative:

        sage: get_PF_root(x^3+x^2-x+1, CDF, False)
        -1.83928675521

    This one passes the tests in the orientable case, but not in the
    non-orientable one::

        sage: get_PF_root(x^4-x^3-x^2-x+1, CDF, True)
        1.72208380574

        sage: get_PF_root(x^4-x^3-x^2-x+1, CDF, False)

    """
    root_tuples = sorted(poly.roots(ring=ring),
                         key = lambda root_tuple : CDF(abs(root_tuple[0])))
    abs_roots = [abs(x[0]) for x in root_tuples]
    largest_root_tuple = root_tuples[-1]
    root = abs_roots[-1]

    # If ring = CDF, and there is floating point error, it might be
    # that a real root is returned as a complex root with small
    # imaginary part. That imaginary part should be too big.
    if abs(largest_root_tuple[0].imag()) >= TOLERANCE:
        # print 1
        return None

    # If ring = QQbar, usually multiplicities are calculated
    # correctly. E.g., the output of (x^2-2*x+1).roots(ring=QQbar) is
    # [(1,2)]. In this case root_tuples is a list of roots WITHOUT
    # multiplicity, so we need to check that the multiplicity of the
    # largest root is 1.
    if largest_root_tuple[1] > 1:
        # print 2
        return None

    # If ring = CDF, multiplicities are not
    # calculated correctly, so we have to check whether the largest
    # and second largest roots are too close.
    if root - abs_roots[-2] <= EPSILON:
        # print 3, abs_roots
        return None

    if 1/root - abs_roots[0] >= EPSILON:
        # print 4, poly
        return None

    if not is_orientable and 1/root - abs_roots[0] >= -EPSILON:
        # print 5
        return None

    return largest_root_tuple[0].real()


def find_polynomial_candidates(degree, largest_root_bound, is_orientable,
                               is_orientation_rev=False, first_trace=None):
    """
    Search for polynomials with Perron root less than the specified bound.

    INPUT:

    - ``degree`` -- the degree of the polynomial
    - ``largest_root_bound`` -- the upper bound on the Perron root
    - ``is_orientable`` -- whether the surface is orientable
    - ``is_orientation_rev`` -- whether the map is orientation-reversing (when the surface is orientable)
    - ``first_trace`` (default: None) -- if None, then the possible traces are determined automatically. If specified, then only polynomial with the specified first trace are searched. This is useful for parallelizing searches.

    EXAMPLE:

    sage: from polynomial_filtering import find_polynomial_candidates
    sage: find_polynomial_candidates(3, 2, False)
    [(x^3 - x^2 - x - 1, 1.839286755214161?)]

    sage: find_polynomial_candidates(4, 2, True)
    [(x^4 - x^3 - x^2 - x + 1, 1.722083805739043?),
     (x^4 - 2*x^3 + x^2 - 2*x + 1, 1.883203505913526?)]

    """
    # The number of coefficients to store.
    num_coeffs = degree if not is_orientable else degree // 2

    # For a polynomial x^d+a_{d-1}x^{d-1}+...+a_0$, this is the sequence
    # a_{d-1}, ..., a_0.
    coeffs = np.zeros(num_coeffs, dtype=np.int)

    # The traces p_1 = z_1+..._z_d, p_2 = z_1^2+...+z_d^2, etc.
    traces = np.zeros(num_coeffs, dtype=np.int)

    trace_bounds = np.zeros((2,num_coeffs), dtype=np.int)
    for typ in [LOWER, UPPER]:
        for i in range(num_coeffs):
            trace_bounds[typ,i] = get_trace_bound(degree, largest_root_bound,
                                                  i+1, is_orientable, 
                                                  is_orientation_rev, typ)
    good_coeffs = []
    find_polynomials_recursive(len(coeffs), is_orientable, is_orientation_rev,
                               coeffs, traces, trace_bounds, 0,
                               good_coeffs, MAXINT if first_trace is None else first_trace)
    good_polys = []
    for coeffs in good_coeffs:
        poly = construct_poly(coeffs, is_orientable, is_orientation_rev)
        if is_good_poly(poly, largest_root_bound, is_orientable, CDF):
            if is_good_poly(poly, largest_root_bound, is_orientable, QQbar):
                pf_root = get_PF_root(poly, QQbar, is_orientable)
                good_polys.append((poly, pf_root))
    return sorted(good_polys, key = lambda x: x[1])


def to_latex_table(candidates):
    r"""
    Create a latex table of stretch factors and polynomials.

    INPUT:

    - ``candidates`` -- a list of tuples (polynomial, stretch factor)

    EXAMPLE:

        sage: load('data.sage')
        sage: from polynomial_filtering import to_latex_table
        sage: to_latex_table(results['nonor'][3])

        \begin{tabular}{c|c}
        Largest root & Characteristic polynomial of action on cohomology \\ \hline
        1.83928675521416 & $ x^{3} - x^{2} - x - 1 $ \\ \hline
        \end{tabular}

    """
    s = "\\begin{tabular}{c|c}\n"
    s += "Largest root & Characteristic polynomial of action "\
         "on cohomology \\\\ \\hline\n"
    for cand in candidates:
        poly = cand[0]
        s += latex(cand[1]) + " & $" + latex(poly)
        # if not poly.is_irreducible():
            # s += "=" + latex(poly.factor())
        s += "$ \\\\ \\hline \n"
    s += "\\end{tabular}"
    return s

# def roots(poly):
#     return [x[0] for x in poly.roots(ring=QQbar)]






