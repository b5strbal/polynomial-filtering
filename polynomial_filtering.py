
#*****************************************************************************
#       Copyright (C) 2017 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



import numpy as np
from sage.all import *
import pyximport
pyximport.install()
from recursive import find_polynomials_recursive
MAXINT = 1000000

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
    # by Penner's lower bound for the dilatation for
    # pseudo-Anosovs, dilatations smaller than 1.001 may occur
    # only if the genus is larger than 50, so it is a safe choice.
    LOWER_BOUND_FOR_PERRON_ROOT = 1.001

    # the number of iterations before raising an error. The test
    # should always finish in just a few iterations.
    MAX_ITERATIONS = 100

    derivative = polynomial.derivative()
    x = largest_root_bound
    for i in xrange(MAX_ITERATIONS):
        prev_x = x

        # the function value must be strictly positive
        y = polynomial(x)
        if y < 0:
            return False

        # the slope must be strictly positive (even at the largest root)
        m = derivative(x)
        if m <= EPSILON:
            return False

        x = x - y/m
        if x < LOWER_BOUND_FOR_PERRON_ROOT:
            return False
        # print x, y, m
        # when two subsequent x-values get too close to each other
        # without failing before, the polynomial passes
        if abs(x - prev_x) < 0.0001:
            return True

    # generally the test finishes in just a few iterations, so
    # this should not happen.
    assert(False)


def is_good_poly(poly, largest_root_bound, is_orientable, ring=CDF):
    if not passes_newton_test(poly, largest_root_bound):
        return False
    pf_root = get_PF_root(poly, ring, is_orientable)
    if pf_root == None:
        return False
    if not abs(pf_root) <= largest_root_bound:
        return False
    if pf_root == None or pf_root < 0:
        return False
    return True


def find_polynomial_candidates(degree, largest_root_bound, is_orientable,
                               is_orientation_rev=False, first_trace=None):
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
    # print [trace_bounds[UPPER,i]-trace_bounds[LOWER,i]+1 for i in range(len(trace_bounds[UPPER]))]
    good_coeffs = []
    # count = np.zeros(1, dtype=np.int)
    find_polynomials_recursive(len(coeffs), is_orientable, is_orientation_rev,
                               coeffs, traces, trace_bounds, 0,
                               good_coeffs, MAXINT if first_trace is None else first_trace)
    # print good_coeffs
    # print [-1,0,1,-1,-1,1,-1,0,1] in good_coeffs
    good_polys = []
    for coeffs in good_coeffs:
        poly = construct_poly(coeffs, is_orientable, is_orientation_rev)
        if is_good_poly(poly, largest_root_bound, is_orientable, CDF):
            if is_good_poly(poly, largest_root_bound, is_orientable, QQbar):
                pf_root = get_PF_root(poly, QQbar, is_orientable)
                good_polys.append((poly, pf_root))
    # print len(good_polys)
    return sorted(good_polys, key = lambda x: x[1])



def get_trace_bound(degree, largest_root_bound, power, is_orientable,
                    is_orientation_rev=False,
                    upper_or_lower=UPPER):
    """
    INPUT:
    - ``degree`` -- the degree of the polynomial
    - ``largest_root_bound`` -- it is assumed that all roots of the polynomial
    are at most this bound in absolute value
    - ``power`` -- the power of the matrix whose trace is considered
    - ``is_orientable`` -- if True, then we consider orientable surfaces, if
    False, then nonorientable surfaces.
    - ``is_orientation_rev`` -- if True, the orientation of the orientable surface is reversed, otherwise it is preserved.
    - ``upper_or_lower`` -- whether an upper of lower bound is returned

    REFERENCES:

    For the bounds on orientable surfaces, see Lemma A.1 in

    ```
    Lanneau, Thiffeault: On the minimum dilatation of pseudo-Anosov
    homeromorphisms on surfaces of small genus (2011)
    ```

    For nonorientable surfaces, see Section 4.2 of

    ```
    Liechti, Strenner: Minimal penner dilatations for nonorientable surfaces
    ```
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
        if upper_or_lower == UPPER:
            bound = h*p+h/p + degree % 2
        else:
            bound = min(-(h-2)*p-h/p - degree % 2,
                 2-2*h - degree%2)
    if upper_or_lower == UPPER:
        return floor(bound)
    else:
        return floor(bound + 1)


# class MiddleCoeffError(Exception):
#     pass


def construct_poly(middle_coeffs, is_orientable, is_orientation_rev=False):
    R = PolynomialRing(IntegerRing(), 'x')
    if is_orientable:
        coeffs = [1] + list(middle_coeffs) + list(middle_coeffs)[-2::-1] + [1]
        if is_orientation_rev:
            # if len(middle_coeffs) % 2 == 0:
            #     # The power of the middle coefficient is even, so
            #     # the middle coefficient has to be zero.
            #     if middle_coeffs[-1] != 0:
            #         raise MiddleCoeffError
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


def to_latex_table(candidates):
    """
    Creates a latex table of stretch factors and polynomials

    INPUT:

    - ``candidates`` -- a list of tuples (polynomial, stretch factor)

    """
    s = "\\begin{tabular}{c|c}\n"
    s += "Largest root & Characteristic polynomial of action "\
         "on cohomology \\\\ \\hline\n"
    for cand in candidates:
        poly = cand[0]
        s += latex(cand[1]) + " & $" + latex(poly)
        if not poly.is_irreducible():
            s += "=" + latex(poly.factor())
        s += "$ \\\\ \\hline \n"
    s += "\\end{tabular}"
    return s

def roots(poly):
    return [x[0] for x in poly.roots(ring=QQbar)]






