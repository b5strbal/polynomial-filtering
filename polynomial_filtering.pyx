
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



cdef int cUPPER = 0
cdef int cLOWER = 1
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
                               first_trace=None):
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
                                                  i+1, is_orientable, typ)
    # print [trace_bounds[UPPER,i]-trace_bounds[LOWER,i]+1 for i in range(len(trace_bounds[UPPER]))]
    good_coeffs = []
    # count = np.zeros(1, dtype=np.int)
    find_polynomials_recursive(len(coeffs), is_orientable,
                               coeffs, traces, trace_bounds, 0,
                               good_coeffs, MAXINT if first_trace is None else first_trace)
    # print good_coeffs
    # print [-1,0,1,-1,-1,1,-1,0,1] in good_coeffs
    good_polys = []
    for coeffs in good_coeffs:
        poly = construct_poly(coeffs, is_orientable)
        if is_good_poly(poly, largest_root_bound, is_orientable, CDF):
            if is_good_poly(poly, largest_root_bound, is_orientable, QQbar):
                pf_root = get_PF_root(poly, QQbar, is_orientable)
                good_polys.append((poly, pf_root))
    # print len(good_polys)
    return sorted(good_polys, key = lambda x: x[1])


cdef long [:] backward_traces = np.zeros(30, dtype=np.int)
cdef long MAXINT = 1000000

cdef do_traces_pass_backwards(long[:] coeffs, long[:] trace_abs_bounds):
    cdef long [:] backward_coeffs = coeffs[-2::-1]

    # coeffs[-1] is the constant term of the polynomial. When it is negative,
    # then all coefficients should be taken with a negative sign.
    backward_traces[0] = -coeffs[-1]*backward_coeffs[0]
    if abs(backward_traces[0]) > trace_abs_bounds[0]:
        return False
    cdef int i, j
    for i in range(1, len(coeffs)-1):
        # using Newton's formula to figure out the next trace
        # E.g. p_4 = -a_1p_3 - a_2p_2 - a_3p_1 - 4a_4
        backward_traces[i] = -coeffs[-1] * (i+1) * backward_coeffs[i]
        for j in range(i):
            backward_traces[i] -= coeffs[-1] * backward_coeffs[j] * backward_traces[i-1-j]
        if abs(backward_traces[i]) > trace_abs_bounds[i]:
            return False
    return True


cdef find_polynomials_recursive(int num_coeffs,
                                int is_orientable,
                                long[:] coeffs,
                                long[:] traces,
                                long[:,:]trace_bounds,
                                int next_idx,
                                list good_coeffs,
                                long first_trace):

    """
    INPUT:

    - ``trace_bounds`` -- a list of 2-tuples, each tuple contaning a lower and
    upper bound for the trace. When the recursion starts, the length of this
    list is the number of middle coefficients (in the nonorientable case) and
    about half of the middle coeffs in the orientable case. As we get deeper in
    the recursion, we keep only the part of the tail that we haven't used yet.

    - ``largest_root_bound`` -- the upper bound for the largest root.

    - ``is_orientable`` -- whether the surface is orientable or not
    """
    cdef int j
    if next_idx == num_coeffs -1 and not is_orientable:
        for j in range(-1,2,2):
            coeffs[next_idx] = j
            find_polynomials_recursive(num_coeffs,
                                       is_orientable, coeffs, traces,
                                       trace_bounds, next_idx+1,
                                       good_coeffs, first_trace)
        return

    if next_idx == num_coeffs:
        if not is_orientable and \
           not do_traces_pass_backwards(coeffs, trace_bounds[cUPPER]):
            return
        good_coeffs.append(list(coeffs))
        return

    # next_idx - The index of the coefficient to figure out.
    cdef long trace_upper_bound = trace_bounds[cUPPER, next_idx]
    cdef long trace_lower_bound = trace_bounds[cLOWER, next_idx]

    # Trying to determine the next coefficient from the previous coefficients
    # and traces using Newton's formula, e.g.:
    # 4a_4 = -a_3p_1 - a_2p_2 - a_1p_3 - p_4
    # Here the coefficients are the a_i and the traces of the powers are p_1,
    # p_2, etc.
    coeffs[next_idx] = 0


    for j in range(next_idx):
        # The right hand side of the formula above, without p_4.
        coeffs[next_idx] -= traces[j] * coeffs[next_idx - 1 - j]

    # print 'b'
    # calculating the largest number less than the upper bound for the
    # trace (current_trace_bounds[1]) such that next_coeff minus this number is
    # divisible by n + 1
    traces[next_idx] = trace_upper_bound - (trace_upper_bound -
                                           coeffs[next_idx]) % (next_idx + 1)

    # Completing Newton's formula above to determine the next coefficient
    coeffs[next_idx] -= traces[next_idx]
    coeffs[next_idx] = coeffs[next_idx] / (next_idx + 1)


    while trace_lower_bound <= traces[next_idx] <= trace_upper_bound:
        if next_idx != 0 or first_trace == MAXINT or traces[next_idx] == first_trace:
            if not is_orientable and \
               2 * next_idx >= num_coeffs-1 and \
               coeffs[next_idx] % 2 != coeffs[num_coeffs-next_idx-2] % 2:
                # In the nonorientable case, if at least half of the coefficients
                # have been determined, then we make sure that the polynomial is
                # reciprocal mod 2.
                pass
            else:
                find_polynomials_recursive(num_coeffs,
                                           is_orientable, coeffs, traces,
                                           trace_bounds, next_idx+1,
                                           good_coeffs, first_trace)
        # Choosing a smaller trace and updating the coefficient
        traces[next_idx] -= next_idx + 1
        coeffs[next_idx] += 1

def get_trace_bound(degree, largest_root_bound, power, is_orientable,
                    upper_or_lower=UPPER):
    """
    INPUT:
    - ``degree`` -- the degree of the polynomial
    - ``largest_root_bound`` -- it is assumed that all roots of the polynomial
    are at most this bound in absolute value
    - ``power`` -- the power of the matrix whose trace is considered
    - ``is_orientable`` -- if True, then we consider orientable surfaces, if
    False, then nonorientable surfaces.
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
        if upper_or_lower == UPPER:
            bound = h*p+h/p
        else:
            bound = -(h-2)*p-(h-2)/p
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



def construct_poly(middle_coeffs, is_orientable):
    R = PolynomialRing(IntegerRing(), 'x')
    if is_orientable:
        coeffs = [1] + list(middle_coeffs) + list(middle_coeffs)[-2::-1] + [1]
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


upper_bounds = {3: 1.84, 4: 1.52, 5: 1.43, 6: 1.422, 7: 1.2885, 8:
                1.3568, 9: 1.2173, 10: 1.22262, 11: 1.1743, 12:1.2764, 13:
                1.145507, 14: 1.1875, 15:1.1249, 16:1.1426, 17:1.10939,
                18:1.20515, 19:1.097305, 21:1.087629}


def timed_test_by_first_trace(g, largest_root_bound):
    trace_bounds = [get_trace_bound(g, largest_root_bound, 1,
                                    False, typ) for
                    typ in [LOWER, UPPER]]

    import os
    traces = sorted(range(trace_bounds[0], trace_bounds[1] + 1), key = lambda x: abs(x))
    print traces
    for i in traces:
        # print 'sage polynomial_filtering.pyx ' + str(g) + ' --first_trace ' + str(i)
        os.system('sage polynomial_filtering.pyx ' + str(g) + ' --first_trace ' + str(i))



# ------------------------------------
import time

def timed_test(g, largest_root_bound, first_trace=None):
    """
    Tests and times a genus, and prints the results to a file and to the screen
    as well.
    """
    start = time.time()
    my_str = "\nTesting g=" + str(g) + " with limit on dilatation " +\
                         str(largest_root_bound) + '\n'
    if first_trace is not None:
        my_str += "Assuming that the first trace is " + str(first_trace) + "\n"
    result = find_polynomial_candidates(g,largest_root_bound,False,first_trace)

    for item in result:
        my_str += repr(item) + '\n'

    my_str += "Elapsed time:" + str(time.time() - start) + " seconds.\n"
    print my_str
    if first_trace is None:
        myfile = open("result_{0}.txt".format(g),'w')
    else:
        myfile = open("result_{0}_{1}.txt".format(g, first_trace),'w')
    myfile.write(my_str)
    myfile.close()


# -----------------------------------

import argparse

parser = argparse.ArgumentParser(description='Test polynomials.')
parser.add_argument('degree', type=int)
parser.add_argument('--largest_root_bound', nargs = '?', type=float)
parser.add_argument('--separate_by_first_trace', action='store_true', default=False)
parser.add_argument('--first_trace', nargs='?', type=int)
args = parser.parse_args()

print args

g = args.degree
if g is not None:
    largest_root_bound = args.largest_root_bound
    if largest_root_bound is None:
        largest_root_bound = upper_bounds[g]

    if args.separate_by_first_trace:
        timed_test_by_first_trace(g, largest_root_bound)
    elif args.first_trace is None:
        timed_test(g, largest_root_bound)
    else:
        timed_test(g, largest_root_bound, args.first_trace)
