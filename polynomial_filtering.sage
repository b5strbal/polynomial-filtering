
#*****************************************************************************
#       Copyright (C) 2014 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from itertools import product


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


def get_good_polys(trace_bounds, largest_root_bound, is_orientable,
                   first_coeffs = [],
                   first_traces = []):
    """


    INPUT:

    - ``trace_bounds`` -- a list of 2-tuples, each tuple contaning a lower and
    upper bound for the trace. When the recursion starts, the length of this
    list is the number of middle coefficients (in the nonorientable case) and
    about half of the middle coeffs in the orientable case. As we get deeper in
    the recursion, we keep only the part of the tail that we haven't used yet.

    - ``largest_root_bound`` -- the upper bound for the largest root.

    - ``is_orientable`` -- whether the surface is orientable or not

    - ``first_coeffs`` -- the middle coefficients that have been determined
    earlier in the resursion. Initially, this is empty, and it gets longer by
    one for each additional depth of the recursion.

    - ``first_traces`` -- the traces that have been determined
    earlier in the resursion. Initially, this is empty, and it gets longer by
    one for each additional depth of the recursion.

    """
    # There is no more coefficients to figure out, so we go ahead and construct
    # all the good polynomials from our list of possible coefficients.
    if len(trace_bounds) == 0:
        return good_polys_from_coeffs(first_coeffs,
                   largest_root_bound, is_orientable)

    # Not all coefficients have been determined, so we continue the recursion.

    good_polys = []
    current_trace_bounds = trace_bounds[0]

    # Trying to determine the next coefficient from the previous coefficients
    # and traces using Newton's formula, e.g.:
    # 4a_4 = -a_3p_1 - a_2p_2 - a_1p_3 - p_4
    # Here the coefficients are the a_i and the traces of the powers are p_1,
    # p_2, etc.
    next_coeff = 0
    n = len(first_coeffs)
    for j in range(n):
        # The right hand side of the formula above, without p_4.
        next_coeff -= first_traces[j] * first_coeffs[n - 1 - j]

    # calculating the largest number less than the upper bound for the
    # trace (current_trace_bounds[1]) such that next_coeff minus this number is
    # divisible by n + 1
    next_trace = current_trace_bounds[1] - (current_trace_bounds[1] - next_coeff) % (n + 1)

    # Completing Newton's formula above to determine the next coefficient
    next_coeff -= next_trace
    next_coeff /= n + 1
    while current_trace_bounds[0] <= next_trace <= current_trace_bounds[1]:
        if is_orientable or len(first_coeffs) <\
           len(trace_bounds) or\
           next_coeff % 2 == first_coeffs[len(trace_bounds)-1] % 2:
            good_polys.extend(get_good_polys(trace_bounds[1:],
                                             largest_root_bound,
                                             is_orientable,
                                             first_coeffs + [next_coeff],
                                             first_traces + [next_trace]))

        # Choosing a smaller trace and updating the coefficient
        next_trace -= n + 1
        next_coeff += 1
    return good_polys


def good_polys_from_coeffs(middle_coeffs, largest_root_bound, is_orientable):
    polys = construct_polys(middle_coeffs, is_orientable)
    good_polys = []
    for poly in polys:
        if not passes_newton_test(poly, largest_root_bound):
            continue
        pf_root = get_PF_root(poly, CDF, is_orientable)
        if pf_root == None:
            continue
        if not abs(pf_root) <= largest_root_bound:
            continue
        if pf_root == None or pf_root < 0:
            continue
        good_polys.append((poly, pf_root))
    return good_polys


def find_polynomial_candidates(degree, largest_root_bound,
                               is_orientable, first_trace=None):
    trace_bounds = get_trace_bounds(degree, largest_root_bound,
                               is_orientable)
    if first_trace != None:
        trace_bounds[0] = (first_trace,first_trace)
    # print "Number of possibilities for traces: ", \
    #         prod(tb[1] - tb[0] + 1 for tb in trace_bounds)

    good_polys = get_good_polys(trace_bounds, largest_root_bound,
                                is_orientable)

    return sorted(good_polys, key = lambda x: x[1])

def get_trace_bounds(degree, largest_root_bound, is_orientable):
    """
    Return the intervals of possible traces of the first few powers.

    INPUT:
    - ``degree`` -- the degree of the polynomial
    - ``largest_root_bound`` -- it is assumed that all roots of the polynomial
    are at most this bound in absolute value
    - ``is_orientable`` -- if True, then we consider orientable surfaces, if
    False, then nonorientable surfaces.

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
    trace_bounds = []
    if is_orientable:
        assert(degree % 2 == 0 and degree >= 4)
        h = degree/2
        for i in range(1, degree / 2 + 1):
            p = largest_root_bound^i
            trace_bounds.append((-(h-2)*p-(h-2)/p,h*p+h/p))
    else:
        h = degree // 2
        for i in range(1, degree):
            p = largest_root_bound^i
            trace_bounds.append(
                (min(-(h-2)*p-h/p - degree % 2,
                     2-2*h - degree%2),
                 h*p+h/p + degree % 2)
            )

    return [(floor(tb[0] +1),floor(tb[1])) for tb in trace_bounds]


def construct_polys(middle_coeffs, is_orientable):
    R.<x> = ZZ[]
    if is_orientable:
        coeffs_list = [[1] + middle_coeffs + middle_coeffs[-2::-1] + [1]]
    else:
        coeffs_list = [[1] + middle_coeffs + [1], [1] + middle_coeffs + [-1]]
    # coeffs_list.extend([[coeffs[k] * (-1)^k for k in
    #                      range(len(coeffs))] for coeffs in coeffs_list
    #                     if coeffs[1] > 0])
    print coeffs_list
    return [R(coeffs[::-1]) for coeffs in coeffs_list]





# def is_symmetric_mod2(polynomial):
#     p = polynomial
#     coeffs = p.coeffs()
#     d = p.degree()
#     return all((coeffs[i] - coeffs[d-i]) % 2 == 0
#                for i in range((d+1)//2))


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
                         key = lambda root_tuple : abs(root_tuple[0]))
    abs_roots = [abs(x[0]) for x in root_tuples]
    largest_root_tuple = root_tuples[-1]
    root = abs_roots[-1]

    # If ring = CDF, and there is floating point error, it might be
    # that a real root is returned as a complex root with small
    # imaginary part. That imaginary part should be too big.
    if abs(largest_root_tuple[0].imag()) >= TOLERANCE:
        print 1
        return None

    # If ring = QQbar, usually multiplicities are calculated
    # correctly. E.g., the output of (x^2-2*x+1).roots(ring=QQbar) is
    # [(1,2)]. In this case root_tuples is a list of roots WITHOUT
    # multiplicity, so we need to check that the multiplicity of the
    # largest root is 1.
    if largest_root_tuple[1] > 1:
        print 2
        return None

    # If ring = CDF, multiplicities are not
    # calculated correctly, so we have to check whether the largest
    # and second largest roots are too close.
    if root - abs_roots[-2] <= EPSILON:
        print 3
        return None

    if 1/root - abs_roots[0] >= EPSILON:
        print 4
        return None

    if not is_orientable and 1/root - abs_roots[0] >= -EPSILON:
        print 5
        return None

    return largest_root_tuple[0].real()


def to_latex_table(candidates):
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


#def find_min_dilatation_candidates(genus, upper_bound, is_orientable):
