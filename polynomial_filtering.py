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

from helpers import get_trace_bound
import numpy as np
from sage.all import *
import pyximport
pyximport.install()
from recursive import find_polynomials_recursive

MAXINT = 1000000 # just some big integer



UPPER = 0
LOWER = 1






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
    find_polynomials_recursive(len(coeffs), largest_root_bound, 
                                is_orientable, is_orientation_rev,
                               coeffs, traces, trace_bounds, 0,
                               good_coeffs, MAXINT if first_trace is None else first_trace)
    # good_polys = []
    # for coeffs in good_coeffs:
    #     poly = construct_poly(coeffs, is_orientable, is_orientation_rev)
    #     if is_good_poly(poly, largest_root_bound, is_orientable, CDF):
    #         if is_good_poly(poly, largest_root_bound, is_orientable, QQbar):
    #             pf_root = get_PF_root(poly, QQbar, is_orientable)
    #             good_polys.append((poly, pf_root))
    return sorted(good_coeffs, key = lambda x: x[1])


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






