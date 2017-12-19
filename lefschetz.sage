import itertools
import numpy as np
from collections import namedtuple

x = polygen(ZZ, 'x')

poly_candidates = {
    7: [
        x^7 - x^4 - x^3 - 1
    ],

    8: [
    x^8 - x^7 + x^6 - x^5 - x^4 + x^3 - x^2 + x - 1, # LEF 32
    x^8 + x^7 - x^5 - 2*x^4 - x^3 - x - 1,
    x^8 - x^7 - x^5 + x^3 - x + 1, # 7 4-pronged singularities that are
        # permuted. Multiplying by (x^6+x^5+x^4+x^3+x^2+x+1), it is likely that
        # the Teichmuller polynomial specializes to x^14 - x^11 - x^10 - 2*x^7
        # + x^4 + x^3 + 1. This is also (x - 1) * (x^7 - x^4 - x^3 - 1).
    x^8 + 2*x^7 + x^6 - x^5 - 2*x^4 - 3*x^3 - 3*x^2 - 2*x - 1, # LEF 28
    x^8 - x^6 - x^5 - x^3 + x^2 + 1,  # LEF 24
    x^8 + x^7 - x^6 - x^5 - x^3 - x^2 - x - 1,
    x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1, # 7 4-pronged singularities that are
        # permuted. Multiplying by (x^6+x^5+x^4+x^3+x^2+x+1), it is likely that
        # the Teichmuller polynomial specializes to x^14 - x^12 - x^9 - 2*x^7 +
        # x^5 + x^2 + 1. This is also (x - 1) * (x^7 - x^5 - x^2 - 1)
    x^8 - 2*x^5 - 1, # LEF 28
    x^8 - x^6 - x^5 - x^4 - x^3 + x^2 + 2*x + 1,
    x^8 + x^7 - 2*x^3 - 4*x^2 - 3*x - 1, # LEF 18
    x^8 - x^7 - 2*x^5 + 2*x^4 - x + 1, # LEF 22
    x^8 - x^7 - 2*x^3 + x + 1, # LEF 14 --- 14 fixed 4-pronged singularities
    # for f^9, no regular fixed points, pretty sure it's impossible
    x^8 - x^6 - x^5 - x^4 + x^3 + x^2 - 1, # LEF 14 --- same
    x^8 + x^7 - 2*x^5 - 2*x^4 - x - 1, # LEF 26
    x^8 - x^6 - x^5 + x^3 - x^2 - 1, # LEF 28
    x^8 + x^7 - x^5 - 3*x^4 - x^3 - x - 1, # LEF 22
    x^8 + x^6 - x^4 - 2*x^3 - 3*x^2 - 2*x - 1, # LEF 18
    x^8 - x^5 - x^4 - x^3 - 1],

    10: [
        x^10 - x^7 - x^6 - x^5 + x^4 + x^3 - 1,
        x^10 - x^9 + x^8 - x^7 - x^5 + x^3 - x^2 + x - 1,
        x^10 + x^9 - x^6 - 2*x^5 - x^4 - x - 1,
        x^10 - x^9 - x^6 + x^4 - x + 1,
        x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1
    ],

    12: [x^12 - x^11 + x^10 - x^9 + x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^11 - x^7 + x^5 - x + 1,
         x^12 + x^11 - x^7 - 2*x^6 - x^5 - x - 1,
         x^12 + x^11 - x^8 - x^7 - x^5 - x^4 - x - 1,
         x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - x + 1,
         x^12 + 2*x^11 + 2*x^10 + x^9 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 3*x^3 - 2*x^2 - 2*x - 1,
         x^12 - x^9 - x^7 - x^5 + x^3 + 1,
         x^12 - 2*x^7 - 1,
         x^12 + x^11 - 2*x^7 - 2*x^6 - x - 1,
         x^12 - x^11 - 2*x^7 + 2*x^6 - x + 1,
         x^12 + x^11 - x^9 - x^8 - x^4 - x^3 - x - 1,
         x^12 - x^11 - x^9 + x^8 - x^4 + x^3 - x + 1,
         x^12 - x^11 - x^10 + 2*x^9 - x^8 + x^6 - 2*x^5 + x^4 - x^2 + x - 1,
         x^12 + x^11 - x^10 - 2*x^9 - x^8 + x^6 - x^4 + x^2 + x + 1,
         x^12 - x^10 - x^9 + x^8 + x^7 - x^6 - x^5 - x^4 + x^3 + x^2 - 1,
         x^12 - x^11 + x^10 - x^9 - x^6 + x^3 - x^2 + x - 1,
         x^12 - x^11 + x^9 - 2*x^7 - x^3 + x + 1,
         x^12 - x^11 - x^9 + 2*x^8 - 2*x^7 + 2*x^4 - x^3 - x + 1,
         x^12 + x^11 - x^9 - 2*x^6 - 2*x^5 + x^3 - x - 1,
         x^12 + x^11 + x^9 + 2*x^8 - 2*x^6 - 2*x^5 - 2*x^4 - 3*x^3 - 4*x^2 - 3*x - 1,
         x^12 + x^11 + x^10 - x^9 - 2*x^8 - 3*x^7 - x^6 + x^5 + 2*x^4 + x^3 - x^2 - x - 1,
         x^12 - 2*x^11 + x^10 - x^9 + x^8 + 2*x^6 - 2*x^5 - x^4 + x^3 - x^2 + 2*x - 1,
         x^12 + x^10 - x^9 - x^8 - 2*x^7 + x^4 + x^3 - x^2 - 1,
         x^12 - x^10 - x^9 - x^8 + 2*x^6 + 2*x^5 - x^4 - x^3 - x^2 + 1,
         x^12 - x^11 + x^10 - x^9 - x^7 + x^6 - x^5 + x^3 - x^2 + x - 1,
         x^12 + 2*x^11 + x^10 - x^9 - 3*x^8 - 4*x^7 - 2*x^6 + 2*x^5 + 3*x^4 + x^3 - x^2 - 2*x - 1,
         x^12 - x^11 + x^9 - x^8 - x^6 + x^4 - x^3 + x - 1,
         x^12 + x^11 - x^10 - x^9 - x^3 - x^2 - x - 1,
         x^12 - x^11 - x^10 + x^9 - x^3 + x^2 - x + 1,
         x^12 - x^11 + x^9 - 2*x^8 + x^6 - 2*x^5 + x^3 + x + 1,
         x^12 - x^11 + x^9 - 2*x^8 + 2*x^7 - x^6 - 2*x^5 + 2*x^4 - x^3 + x - 1,
         x^12 - 2*x^7 - x^6 + 1,
         x^12 - x^11 + x^8 - x^7 - x^6 + x^5 - x^4 + x - 1,
         x^12 - x^11 + x^8 - 2*x^7 + x^6 - x^4 + x - 1,
         x^12 - x^11 + x^10 - x^9 - x^8 + x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^10 - x^9 + x^7 + x^6 - x^5 - 2*x^4 + x^3 + x^2 - 1,
         x^12 - 3*x^11 + 4*x^10 - 4*x^9 + 3*x^8 - 2*x^7 + 2*x^6 - 3*x^4 + 4*x^3 - 4*x^2 + 3*x - 1,
         x^12 - x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^7 - x^4 + 2*x^3 - 2*x^2 + x - 1,
         x^12 + x^11 - x^8 - 2*x^7 - 2*x^6 + x^4 - x - 1,
         x^12 - x^11 - x^8 + 2*x^5 - x^4 - x + 1,
         x^12 - 2*x^11 + 3*x^10 - 3*x^9 + 2*x^8 - 2*x^7 + x^6 - 2*x^4 + 3*x^3 - 3*x^2 + 2*x - 1,
         x^12 + x^10 - x^9 - 2*x^7 - x^6 + x^3 - x^2 - 1,
         x^12 - 3*x^11 + 3*x^10 - x^9 - x^8 + 2*x^7 - 2*x^5 + x^4 - x^3 + 3*x^2 - 3*x + 1,
         x^12 + x^11 + x^10 + x^9 - x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^3 - x^2 - x - 1,
         x^12 + 3*x^11 + 3*x^10 + x^9 - x^8 - 4*x^7 - 6*x^6 - 4*x^5 - x^4 - x^3 - 3*x^2 - 3*x - 1,
         x^12 - x^9 - x^8 - x^7 + x^5 + x^4 - x^3 + 1,
         x^12 + 2*x^11 + 2*x^10 + x^9 - x^8 - 3*x^7 - 4*x^6 - 3*x^5 - x^4 - x^3 - 2*x^2 - 2*x - 1,
         x^12 - x^11 + x^10 - x^9 - x^8 + x^4 - x^3 + x^2 - x + 1,
         x^12 + x^9 - x^8 - x^7 - x^5 - x^4 - x^3 - 1,
         x^12 + x^11 - x^10 - x^9 - x^8 - 2*x^7 + 2*x^5 + x^4 - x^3 - x^2 + x + 1,
         x^12 - 2*x^11 + 2*x^10 - x^9 - x^8 + x^7 - x^5 + x^4 - x^3 + 2*x^2 - 2*x + 1,
         x^12 - x^11 - x^10 + x^9 - x^8 + 2*x^6 - x^4 - x^3 + x^2 + x - 1,
         x^12 - x^10 + x^9 - 2*x^5 - x^3 - x^2 - 1,
         x^12 - 2*x^11 + x^10 + x^9 - 2*x^8 + 2*x^7 - 2*x^6 + 2*x^4 - 3*x^3 + 3*x^2 - 2*x + 1,
         x^12 - x^7 - x^6 - x^5 - 1,
         x^12 + x^11 - x^10 + x^8 - 2*x^7 - 2*x^6 - x^4 + x^2 - x - 1,
         x^12 - 3*x^11 + 3*x^10 - 3*x^8 + 2*x^7 + 2*x^6 - 4*x^5 + 3*x^4 - 3*x^2 + 3*x - 1,
         x^12 - x^11 + x^10 - x^8 - 2*x^5 + x^4 - x^2 + x - 1,
         x^12 - x^11 - x^10 + 2*x^9 - x^8 - 2*x^7 + 2*x^6 - x^4 + 2*x^3 - x^2 - x + 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^8 + x^7 + x^6 - 3*x^5 + 2*x^4 - 2*x^2 + 2*x - 1,
         x^12 - 3*x^11 + 4*x^10 - 4*x^9 + 2*x^8 + 2*x^7 - 4*x^6 + 4*x^5 - 4*x^4 + 4*x^3 - 4*x^2 + 3*x - 1,
         x^12 + x^10 - x^9 - x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 - 1,
         x^12 - 2*x^11 + 3*x^10 - 3*x^9 + x^8 + x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1,
         x^12 + x^11 - 2*x^8 - 2*x^7 - x - 1,
         x^12 - x^11 - 2*x^8 + 2*x^7 - x + 1,
         x^12 - x^11 + 2*x^10 - 2*x^9 - 2*x^6 + 2*x^5 - 2*x^4 + 2*x^3 - 2*x^2 + x - 1,
         x^12 + x^11 - x^10 - x^9 - x^7 + x^5 - x^3 - x^2 - x - 1,
         x^12 - x^11 - x^10 + x^9 - x^7 + 2*x^6 - x^5 - x^3 + x^2 - x + 1,
         x^12 - x^11 - x^6 + 2*x^5 - 2*x^4 + x - 1,
         x^12 + x^11 - x^9 - x^8 - 2*x^7 - x^6 + x^4 + x^3 - x - 1,
         x^12 - x^9 - x^8 - x^6 + x^4 + x^3 - 1,
         x^12 - x^8 - x^7 + x^6 - x^5 - x^4 - 1,
         x^12 - 3*x^11 + 3*x^10 - x^9 + x^7 - 4*x^6 + 5*x^5 - 2*x^4 + x^3 - 3*x^2 + 3*x - 1,
         x^12 - x^11 + x^10 - x^9 + x^7 - 2*x^6 + x^5 - 2*x^4 + x^3 - x^2 + x - 1,
         x^12 + x^11 - x^10 - x^9 + x^7 - 3*x^5 - 2*x^4 + x^3 + x^2 - x - 1,
         x^12 - x^11 - x^10 + x^9 + x^7 - 2*x^6 - x^5 + 2*x^4 + x^3 - x^2 - x + 1,
         x^12 - x^9 + x^7 - x^6 - x^5 - 2*x^4 + x^3 - 1,
         x^12 - 2*x^11 + 2*x^10 - x^9 + x^7 - 3*x^6 + 3*x^5 - 2*x^4 + x^3 - 2*x^2 + 2*x - 1,
         x^12 - x^11 + x^10 - 2*x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^2 + x - 1,
         x^12 - x^11 - x^7 + x^6 + x^5 - 2*x^4 + x - 1,
         x^12 - 2*x^11 + x^9 + x^8 - 2*x^7 + 2*x^6 - x^4 - x^3 + 2*x^2 - 2*x + 1,
         x^12 - 2*x^11 + 3*x^9 - 3*x^8 + 2*x^6 - 2*x^5 + x^4 + x^3 - 2*x^2 + 2*x - 1,
         x^12 - 2*x^10 - x^9 + x^8 + 2*x^5 + x^4 - x^3 - 1,
         x^12 - x^9 - x^8 - 2*x^7 + x^4 + x^3 + 2*x^2 + 1,
         x^12 + 2*x^11 - 3*x^9 - 3*x^8 - 2*x^7 - 2*x^6 + 3*x^4 + 3*x^3 + 2*x^2 + 2*x + 1,
         x^12 + x^9 - x^8 - 2*x^5 - x^4 - x^3 - 2*x^2 - 1,
         x^12 + x^11 - x^6 - 2*x^5 - 2*x^4 - 2*x^3 - 2*x^2 - x - 1,
         x^12 + 2*x^11 - x^9 + x^8 - 2*x^6 - 2*x^5 - 3*x^4 - 3*x^3 - 2*x^2 - 2*x - 1,
         x^12 - x^11 + 2*x^9 - 2*x^8 + x^6 - 2*x^5 - 2*x^2 + x - 1,
         x^12 - 2*x^10 + x^9 + x^8 - 2*x^7 - x^4 + x^3 + 1,
         x^12 - x^11 - 2*x^7 + x^6 + 2*x^2 - x + 1,
         x^12 + x^11 - 2*x^9 - 2*x^8 - 2*x^7 - x^6 + 2*x^4 + 2*x^3 + 2*x^2 + x + 1,
         x^12 - 4*x^11 + 8*x^10 - 11*x^9 + 10*x^8 - 4*x^7 - 4*x^6 + 10*x^5 - 12*x^4 + 11*x^3 - 8*x^2 + 4*x - 1,
         x^12 - 3*x^11 + 3*x^10 - x^9 - 2*x^8 + 6*x^7 - 6*x^6 + 2*x^5 - x^3 + 3*x^2 - 3*x + 1,
         x^12 - x^11 + 3*x^10 - 3*x^9 + 2*x^8 - 2*x^7 - 2*x^6 + 2*x^5 - 4*x^4 + 3*x^3 - 3*x^2 + x - 1,
         x^12 + 3*x^11 + 3*x^10 + x^9 - 2*x^8 - 6*x^7 - 6*x^6 - 2*x^5 - x^3 - 3*x^2 - 3*x - 1,
         x^12 - x^11 - x^8 + x^7 + x^6 - x^5 - x^4 + x - 1,
         x^12 - 2*x^11 + 2*x^10 - x^9 - 2*x^8 + 4*x^7 - 4*x^6 + 2*x^5 - x^3 + 2*x^2 - 2*x + 1,
         x^12 - x^11 + 2*x^10 - 2*x^9 + x^8 - x^7 - x^6 + x^5 - 3*x^4 + 2*x^3 - 2*x^2 + x - 1,
         x^12 - x^11 - x^10 + x^9 - 2*x^8 + 2*x^7 + 2*x^6 - 2*x^5 - x^3 + x^2 + x - 1,
         x^12 + x^11 - x^10 - x^9 - 2*x^8 - 2*x^7 + 2*x^6 + 2*x^5 - x^3 - x^2 + x + 1,
         x^12 + x^10 - x^8 - x^7 - x^6 - x^5 - x^4 - x^2 - 1,
         x^12 + 2*x^11 + 2*x^10 + x^9 - 2*x^8 - 4*x^7 - 4*x^6 - 2*x^5 - x^3 - 2*x^2 - 2*x - 1,
         x^12 - 2*x^11 + 4*x^10 - 5*x^9 + 4*x^8 - 2*x^7 - 2*x^6 + 4*x^5 - 6*x^4 + 5*x^3 - 4*x^2 + 2*x - 1,
         x^12 - x^9 - 2*x^8 + 2*x^5 - x^3 + 1,
         x^12 + x^11 + x^10 + x^9 - 2*x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^3 - x^2 - x - 1,
         x^12 + 2*x^10 - x^9 - 2*x^7 - 2*x^6 - 2*x^4 + x^3 - 2*x^2 - 1,
         x^12 - 3*x^11 + 6*x^10 - 8*x^9 + 7*x^8 - 3*x^7 - 3*x^6 + 7*x^5 - 9*x^4 + 8*x^3 - 6*x^2 + 3*x - 1,
         x^12 - x^11 + x^10 - x^9 - 2*x^8 + 2*x^7 - 2*x^6 + 2*x^5 - x^3 + x^2 - x + 1,
         x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 3*x^8 - x^7 - x^6 + 3*x^5 - 5*x^4 + 4*x^3 - 3*x^2 + 2*x - 1,
         x^12 - x^11 + x^10 - x^9 - 2*x^4 + x^3 - x^2 + x - 1,
         x^12 + x^9 - 2*x^8 - 2*x^5 - x^3 - 1,
         x^12 - 2*x^11 + 2*x^10 - 3*x^9 + 2*x^8 + 2*x^5 - 4*x^4 + 3*x^3 - 2*x^2 + 2*x - 1,
         x^12 + x^11 + 2*x^10 - x^8 - 3*x^7 - 3*x^6 - x^5 - x^4 - 2*x^2 - x - 1,
         x^12 - 3*x^11 + 5*x^10 - 7*x^9 + 6*x^8 - 2*x^7 - 2*x^6 + 6*x^5 - 8*x^4 + 7*x^3 - 5*x^2 + 3*x - 1,
         x^12 - 5*x^11 + 11*x^10 - 15*x^9 + 14*x^8 - 6*x^7 - 6*x^6 + 14*x^5 - 16*x^4 + 15*x^3 - 11*x^2 + 5*x - 1,
         x^12 - x^9 - x^7 - 2*x^6 - x^5 + x^3 + 2*x^2 + 2*x + 1,
         x^12 + x^11 - 2*x^10 - 3*x^9 + x^8 + 3*x^7 - x^5 - x^4 - x^3 + x + 1,
         x^12 + x^11 - x^9 - x^8 - x^7 + x^5 - x^4 - x^3 - 2*x^2 - x - 1,
         x^12 - x^11 - 2*x^10 + x^9 + 3*x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 + x - 1,
         x^12 + 3*x^11 + 2*x^10 - 3*x^9 - 5*x^8 - x^7 + 2*x^6 + x^5 - x^4 - 3*x^3 - 4*x^2 - 3*x - 1,
         x^12 - x^10 + x^8 - x^7 - x^6 + x^5 - x^4 - x^2 - 1,
         x^12 + 2*x^11 + x^10 - 2*x^9 - 3*x^8 - x^7 + x^6 + x^5 - x^4 - 2*x^3 - 3*x^2 - 2*x - 1,
         x^12 - x^10 + x^9 - 2*x^7 + x^6 - 2*x^4 - x^3 + x^2 - 1,
         x^12 + x^11 + x^8 - 2*x^6 - 2*x^5 - 3*x^4 - 2*x^3 - 2*x^2 - x - 1,
         x^12 + x^11 - x^9 - 2*x^8 - x^7 + x^5 - x^3 - x - 1,
         x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 4*x^8 - 4*x^7 + 3*x^6 - 2*x^5 + 2*x^3 - 3*x^2 + 2*x - 1,
         x^12 + x^10 - 2*x^9 - 2*x^7 + x^6 - x^2 - 1,
         x^12 - x^11 + 2*x^10 - 3*x^9 + 2*x^8 - 3*x^7 + 2*x^6 - x^5 + x^3 - 2*x^2 + x - 1,
         x^12 - x^11 - x^9 + x^7 + x^5 - 2*x^4 + x^3 - x + 1,
         x^12 - 3*x^11 + 4*x^10 - 5*x^9 + 6*x^8 - 5*x^7 + 4*x^6 - 3*x^5 + 3*x^3 - 4*x^2 + 3*x - 1,
         x^12 - x^9 - 2*x^7 + x^6 - x^3 + 1,
         x^12 - x^8 - 2*x^7 + x^4 - 1,
         x^12 + x^11 + x^10 - x^9 - 2*x^8 - 3*x^7 - x^6 - x^5 + x^3 + x^2 + x + 1,
         x^12 - x^10 - x^7 - x^6 + x^5 + x^2 - 1,
         x^12 - x^8 - x^7 - x^6 - x^5 + x^4 + 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 3*x^8 - 3*x^7 + x^6 - x^5 + x^4 - 2*x^2 + 2*x - 1,
         x^12 - x^11 - x^10 + x^9 + x^8 - 2*x^7 + x^4 + x^3 - x^2 - x + 1,
         x^12 + x^11 - x^10 - x^9 + x^8 - 2*x^6 - 2*x^5 - x^4 + x^3 + x^2 - x - 1,
         x^12 - x^11 - 2*x^10 + 2*x^9 + x^8 - x^7 - x^5 + x^4 + x - 1,
         x^12 - x^10 - x^9 + x^7 - x^5 - x^3 + x^2 + 1,
         x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - 2*x^3 + 2*x^2 - x + 1,
         x^12 + x^11 - x^8 - x^7 - x^5 - x^4 - 2*x^3 - 2*x^2 - x - 1,
         x^12 + 2*x^11 + x^10 - x^9 - 2*x^8 - x^7 - x^5 - 2*x^4 - 3*x^3 - 3*x^2 - 2*x - 1,
         x^12 - 2*x^11 + x^10 + x^9 - 2*x^8 + x^7 - x^5 + 2*x^4 - 3*x^3 + 3*x^2 - 2*x + 1,
         x^12 - x^10 + x^9 - x^7 - x^5 - x^3 - x^2 - 1,
         x^12 + 3*x^11 + 2*x^10 - 2*x^9 - 3*x^8 - x^7 - x^5 - 3*x^4 - 4*x^3 - 4*x^2 - 3*x - 1,
         x^12 + x^11 - 2*x^10 - 2*x^9 + x^8 + x^7 - x^5 - x^4 + x + 1,
         x^12 - 3*x^11 + 2*x^10 + 2*x^9 - 3*x^8 + x^7 - x^5 + 3*x^4 - 4*x^3 + 4*x^2 - 3*x + 1,
         x^12 - x^9 - x^7 - x^5 + x^3 - 1,
         x^12 - 2*x^7 - x^6 - 1,
         x^12 - x^11 - x^10 + x^9 + x^8 - 3*x^6 + 2*x^5 + x^4 - 3*x^3 + x^2 + x - 1,
         x^12 - 2*x^11 + x^10 - x^2 + 1,
         x^12 - x^10 - x^2 - 2*x - 1,
         x^12 + x^11 - x^7 - x^6 - x^5 - 2*x^4 - 2*x^3 - 2*x^2 - x - 1,
         x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 4*x^8 - 5*x^7 + 5*x^6 - 3*x^5 + 2*x^3 - 3*x^2 + 2*x - 1,
         x^12 - x^8 - x^6 - 2*x^5 + x^4 - 1,
         x^12 - x^11 + x^10 - x^9 + 2*x^8 - 2*x^7 - x^6 - x^3 - x^2 + x - 1,
         x^12 - 2*x^11 + x^10 + x^8 - 2*x^7 + 2*x^6 - 4*x^5 + 5*x^4 - 4*x^3 + 3*x^2 - 2*x + 1,
         x^12 - x^10 + x^8 - 2*x^5 - x^4 - x^2 - 1,
         x^12 + x^11 - x^9 - x^8 - x^6 - 2*x^5 - x^4 + x^3 - x - 1,
         x^12 - x^10 + x^9 + x^8 - x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - 1,
         x^12 - x^11 + x^10 - 2*x^9 + x^8 - 2*x^7 + 2*x^6 - 2*x^5 + 3*x^4 - 2*x^3 + x^2 - x + 1,
         x^12 + x^11 + x^10 - x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^2 - x - 1,
         x^12 - 3*x^11 + 5*x^10 - 6*x^9 + 5*x^8 - 4*x^7 + 4*x^6 - 4*x^5 + 3*x^4 - 2*x^3 + x^2 - x + 1,
         x^12 - x^11 + x^10 - x^8 - x^4 - x^2 - x - 1,
         x^12 - x^10 - x^7 - x^5 + x^2 + 1,
         x^12 + 2*x^11 + x^10 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 4*x^3 - 3*x^2 - 2*x - 1,
         x^12 - x^10 - x^9 + 2*x^7 - x^6 - 2*x^5 + x^3 + x^2 - 1,
         x^12 - x^11 + x^8 - x^7 - x^6 - x^5 + x^4 + x - 1,
         x^12 - x^11 + x^9 - x^8 - x^7 + x^6 + x^5 - x^4 - x^3 - x - 1,
         x^12 + x^11 - x^10 - x^9 - 2*x^5 - 2*x^4 + x^3 + x^2 - x - 1,
         x^12 - x^11 - x^10 + x^9 - 2*x^5 + 2*x^4 + x^3 - x^2 - x + 1,
         x^12 + x^11 - 3*x^7 - 3*x^6 - x^5 + x + 1,
         x^12 - x^11 + x^10 - x^9 - x^7 - x^6 + x^5 + x^3 - x^2 + x - 1,
         x^12 - x^10 + x^9 + x^8 - 3*x^7 - x^6 + x^5 - x^4 - x^3 + x^2 - 1,
         x^12 - x^11 - x^8 + 2*x^6 - x^4 - x - 1,
         x^12 + x^9 - x^7 - x^6 - x^5 - 2*x^4 - x^3 - 2*x^2 - 1,
         x^12 + x^11 - x^10 - x^9 - x^7 - 2*x^6 - x^5 + x^3 + x^2 + x + 1,
         x^12 + 2*x^11 + 2*x^10 + x^9 - x^7 - 3*x^6 - 5*x^5 - 6*x^4 - 5*x^3 - 4*x^2 - 2*x - 1,
         x^12 + 3*x^11 + 3*x^10 + x^9 - x^7 - 4*x^6 - 7*x^5 - 8*x^4 - 7*x^3 - 5*x^2 - 3*x - 1,
         x^12 - x^11 - x^10 + x^9 - x^7 + x^5 + x^3 - x^2 + x - 1,
         x^12 + x^11 + x^10 + x^9 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 3*x^3 - 3*x^2 - x - 1,
         x^12 + x^11 - x^9 - x^8 - x^7 - 2*x^6 - x^5 + x^4 + x^3 - x - 1,
         x^12 - x^11 - x^9 + x^8 - x^7 + x^5 + x^4 - x^3 - x + 1,
         x^12 - x^10 - x^6 - 2*x^5 + 2*x^3 + x^2 - 1,
         x^12 - x^8 - x^7 - x^6 + x^5 - x^4 - 1,
         x^12 + 2*x^11 + x^10 - 2*x^9 - 3*x^8 - x^7 - x^5 - x^4 - x^2 - 2*x - 1,
         x^12 - x^10 - 2*x^9 + x^8 + x^7 - x^5 + x^4 - x^2 + 1,
         x^12 - x^10 + x^8 - 2*x^7 - x^6 + x^4 - x^2 + 1,
         x^12 - x^9 - x^8 - x^4 + x^3 - 1,
         x^12 - 2*x^8 - x^7 + x^6 - x^5 + 1,
         x^12 - x^11 + x^10 - x^9 + x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^11 + x^10 - x^9 - 2*x^7 + x^6 + x^3 - x^2 + x - 1,
         x^12 + x^10 - 2*x^9 - x^8 - 2*x^7 - x^6 + 2*x^5 + x^4 + 2*x^3 - x^2 - 1,
         x^12 + x^11 - 2*x^9 - 2*x^8 - x^7 + x^6 + x^5 - x - 1,
         x^12 - x^8 - 2*x^5 - x^4 - 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^5 + 3*x^4 - 2*x^3 + 2*x^2 - 2*x + 1,
         x^12 + 2*x^11 + x^10 - 2*x^8 - 5*x^7 - 4*x^6 - x^5 + x^2 + 2*x + 1,
         x^12 - x^11 + x^10 - 2*x^8 + x^7 - x^6 - x^5 + x^2 - x + 1,
         x^12 - x^10 - 2*x^8 - x^7 + 2*x^6 + x^5 + x^2 - 1,
         x^12 - 2*x^11 + x^10 - 2*x^8 + 3*x^7 - x^5 + x^2 - 2*x + 1,
         x^12 + x^11 + x^10 - 2*x^8 - 3*x^7 - 3*x^6 - x^5 + x^2 + x + 1,
         x^12 + x^10 - 2*x^8 - x^7 - 2*x^6 - x^5 + x^2 + 1,
         x^12 + x^11 - 2*x^9 - x^8 - 2*x^5 - x^4 + x + 1,
         x^12 + x^9 - 2*x^8 - x^7 - 3*x^5 + x^3 + 1,
         x^12 - x^11 + x^10 - 2*x^9 + x^8 - x^6 - x^4 + 2*x^3 - x^2 + x - 1,
         x^12 - 2*x^11 + 2*x^10 - x^9 - x^8 + 3*x^7 - 2*x^6 - x^5 + x^4 - x^3 - 1,
         x^12 - x^11 + x^10 - x^9 + x^8 - 2*x^7 - x^4 + x^3 - x^2 + x - 1,
         x^12 + x^10 - x^9 - x^8 - 2*x^7 - x^6 - x^4 + x^3 + x^2 + 1,
         x^12 + x^11 - 2*x^7 - 2*x^6 - 2*x^4 - 2*x^3 - x - 1,
         x^12 - x^11 - 2*x^7 + 2*x^6 - 2*x^4 + 2*x^3 - x + 1,
         x^12 + x^11 - 3*x^9 - 2*x^8 + 2*x^6 + x^3 - x - 1,
         x^12 - x^11 - x^9 + x^8 - x^6 + x^4 + x^3 - x - 1,
         x^12 + x^11 - x^9 - 2*x^8 - x^7 - x^5 + x^3 - x - 1,
         x^12 - x^11 - x^9 + x^7 - x^5 + 2*x^4 - x^3 - x + 1,
         x^12 - 2*x^10 - x^9 + 2*x^8 - x^6 + 2*x^5 - x^3 - 1,
         x^12 - 2*x^11 + x^10 + x^9 - x^8 - x^6 + 2*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1,
         x^12 - x^11 - x^10 + x^9 + x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^11 - x^10 + 2*x^9 - x^8 - x^7 + x^6 - x^5 + x^4 - x^2 + x - 1,
         x^12 - x^10 - x^8 - x^7 + x^6 + x^5 - x^4 + x^2 - 1,
         x^12 - x^11 - 2*x^8 + 3*x^7 - x^6 - x^5 + 2*x^4 - 2*x^3 + x - 1,
         x^12 - x^10 + x^9 - 3*x^7 + x^6 + x^5 - 2*x^4 - x^3 + x^2 - 1,
         x^12 - x^11 + x^7 - 3*x^6 + x^5 + x - 1,
         x^12 - x^11 + x^9 - 2*x^8 + x^7 - x^6 - x^5 + 2*x^4 - x^3 + x - 1,
         x^12 - x^11 + x^10 - x^9 - x^8 + x^7 - 3*x^6 + 3*x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^10 - x^7 + x^5 - 2*x^4 + x^2 - 1,
         x^12 - x^11 - x^10 + 2*x^9 - x^7 - x^5 - x^2 + x - 1,
         x^12 - x^11 + x^10 - x^9 - x^6 + x^3 - x^2 - x - 1,
         x^12 - x^11 + x^9 - x^7 - x^6 + x^5 - 2*x^4 + x^3 - 2*x^2 + x - 1,
         x^12 - x^11 + x^10 - x^8 + x^7 - 3*x^6 + x^5 - x^4 - x^2 + x - 1,
         x^12 - x^8 - 2*x^6 + x^4 - 2*x^3 - 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^6 + 4*x^5 - 3*x^4 + 2*x^2 - 2*x + 1,
         x^12 + x^10 - 2*x^9 + x^8 - 2*x^7 + x^6 - 2*x^5 - x^4 - x^2 - 1,
         x^12 - x^11 + x^9 - 2*x^7 - x^3 + x - 1,
         x^12 - x^11 + x^8 - 2*x^7 - x^6 + 2*x^5 - x^4 + x - 1,
         x^12 + x^11 + 2*x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - 3*x^4 - 2*x^3 - 2*x^2 - x - 1,
         x^12 + x^11 - x^9 - x^8 - x^7 + x^6 - x^5 - 3*x^4 - x^3 - x - 1,
         x^12 + x^11 + x^10 - x^9 - 2*x^8 - 2*x^7 - x^6 - x^3 - x^2 - x - 1,
         x^12 - x^9 - x^7 - x^6 + x^5 - x^3 - 1,
         x^12 - 2*x^11 + x^10 - x^8 + x^7 + x^5 - x^4 - x^2 + 2*x - 1,
         x^12 - x^10 - x^8 - x^7 + x^5 + x^4 - x^2 + 1,
         x^12 + x^10 - x^8 - x^7 - 2*x^6 - x^5 - x^4 - x^2 - 1,
         x^12 - x^11 + x^10 - x^8 - x^6 - x^4 - x^2 + x - 1,
         x^12 + 2*x^11 + x^10 - x^8 - 3*x^7 - 4*x^6 - 3*x^5 - x^4 - x^2 - 2*x - 1,
         x^12 + x^11 + x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - x^4 - x^2 - x - 1,
         x^12 - x^8 - x^7 + x^6 + x^5 - x^4 - 2*x^3 - 2*x^2 - 2*x - 1,
         x^12 - x^9 - x^8 - x^7 + x^5 + x^4 - x^3 - 1,
         x^12 - 2*x^11 + x^10 - x^9 + x^8 + x^7 - 2*x^6 + 3*x^5 - 3*x^4 + x^3 - x^2 + 2*x - 1,
         x^12 + x^11 + x^10 - x^9 - 2*x^8 - 2*x^7 - 2*x^6 + x^3 - x^2 - x - 1,
         x^12 + x^10 - x^9 - x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 - x^2 - 1,
         x^12 - x^10 - x^9 - x^8 + x^7 + x^5 + x^4 - x^3 - x^2 + 1,
         x^12 - x^11 + x^10 - x^9 - 2*x^6 + 2*x^5 - 2*x^4 + x^3 - x^2 + x - 1,
         x^12 + 2*x^11 + x^10 - x^9 - 3*x^8 - 3*x^7 - 2*x^6 - x^5 + x^4 + x^3 - x^2 - 2*x - 1,
         x^12 - x^7 - 2*x^6 - x^5 - 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 3*x^7 + 2*x^6 - x^5 + 2*x^4 - 2*x^3 + 2*x^2 - 2*x + 1,
         x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 + x^7 - 2*x^6 + x^5 - x^4 + 1,
         x^12 - x^8 + x^7 - x^5 - x^4 - 2*x^3 - 2*x^2 - 2*x - 1,
         x^12 - x^8 - x^7 - x^6 - x^5 + x^4 - 1,
         x^12 - x^10 + x^8 - x^7 - x^5 - x^4 - x^2 - 1,
         x^12 - 2*x^11 + x^10 + x^8 - 3*x^7 + 4*x^6 - 5*x^5 + 5*x^4 - 4*x^3 + 3*x^2 - 2*x + 1,
         x^12 - x^11 - x^10 + 2*x^9 - 2*x^7 + 2*x^6 - 2*x^5 - x^2 + x - 1,
         x^12 + x^11 - x^9 - x^8 - 2*x^7 - 2*x^6 + x^4 + x^3 - x - 1,
         x^12 - x^11 - x^9 + x^8 - 2*x^7 + 2*x^6 + x^4 - x^3 - x + 1,
         x^12 + x^11 - x^10 - 2*x^9 + x^7 - x^5 - 2*x^4 + x^2 - x - 1,
         x^12 - x^11 - x^10 + 2*x^8 - x^7 - x^5 + 2*x^3 - x^2 - x + 1,
         x^12 + x^11 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 + x^4 + 2*x^3 - x - 1,
         x^12 - x^11 + x^10 - x^9 - x^8 - x^6 + 2*x^5 - x^4 + x^3 - x^2 + x - 1,
         x^12 - x^11 + x^10 - x^9 - x^6 - 2*x^5 + 2*x^4 - x^3 + x^2 - x + 1,
         x^12 + x^11 - 2*x^9 - x^8 - x^7 - x^6 - x^5 + x^4 + x + 1,
         x^12 - 2*x^7 - 2*x^6 + 1,
         x^12 - x^10 - x^7 + x^6 - x^5 - x^2 - 1,
         x^12 - x^9 - x^7 - 2*x^6 + x^5 + x^3 - 1,
         x^12 - x^11 - x^6 + x - 1],

    14: [x^14 + x^13 - x^8 - 2*x^7 - x^6 - x - 1,
         x^14 - x^13 - x^8 + x^6 - x + 1,
         x^14 - x^13 - x^9 + x^8 - x^6 + x^5 - x + 1,
         x^14 + x^13 - x^9 - x^8 - x^6 - x^5 - x - 1,
         x^14 - x^13 - x^10 + x^9 - x^5 + x^4 - x + 1,
         x^14 + x^13 - x^10 - x^9 - x^5 - x^4 - x - 1,
         x^14 - 3*x^13 + 4*x^12 - 4*x^11 + 4*x^10 - 4*x^9 + 2*x^8 + 2*x^7 - 4*x^6 + 4*x^5 - 4*x^4 + 4*x^3 - 4*x^2 + 3*x - 1,
         x^14 - 2*x^13 + 3*x^12 - 3*x^11 + 3*x^10 - 3*x^9 + x^8 + x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1,
         x^14 + x^13 - 2*x^8 - 2*x^7 - x - 1,
         x^14 - x^13 - 2*x^8 + 2*x^7 - x + 1,
         x^14 - x^13 + 2*x^12 - 2*x^11 + 2*x^10 - 2*x^9 - 2*x^6 + 2*x^5 - 2*x^4 + 2*x^3 - 2*x^2 + x - 1,
         x^14 + x^12 - x^11 + x^10 - x^9 - x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 - 1,
         x^14 - x^11 - x^10 - x^7 + x^4 + x^3 + 1,
         x^14 - x^11 - x^10 + x^7 - x^4 + x^3 - 1,
         x^14 - x^13 - x^11 + x^10 - x^4 + x^3 - x + 1,
         x^14 + x^13 - x^11 - x^10 - x^4 - x^3 - x - 1,
         x^14 - x^10 - x^9 - x^7 + x^5 + x^4 - 1,
         x^14 - x^9 - x^8 - x^6 + x^5 - 1,
         x^14 - 2*x^8 - x^7 + 1,
         x^14 - 3*x^13 + 3*x^12 - x^11 - x^9 + 2*x^8 - 2*x^6 + x^5 - x^3 + 3*x^2 - 3*x + 1,
         x^14 + 2*x^13 + 2*x^12 + x^11 - x^9 - 3*x^8 - 4*x^7 - 3*x^6 - x^5 - x^3 - 2*x^2 - 2*x - 1,
         x^14 - 2*x^13 + 2*x^12 - x^11 - x^9 + x^8 - x^6 + x^5 - x^3 + 2*x^2 - 2*x + 1,
         x^14 - x^11 - x^9 - x^8 + x^6 + x^5 - x^3 + 1,
         x^14 + 3*x^13 + 3*x^12 + x^11 - x^9 - 4*x^8 - 6*x^7 - 4*x^6 - x^5 - x^3 - 3*x^2 - 3*x - 1,
         x^14 + x^11 - x^9 - x^8 - x^6 - x^5 - x^3 - 1,
         x^14 + x^13 - x^12 - x^11 - x^9 - 2*x^8 + 2*x^6 + x^5 - x^3 - x^2 + x + 1,
         x^14 - x^13 - x^12 + x^11 - x^9 + 2*x^7 - x^5 - x^3 + x^2 + x - 1,
         x^14 - x^13 + x^12 - x^11 - x^9 + x^5 - x^3 + x^2 - x + 1,
         x^14 + x^13 + x^12 + x^11 - x^9 - 2*x^8 - 2*x^7 - 2*x^6 - x^5 - x^3 - x^2 - x - 1,
         x^14 - 2*x^13 + 2*x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^5 - 2*x^4 + 2*x^2 - 2*x + 1,
         x^14 - 2*x^6 - 2*x^3 - 1,
         x^14 - 2*x^9 - 1,
         x^14 + x^13 - x^9 - 2*x^8 - 2*x^7 + x^5 - x - 1,
         x^14 - x^13 - x^9 + 2*x^6 - x^5 - x + 1,
         x^14 - x^9 - x^8 - x^7 + x^6 - x^5 + 1,
         x^14 + x^13 - x^11 - x^10 - x^7 - 2*x^6 + x^4 + x^3 - x - 1,
         x^14 - x^13 - 2*x^9 + 2*x^8 - x + 1,
         x^14 + x^13 - 2*x^9 - 2*x^8 - x - 1,
         x^14 - 3*x^13 + 3*x^12 - x^11 - x^10 + 3*x^9 - 3*x^8 + 3*x^6 - 3*x^5 + x^4 - x^3 + 3*x^2 - 3*x + 1,
         x^14 - x^13 + x^12 - x^11 - x^10 + x^9 - x^8 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1,
         x^14 - x^13 - x^12 + x^11 - x^10 + x^9 + x^8 - 2*x^7 + x^6 + x^5 - x^4 - x^3 + x^2 + x - 1,
         x^14 + 3*x^13 + 3*x^12 + x^11 - x^10 - 3*x^9 - 3*x^8 - 2*x^7 - 3*x^6 - 3*x^5 - x^4 - x^3 - 3*x^2 - 3*x - 1,
         x^14 + x^13 - x^12 - x^11 - x^10 - x^9 + x^8 - x^6 + x^5 + x^4 - x^3 - x^2 + x + 1,
         x^14 - 2*x^13 + 2*x^12 - x^11 - x^10 + 2*x^9 - 2*x^8 + 2*x^6 - 2*x^5 + x^4 - x^3 + 2*x^2 - 2*x + 1,
         x^14 + x^11 - x^10 - 2*x^7 - x^4 - x^3 - 1,
         x^14 + 2*x^13 + 2*x^12 + x^11 - x^10 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^3 - 2*x^2 - 2*x - 1,
         x^14 + x^13 + x^12 + x^11 - x^10 - x^9 - x^8 - 2*x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - x - 1,
         x^14 - x^11 - x^10 + x^4 - x^3 + 1,
         x^14 - 2*x^13 + x^12 - x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^6 - 2*x^4 + x^3 + x^2 - 2*x + 1,
         x^14 - x^12 - x^11 + 2*x^5 - x^3 + x^2 - 1,
         x^14 + 3*x^13 + 5*x^12 + 5*x^11 + 3*x^10 - 3*x^8 - 6*x^7 - 9*x^6 - 10*x^5 - 9*x^4 - 7*x^3 - 5*x^2 - 3*x - 1,
         x^14 + 4*x^13 + 7*x^12 + 7*x^11 + 4*x^10 - 4*x^8 - 8*x^7 - 12*x^6 - 14*x^5 - 12*x^4 - 9*x^3 - 7*x^2 - 4*x - 1,
         x^14 + 2*x^13 + 3*x^12 + 3*x^11 + 2*x^10 - 2*x^8 - 4*x^7 - 6*x^6 - 6*x^5 - 6*x^4 - 5*x^3 - 3*x^2 - 2*x - 1,
         x^14 + x^13 + x^12 + x^11 + x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - 3*x^4 - 3*x^3 - x^2 - x - 1,
         x^14 - x^13 + x^12 - x^11 + x^10 - 2*x^9 + x^8 - 2*x^7 + x^6 - x^4 + x^3 + x^2 - x + 1,
         x^14 + x^13 + x^12 - x^11 - x^10 - 2*x^9 - x^8 - 2*x^7 - x^6 + x^4 + x^3 + x^2 + x + 1,
         x^14 + x^12 - x^11 - 2*x^9 - 2*x^7 + x^3 + x^2 + 1,
         x^14 + 2*x^13 + x^12 - x^11 - 2*x^10 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 + 2*x^4 + x^3 + x^2 + 2*x + 1,
         x^14 - x^13 - x^12 + x^11 - x^3 + x^2 - x + 1,
         x^14 + x^13 - x^12 - x^11 - x^3 - x^2 - x - 1,
         x^14 + 2*x^13 + x^12 - 2*x^9 - 4*x^8 - 2*x^7 - x^2 - 2*x - 1,
         x^14 - 2*x^13 + x^12 - 2*x^9 + 4*x^8 - 2*x^7 - x^2 + 2*x - 1,
         x^14 + x^13 + x^12 - 2*x^9 - 2*x^8 - 2*x^7 - x^2 - x - 1,
         x^14 + x^12 - 2*x^9 - 2*x^7 - x^2 - 1,
         x^14 - x^13 + x^12 - 2*x^9 + 2*x^8 - 2*x^7 - x^2 + x - 1,
         x^14 - x^12 - 2*x^9 + 2*x^7 - x^2 + 1,
         x^14 - x^13 + x^12 - x^11 - x^7 + x^3 - x^2 + x - 1,
         x^14 - x^10 - x^9 + x^5 - x^4 - 1,
         x^14 - x^12 - x^11 + x^9 + x^8 - 2*x^7 - x^6 + x^5 + 2*x^4 - x^3 - x^2 + 1,
         x^14 + 2*x^13 + x^12 - x^11 - 2*x^10 - x^9 + x^8 - 3*x^6 - 3*x^5 + x^3 - x^2 - 2*x - 1,
         x^14 + x^13 + x^12 - 2*x^10 - 2*x^9 - 2*x^8 - x^7 + x^2 + x + 1,
         x^14 - x^13 + x^12 - x^11 - x^8 + x^7 - x^6 + x^3 - x^2 + x - 1,
         x^14 - 3*x^13 + 3*x^12 - x^11 - 2*x^9 + 6*x^8 - 6*x^7 + 2*x^6 - x^3 + 3*x^2 - 3*x + 1,
         x^14 + x^13 - x^12 - x^11 - 2*x^9 - 2*x^8 + 2*x^7 + 2*x^6 - x^3 - x^2 + x + 1,
         x^14 - x^13 + x^12 - x^11 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^6 - x^3 + x^2 - x + 1,
         x^14 - 2*x^13 + 2*x^12 - x^11 - 2*x^9 + 4*x^8 - 4*x^7 + 2*x^6 - x^3 + 2*x^2 - 2*x + 1,
         x^14 + x^13 + x^12 + x^11 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 - x^3 - x^2 - x - 1,
         x^14 - x^11 - 2*x^9 + 2*x^6 - x^3 + 1,
         x^14 + x^11 - 2*x^9 - 2*x^6 - x^3 - 1,
         x^14 + 2*x^13 + 2*x^12 + x^11 - 2*x^9 - 4*x^8 - 4*x^7 - 2*x^6 - x^3 - 2*x^2 - 2*x - 1,
         x^14 + 3*x^13 + 3*x^12 + x^11 - 2*x^9 - 6*x^8 - 6*x^7 - 2*x^6 - x^3 - 3*x^2 - 3*x - 1,
         x^14 - x^13 - x^12 + x^11 - 2*x^9 + 2*x^8 + 2*x^7 - 2*x^6 - x^3 + x^2 + x - 1,
         x^14 - x^8 - x^7 - x^6 - 1]
}
# deg_8_short = [
#     x^8 - x^7 - x^5 + x^3 - x + 1,
#     x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1,
#     x^8 - x^5 - x^4 - x^3 - 1
# ]



def partitions_compact(n):
    temp = Partitions(n).list()
    final = []
    for partition in temp:
        s = set(partition)
        final.append([(k,list(partition).count(k)) for k in s])
    return final

def strata(genus):
    """
    Return the possible strata for the genus g orientable surface that lift
    from a nonorientable surface

    sage: strata(6)

    [[(12, 2)],
     [(4, 2), (10, 2)],
     [(6, 2), (8, 2)],
     [(4, 4), (8, 2)],
     [(4, 2), (6, 4)],
     [(4, 6), (6, 2)],
     [(4, 10)]]
    """
    return [[(2*item[0]+2,2*item[1]) for item in partition] for partition in partitions_compact(genus-1)]
    temp = Partitions(genus-1).list()

# def new_orbit_patterns(set_size):
#     """
#     Return a list whose elements describe the possible orbits of
#     ``set_size`` singularities of the same degree. Every element of the
#     list is also a list and it describes one scenario. Elements of each list
#     are tuples (``orbit_length``, ``num_orbits``). Every orbit length has to be
#     even unless it is one. The number of orbits of length 1 has to be even, so
#     in that case we return (1, ``num_orbits``/2).
#     """
#     assert(set_size % 2) == 0
#     half_mult = set_size // 2
#     temp = [[(2*k,mult) for (k,mult) in partition]
#             for partition in partitions_compact(half_mult)]
#     # return temp
#     final = []
#     for partition in temp:
#         final.append(partition)
#         if partition[0][0] != 2:
#             continue
#         for i in range(partition[0][1]):
#             start = [(1,1*i+1)]  # the true multiplicity is double, but this is
#             # convenient, since these come in pairs
#             if partition[0][1]-i-1 > 0:
#                 start.append((2,partition[0][1]-i-1))
#             final.append(start + partition[1:])
#     return final


def is_symmetric(partition):
    return all(orbit%2 == 0 or mult%2 ==0 for (orbit, mult) in partition)

def orbit_patterns(set_size):
    """
    sage: orbit_patterns(4)
    [[(4, 1)], [(2, 2)], [(1, 2), (2, 1)], [(1, 4)]]
    sage: orbit_patterns(6)
    [[(6, 1)],
    [(2, 1), (4, 1)],
    [(1, 2), (4, 1)],
    [(3, 2)],
    [(2, 3)],
    [(1, 2), (2, 2)],
    [(1, 4), (2, 1)],
    [(1, 6)]]
    """
    assert(set_size % 2) == 0
    return filter(is_symmetric, partitions_compact(set_size))

def append_options(list_with_mult, options):
    """
    INPUT:
    - ``list_with_mult`` -- list of tuples (data, multiplicity)
    - ``options`` -- list of options to append on the pieces of data

    OUTPUT:
    A list of possible appendings. Each appending is a list, containing tuples
    (data, multiplicity, option), counting how many pieces of data receives a
    specific option.

    EXAMPLE:

    sage: append_options([(10,1), (13,3)],[1,5])
    [[(10, 1, 1), (13, 3, 1)],
     [(10, 1, 1), (13, 2, 1), (13, 1, 5)],
     [(10, 1, 1), (13, 1, 1), (13, 2, 5)],
     [(10, 1, 1), (13, 3, 5)],
     [(10, 1, 5), (13, 3, 1)],
     [(10, 1, 5), (13, 2, 1), (13, 1, 5)],
     [(10, 1, 5), (13, 1, 1), (13, 2, 5)],
     [(10, 1, 5), (13, 3, 5)]]
    """
    result = []
    append_options_recursive(0, list_with_mult, options, [], result)
    return result

def append_options_recursive(current_idx, list_with_mult, options, list_so_far, result):
    if current_idx == len(list_with_mult):
        result.append(list_so_far)
        return
    num_options = len(options)
    data, multiplicity = list_with_mult[current_idx]
    if data % 2 == 1:
        assert multiplicity % 2 == 0
        multiplicity /= 2
    comps = Compositions(num_options+multiplicity, length=num_options).list()
    comps = [[i-1 for i in comp] for comp in comps]
    for comp in comps:
        new_list = list(list_so_far)
        for k in range(num_options):
            mult = comp[k] if data % 2 == 0 else 2*comp[k]
            if mult > 0:
                new_list.append((data, mult, options[k]))
        append_options_recursive(current_idx + 1, list_with_mult, options,
                                 new_list, result)


Orbit = namedtuple("Orbit", "num_prongs orbit_length num_orbits order_after_first_return")

def orbit_combinations_fixed_num_prongs(num_prongs, multiplicity):
    """

    OUTPUT:

    list of list of Orbits.

    EXAMPLE:

    sage:     sage: orbit_combinations_fixed_num_prongs(6, 2)

    [[Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=1)],
     [Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3)],
     [Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=1)],
     [Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=3)]]

    sage: orbit_combinations_fixed_num_prongs(2, 4)
    [[Orbit(num_prongs=2, orbit_length=4, num_orbits=1, order_after_first_return=1)],
     [Orbit(num_prongs=2, orbit_length=2, num_orbits=2, order_after_first_return=1)],
     [Orbit(num_prongs=2, orbit_length=1, num_orbits=2, order_after_first_return=1),
      Orbit(num_prongs=2, orbit_length=2, num_orbits=1, order_after_first_return=1)],
     [Orbit(num_prongs=2, orbit_length=1, num_orbits=4,
        order_after_first_return=1)]]

    """
    result = []
    for pattern in orbit_patterns(multiplicity):
        result.extend(append_options(pattern, (num_prongs//2).divisors()))

    return [[Orbit(num_prongs,*item) for item in ls] for ls in result]


def orbit_combinations_in_stratum(stratum):
    for item in itertools.product(*[orbit_combinations_fixed_num_prongs(*x) for x in
    stratum]):
        yield list(itertools.chain(*item))

def lefschetz_contribution(orbit, power, lambda_pos=True):
    """
    INPUT:

    - ``orbit`` -- an Orbit
    - ``power`` -- the power of the mapping class to be considered

    TESTS:

    sage: lefschetz_contribution(18, 2, 9, 18)
    -34
    sage: lefschetz_contribution(18, 2, 9, 16)
    2
    sage: lefschetz_contribution(18, 2, 9, 17)
    0
    sage: lefschetz_contribution(18, 1, 9, 16)
    1
    sage: lefschetz_contribution(18, 1, 9, 17)
    1
    sage: lefschetz_contribution(18, 1, 9, 18)
    -17
    sage: lefschetz_contribution(10, 1, 5, 5, False)
    1
    """
    result = 0
    if power % orbit.orbit_length != 0:
        return 0
    stabilizer_power = power // orbit.orbit_length
    if lambda_pos or power%2 == 0:
        if stabilizer_power % orbit.order_after_first_return == 0:
            result = orbit.orbit_length * (1 - orbit.num_prongs)
        else:
            result = orbit.orbit_length
    else:
        result = orbit.orbit_length
    return orbit.num_orbits * result

# def total_lefschetz_contribution(num_prongs,orbit_combinations_fixed_num_prongs,k):
#     """
#     EXAMPLE:

#     sage: total_lefschetz_contribution([4,6], ([(2, 1, 1)], [(2, 1, 1)]), 1)

#     """
#     contribution = 0
#     for i in range(len(num_prongs)):
#         ops = orbit_combinations_fixed_num_prongs[i]
#         for op in ops:
#             contribution += op[1] * lefschetz_contribution(num_prongs[i], op[0],
#                                                            op[2], k)
#     return contribution



def max_num_singularities(genus):
    return 2*genus-2

def lefschetz_numbers(poly, max_power=10, is_orientable=True):
    mat = companion_matrix(poly)
    base = 2 if is_orientable else 1
    return [base-(mat**k).trace() for k in range(1, max_power+1)]

def reciprocalize(poly):
    """
    Multiply a polynomial by its reciprocal to make it reciprocal.
    """
    p = poly * poly.reverse()
    return p * sign(p[0])

def create_summary(recip_poly, orbits, max_power,
                   lambda_pos=True, debug=False, lef_numbers=None):
    """

    EXAMPLE:

    sage: create_summary((-1-x)^2*(1-x+x^2-x^3+x^4)*(1-x^2+2*x^5-x^8+x^10),
    [Orbit(18, 2, 1, 9)],20,False)

    sage: create_summary((-1-x)^2*(1-x+x^2-x^3+x^4)*(1-x^2+2*x^5-x^8+x^10), [Orbit(18,1, 2, 9)],20,False)

    sage: create_summary((-1+x)^4*(1+x)^2*(1-x+x^2-x^3+x^4-4*x^5+x^6-x^7+x^8-x^9+x^10), [Orbit(10,1, 1, 5), Orbit(6,5, 1, 3)], 20, True)
    """
    if debug:
        print_summary(recip_poly.factor(), orbits, lambda_pos)
    summary_array = np.zeros((3+len(orbits), max_power), dtype=int)
    if lef_numbers is not None:
        lef_nums = lef_numbers
    else:
        lef_nums = lefschetz_numbers(recip_poly,max_power,True)
    summary_array[0] = lef_nums

    for i in range(len(orbits)):
        orbit = orbits[i]
        for power in range(1,max_power+1):
            summary_array[i+1][power-1] = \
                lefschetz_contribution(orbit, power, lambda_pos)

    # filling in regular contributions
    for power in range(1, max_power+1):
        summary_array[-2][power-1] = summary_array[0][power-1]-\
                                     sum(summary_array[1:-2,power-1])

    # determining regular orbits
    for power in range(1, max_power+1):
        regular_contr = lefschetz_contribution(Orbit(2, power, 1, 1), power, lambda_pos)
        # print regular_contr
        previous_contr = 0
        for k in range(1, power):
            previous_contr += summary_array[-1][k-1] * \
                              lefschetz_contribution(Orbit(2, k, 1, 1), power, lambda_pos)
        diff = summary_array[-2][power-1]-previous_contr

        def print_data():
            print "================="
            print "Power: ", power
            print "Difference: ", diff
            print "Contribution of each regular orbit:", regular_contr
            print "================="
            print matrix(summary_array)
            print "================="

        if diff % regular_contr != 0:
            if debug:
                print_data()
            raise RuntimeError("Divisibility issue for regular orbits")
        quotient = diff // regular_contr
        if quotient < 0:
            if debug:
                print_data()
            raise RuntimeError("Lefschetz number cannot be fixed with "
                               "regular fixed points.")
        summary_array[-1][power-1] = diff / regular_contr
    return summary_array

def get_compatible_strata(poly, orbit_combinations, max_power=20):
    rec_poly = reciprocalize(poly)
    candidate_strata = []

    # testing all strata
    print "-------------------------------------"
    print "Testing ", poly

    lef_nums = lefschetz_numbers(rec_poly,max_power,True)
    max_lefschetz = max_num_singularities(poly.degree())
    for k in lef_nums:
        if k > max_lefschetz:
            # print k, "is too big Lefschetz number:"
            # print poly
            # break
            print "-------------------------------------"
            return []
    # print poly, poly.factor()


    for orbits in orbit_combinations:
        for lambda_pos in [True, False]:
            try:
                summary = create_summary(rec_poly, orbits,
                                         max_power, lambda_pos,
                                         lef_numbers=lef_nums)
                print_summary(rec_poly.factor(), orbits, lambda_pos)
                print matrix(summary)
            except RuntimeError as ex:
                continue
                # print str(ex)
            candidate_strata.append([orbits, lambda_pos])
    # print lef_nums
    print "-------------------------------------"
    return candidate_strata
# for poly in deg_8_candidates:
#     test(poly)
# test(deg_8_candidates[0])

def print_summary(poly, orbits, lambda_pos):
    print "Polynomial:", poly
    # print "Stratum:", stratum
    # print "Number of prongs:", num_prongs
    print "Orbit patterns:", orbits
    print "The stretch factor is: ", "POSITIVE" if lambda_pos else "NEGATIVE"


def get_possible_polynomials(degree, max_power=50):
    result = {}
    for poly in poly_candidates[degree]:
        orbit_combinations = [orbits for stratum in strata(poly.degree())
                              for orbits in orbit_combinations_in_stratum(stratum)]
        candidate_strata = get_compatible_strata(poly, orbit_combinations,
                                                 max_power)
        if len(candidate_strata) > 0:
            result[poly] = candidate_strata
    return result


# the result fo get_possible_polynomials(8) is:
# candidates_8_old ={x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1: [[[4], [(14, 1, 1)], True],
#                                                            [[4], [(14, 1, 2)], True],
#                                                            [[4], [(7, 2, 2)], True]],
#                x^8 - x^7 - x^5 + x^3 - x + 1: [[[4], [(14, 1, 1)], True],
#                                                [[4], [(14, 1, 2)], True],
#                                                [[4], [(10, 1, 2), (4, 1, 1)], True],
#                                                [[4], [(10, 1, 2), (4, 1, 2)], True],
#                                                [[4], [(10, 1, 2), (2, 2, 1)], True],
#                                                [[4], [(10, 1, 2), (2, 1, 1), (2, 1, 2)], True],
#                                                [[4], [(10, 1, 2), (2, 2, 2)], True],
#                                                [[4], [(7, 2, 2)], True]],
#                x^8 - x^5 - x^4 - x^3 - 1: [[[16], [(1, 2, 8)], True]]}

candidates_8 = {x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1:
                [[[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=1)],
                  True],
                 [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2)],
                  True],
                 [[Orbit(num_prongs=4, orbit_length=7, num_orbits=2, order_after_first_return=2)],
                  True]],
                x^8 - x^7 - x^5 + x^3 - x + 1:
                [[[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=1)],
                  True],
                 [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2)],
                  True],
                 [[Orbit(num_prongs=4, orbit_length=7, num_orbits=2, order_after_first_return=2)],
                  True]],
                x^8 - x^5 - x^4 - x^3 - 1:
                [[[Orbit(num_prongs=16, orbit_length=1, num_orbits=2, order_after_first_return=8)],
                  True]]}


candidates_10 = {x^10 - x^9 - x^6 + x^4 - x + 1:
                  [[[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=1)],
                   True],
                  [[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=2)],
                   True],
                  [[Orbit(num_prongs=4, orbit_length=9, num_orbits=2, order_after_first_return=2)],
                   True]],
                 x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1:
                 [[[Orbit(num_prongs=8, orbit_length=3, num_orbits=2, order_after_first_return=4)],
                   True]]}

candidates_12 = {
    x^12 - 3*x^11 + 2*x^10 + 2*x^9 - 3*x^8 + x^7 - x^5 + 3*x^4 - 4*x^3 + 4*x^2
    - 3*x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - 3*x^11 + 3*x^10 - x^9 - x^8 + 2*x^7 - 2*x^5 + x^4 - x^3 + 3*x^2 -
    3*x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True]],

    x^12 - 2*x^11 + x^9 + x^8 - 2*x^7 + 2*x^6 - x^4 - x^3 + 2*x^2 - 2*x + 1:
    [[[Orbit(num_prongs=24, orbit_length=2, num_orbits=1, order_after_first_return=12)],
      True],
     [[Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5),
       Orbit(num_prongs=16, orbit_length=2, num_orbits=1, order_after_first_return=8)],
      True],
     [[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=12, orbit_length=4, num_orbits=1, order_after_first_return=6)],
      True],
     [[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=8, orbit_length=4, num_orbits=1, order_after_first_return=4),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=1),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3),
       Orbit(num_prongs=8, orbit_length=6, num_orbits=1, order_after_first_return=4)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=8, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=8, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - 2*x^11 + x^10 - x^8 + x^7 + x^5 - x^4 - x^2 + 2*x - 1:
    [[[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=10, orbit_length=2, num_orbits=1, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=4, orbit_length=10, num_orbits=1, order_after_first_return=1),
       Orbit(num_prongs=6, orbit_length=6, num_orbits=1, order_after_first_return=3)],
      True],
     [[Orbit(num_prongs=4, orbit_length=10, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=6, orbit_length=6, num_orbits=1, order_after_first_return=3)],
      True],
     [[Orbit(num_prongs=4, orbit_length=5, num_orbits=2, order_after_first_return=2),
       Orbit(num_prongs=6, orbit_length=6, num_orbits=1, order_after_first_return=3)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - 2*x^11 + x^10 - x^2 + 1:
    [[[Orbit(num_prongs=24, orbit_length=2, num_orbits=1, order_after_first_return=12)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - 2*x^10 + 2*x^9 + x^8 - x^7 - x^5 + x^4 + x - 1:
    [[[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
       Orbit(num_prongs=22, orbit_length=2, num_orbits=1, order_after_first_return=11)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^10 + x^9 - x^8 + 2*x^6 - x^4 - x^3 + x^2 + x - 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^10 + x^9 - x^3 + x^2 - x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^9 + x^8 - x^4 + x^3 - x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - 2*x^3 + 2*x^2 - x + 1:
    [[[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=22, orbit_length=2, num_orbits=1, order_after_first_return=11)],
      True]],

    x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^7 + x^5 - x + 1:
    [[[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=1)],
      True],
     [[Orbit(num_prongs=4, orbit_length=22, num_orbits=1, order_after_first_return=2)],
      True],
     [[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - x^11 - x^6 + x - 1:
    [[[Orbit(num_prongs=4, orbit_length=11, num_orbits=2, order_after_first_return=2)],
      True]],

    x^12 - 2*x^10 - x^9 + x^8 + 2*x^5 + x^4 - x^3 - 1:
    [[[Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5),
       Orbit(num_prongs=16, orbit_length=2, num_orbits=1, order_after_first_return=8)],
      True],
     [[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5),
       Orbit(num_prongs=14, orbit_length=2, num_orbits=1, order_after_first_return=7)],
      True],
     [[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=8, orbit_length=4, num_orbits=1, order_after_first_return=4),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=1),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True],
     [[Orbit(num_prongs=4, orbit_length=14, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
      True]],

    x^12 - x^9 - x^8 - 2*x^7 + x^4 + x^3 + 2*x^2 + 1:
    [[[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
       Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5),
       Orbit(num_prongs=14, orbit_length=2, num_orbits=1, order_after_first_return=7)],
      True]],

    x^12 - 2*x^7 - 2*x^6 + 1:
    [[[Orbit(num_prongs=8, orbit_length=4, num_orbits=1, order_after_first_return=4),
       Orbit(num_prongs=12, orbit_length=1, num_orbits=2, order_after_first_return=6)],
      True]],

    x^12 - x^7 - 2*x^6 - x^5 - 1:
    [[[Orbit(num_prongs=24, orbit_length=1, num_orbits=2, order_after_first_return=12)],
      True]],

    x^12 - x^7 - x^6 - x^5 - 1:
    [[[Orbit(num_prongs=24, orbit_length=1, num_orbits=2, order_after_first_return=12)],
      True]]}

# {x^10 - x^9 - x^6 + x^4 - x + 1: [[[8], [(6, 1, 4)], True],
#                                                   [[4], [(18, 1, 1)], True],
#                                                   [[4], [(18, 1, 2)], True],
#                                                   [[4], [(16, 1, 2), (2, 1, 1)], True],
#                                                   [[4], [(16, 1, 2), (2, 1, 2)], True],
#                                                   [[4], [(16, 1, 2), (1, 2, 1)], True],
#                                                   [[4], [(16, 1, 2), (1, 2, 2)], True],
#                                                   [[4], [(12, 1, 2), (6, 1, 1)], True],
#                                                   [[4], [(12, 1, 2), (6, 1, 2)], True],
#                                                   [[4], [(9, 2, 2)], True]],
#                  x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1: [[[8], [(3, 2, 4)], True]]}



def double_check_candidates(candidates, max_power=100, debug=False):
    for poly in candidates.keys():
        rec_poly = reciprocalize(poly)
        for data in candidates_8[poly]:
            orbits, lambda_pos = data
            try:
                summary = create_summary(rec_poly, orbits, max_power,
                                         lambda_pos, debug)
                print "-------------------------------------------"
                print poly, data, "SEEMS OKAY."
                print "-------------------------------------------"
                print_summary(rec_poly, orbits, lambda_pos)
                print summary
            except RuntimeError as ex:
                print "-------------------------------------------"
                print poly, data, "HAS FAILED!!!!!!!"
                print "-------------------------------------------"

# double_check_candidates(candidates_8, 20)

# def double_check_candidate(candidates, max_power=100):
#     result = {}
#     for poly in poly_candidates[degree]:
#         candidate_strata = get_compatible_strata(poly)
#         if len(candidate_strata) > 0:
#             result[poly] = candidate_strata

#                 try:
#                     summary = create_summary(rec_poly, num_prongs, orbits, max_power, lambda_pos)
#                 except RuntimeError as ex:
#                     continue
#                     # print str(ex)
#                 print "Stratum:", stratum
#                 print "Number of prongs:", num_prongs
#                 print "Orbit patterns:", orbits
#                 print "The stretch factor is: ", "POSITIVE" if lambda_pos else "NEGATIVE"
#                 print summary

#     return result
