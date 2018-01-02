"""
This file contains some data.

- ``upper_bounds`` -- This contains efficient upper bounds for the minimal stretch factors. (That is, there exists a stretch factor smaller than this upper bound, but hopefully not many.)

- ``results`` -- the polynomials whose Perron root is at most the upper bound. The data here are results of running the polynomial search algorithm earlier. Other than storing the results at one place, another purpose of this database is to be able to test the polynomial search code in tests.py.

"""

upper_bounds = {'nonor': {3: 1.84, 4: 1.52, 5: 1.43, 6: 1.422, 7: 1.2885, 8:
                1.3568, 9: 1.2173, 10: 1.22262, 11: 1.1743, 12:1.2764, 13:
                1.145507, 14: 1.1875, 15:1.1249, 16:1.1426, 17:1.10939,
                18:1.20515, 19:1.097305, 21:1.087629, 23:1.079704},
                'reversing': {
                    4: 1.62, 6: 1.253, 8: 1.4, 10: 1.2, 12: 1.16
                },
                'preserving': {
                    4: 1.73, 6: 1.402, 8: 1.281, 10: 1.177, 12: 1.177, 14: 1.116, 16: 1.129
                }
}

from sage.all import *
x = var('x')
results = {'nonor': 
{3: [(x^3 - x^2 - x - 1, 1.839286755214161)],
4: [(x^4 - x^3 - x^2 + x - 1, 1.512876396864095)],
5: [(x^5 - x^3 - x^2 - 1, 1.429108319838146)],
6: [(x^6 - x^5 - x^3 + x - 1, 1.421975014306898)],
7: [(x^7 - x^4 - x^3 - 1, 1.288452726276414)],
8: [(x^8 - x^7 + x^6 - x^5 - x^4 + x^3 - x^2 + x - 1, 1.237428905126040),
(x^8 - x^7 - x^5 + x^3 - x + 1, 1.288452726276414),
(x^8 + x^7 - x^5 - 2*x^4 - x^3 - x - 1, 1.288452726276414),
(x^8 - x^6 - x^5 - x^3 + x^2 + 1, 1.289444523350383),
(x^8 + 2*x^7 + x^6 - x^5 - 2*x^4 - 3*x^3 - 3*x^2 - 2*x - 1, 1.289444523350383),
(x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1, 1.307395463062887),
(x^8 + x^7 - x^6 - x^5 - x^3 - x^2 - x - 1, 1.307395463062887),
(x^8 - 2*x^5 - 1, 1.311809917105982),
(x^8 - x^7 - 2*x^5 + 2*x^4 - x + 1, 1.324717957244746),
(x^8 - x^7 - 2*x^3 + x + 1, 1.324717957244746),
(x^8 - x^6 - x^5 - x^4 - x^3 + x^2 + 2*x + 1, 1.324717957244746),
(x^8 - x^6 - x^5 - x^4 + x^3 + x^2 - 1, 1.324717957244746),
(x^8 + x^7 - 2*x^5 - 2*x^4 - x - 1, 1.324717957244746),
(x^8 + x^7 - 2*x^3 - 4*x^2 - 3*x - 1, 1.324717957244746),
(x^8 - x^6 - x^5 + x^3 - x^2 - 1, 1.344489120817385),
(x^8 + x^7 - x^5 - 3*x^4 - x^3 - x - 1, 1.350982326113995),
(x^8 + x^6 - x^4 - 2*x^3 - 3*x^2 - 2*x - 1, 1.355187804953756),
(x^8 - x^5 - x^4 - x^3 - 1, 1.356797155088499)],
9: [(x^9 - x^5 - x^4 - 1, 1.217281181744868)],
11: [(x^11 - x^6 - x^5 - 1, 1.174290884852624)],
13: [(x^13 - x^7 - x^6 - 1, 1.145506428012675)]},
'reversing': {
4: [(x^4 - x^3 - x - 1, 1.618033988749895)],
6: [(x^6 - x^4 - x^3 + x^2 - 1, 1.252073078214420)],
8: [(x^8 - x^5 - x^3 - 1, 1.252073078214420),
 (x^8 - x^7 - x - 1, 1.324717957244746),
 (x^8 - x^6 - x^5 - x^3 + x^2 - 1, 1.397758813869009)],
10: [(x^10 - x^8 + x^6 - x^5 - x^4 + x^2 - 1, 1.159731372540720),
 (x^10 + x^8 - x^7 - x^5 - x^3 - x^2 - 1, 1.193859111321223)],
12: [(x^12 - x^7 - x^5 - 1, 1.159731372540720)],

},
'preserving': {
4: [(x^4 - x^3 - x^2 - x + 1, 1.722083805739043)],
6: [(x^6 + x^5 - x^4 - 3*x^3 - x^2 + x + 1, 1.324717957244746),
 (x^6 - x^4 - x^3 - x^2 + 1, 1.401268367939855)],
8: [(x^8 - x^5 - x^4 - x^3 + 1, 1.280638156267758)],
10: [(x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, 1.176280818259918)],
12: [(x^12 - x^11 - x^10 + x^8 + x^4 - x^2 - x + 1, 1.176280818259918),
 (x^12 - x^7 - x^6 - x^5 + 1, 1.176280818259918),
 (x^12 + x^11 + x^10 - x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 + x^2 + x + 1,
  1.176280818259918),
 (x^12 + 2*x^11 + 2*x^10 - 2*x^8 - 3*x^7 - 3*x^6 - 3*x^5 - 2*x^4 + 2*x^2 + 2*x + 1,
  1.176280818259918),
 (x^12 + 3*x^11 + 3*x^10 - 3*x^8 - 4*x^7 - 4*x^6 - 4*x^5 - 3*x^4 + 3*x^2 + 3*x + 1,
  1.176280818259918)],
14: [(x^14 + x^13 + x^12 - x^10 - x^9 - x^8 - x^7 - x^6 - x^5 - x^4 + x^2 + x + 1,
  1.105400852650789),
 (x^14 + x^13 - x^9 - x^8 - x^7 - x^6 - x^5 + x + 1, 1.115481109456592)]}
}