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

10: [(x^10 - x^7 - x^6 - x^5 + x^4 + x^3 - 1, 1.2037667118876316),
(x^10 - x^9 + x^8 - x^7 - x^5 + x^3 - x^2 + x - 1, 1.217150549717164),
(x^10 + x^9 - x^6 - 2*x^5 - x^4 - x - 1, 1.2172811817448663),
(x^10 - x^9 - x^6 + x^4 - x + 1, 1.2172811817448697),
(x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1, 1.22261519674945)],

11: [(x^11 - x^6 - x^5 - 1, 1.174290884852624)],

12: [(x^12 - x^11 + x^10 - x^9 + x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.1543552117091553),
(x^12 - 2*x^9 - x^6 + 2*x^3 + 1, 1.1739850106977299),
(x^12 - x^11 - x^7 + x^5 - x + 1, 1.1742908848526246),
(x^12 + x^11 - x^7 - 2*x^6 - x^5 - x - 1, 1.1742908848526257),
(x^12 + x^11 - x^8 - x^7 - x^5 - x^4 - x - 1, 1.178375014855083),
(x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - x + 1, 1.1783750148550873),
(x^12 + 2*x^11 + 2*x^10 + x^9 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 3*x^3 - 2*x^2 - 2*x - 1, 1.1785964146670196),
(x^12 - x^9 - x^7 - x^5 + x^3 + 1, 1.1785964146670211),
(x^12 - 2*x^7 - 1, 1.1822598571242353),
(x^12 + x^11 - 2*x^7 - 2*x^6 - x - 1, 1.187093002902818),
(x^12 - x^11 - 2*x^7 + 2*x^6 - x + 1, 1.1870930029028193),
(x^12 + x^11 - x^9 - x^8 - x^4 - x^3 - x - 1, 1.18773435095945),
(x^12 - x^11 - x^9 + x^8 - x^4 + x^3 - x + 1, 1.1877343509594531),
(x^12 - x^11 - x^10 + 2*x^9 - x^8 + x^6 - 2*x^5 + x^4 - x^2 + x - 1, 1.1907127070817989),
(x^12 + x^11 - x^10 - 2*x^9 - x^8 + x^6 - x^4 + x^2 + x + 1, 1.1907127070817989),
(x^12 - x^10 - x^9 + x^8 + x^7 - x^6 - x^5 - x^4 + x^3 + x^2 - 1, 1.1969835607261106),
(x^12 - x^11 + x^10 - x^9 - x^6 + x^3 - x^2 + x - 1, 1.1990828456262594),
(x^12 - x^11 + x^9 - 2*x^7 - x^3 + x + 1, 1.2032160335179485),
(x^12 - x^11 - x^9 + 2*x^8 - 2*x^7 + 2*x^4 - x^3 - x + 1, 1.2032160335179505),
(x^12 + x^11 - x^9 - 2*x^6 - 2*x^5 + x^3 - x - 1, 1.2032160335179505),
(x^12 + x^11 + x^9 + 2*x^8 - 2*x^6 - 2*x^5 - 2*x^4 - 3*x^3 - 4*x^2 - 3*x - 1, 1.203216033517951),
(x^12 + x^11 + x^10 - x^9 - 2*x^8 - 3*x^7 - x^6 + x^5 + 2*x^4 + x^3 - x^2 - x - 1, 1.2037667118876272),
(x^12 - 2*x^11 + x^10 - x^9 + x^8 + 2*x^6 - 2*x^5 - x^4 + x^3 - x^2 + 2*x - 1, 1.203766711887628),
(x^12 + x^10 - x^9 - x^8 - 2*x^7 + x^4 + x^3 - x^2 - 1, 1.2037667118876292),
(x^12 - x^10 - x^9 - x^8 + 2*x^6 + 2*x^5 - x^4 - x^3 - x^2 + 1, 1.2037667118876312),
(x^12 - x^11 + x^10 - x^9 - x^7 + x^6 - x^5 + x^3 - x^2 + x - 1, 1.2037667118876318),
(x^12 + 2*x^11 + x^10 - x^9 - 3*x^8 - 4*x^7 - 2*x^6 + 2*x^5 + 3*x^4 + x^3 - x^2 - 2*x - 1, 1.2037667118876338),
(x^12 - x^11 + x^9 - x^8 - x^6 + x^4 - x^3 + x - 1, 1.2046221387845386),
(x^12 + x^11 - x^10 - x^9 - x^3 - x^2 - x - 1, 1.2060690022944434),
(x^12 - x^11 - x^10 + x^9 - x^3 + x^2 - x + 1, 1.2060690022944445),
(x^12 - x^11 + x^9 - 2*x^8 + x^6 - 2*x^5 + x^3 + x + 1, 1.2065470532989873),
(x^12 - x^11 + x^9 - 2*x^8 + 2*x^7 - x^6 - 2*x^5 + 2*x^4 - x^3 + x - 1, 1.206547053298997),
(x^12 - 2*x^7 - x^6 + 1, 1.2068773951849465),
(x^12 - x^11 + x^8 - x^7 - x^6 + x^5 - x^4 + x - 1, 1.2086326843701178),
(x^12 - x^11 + x^8 - 2*x^7 + x^6 - x^4 + x - 1, 1.214243002891164),
(x^12 - x^11 + x^10 - x^9 - x^8 + x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.2142981907580435),
(x^12 - x^10 - x^9 + x^7 + x^6 - x^5 - 2*x^4 + x^3 + x^2 - 1, 1.2169048242672886),
(x^12 - 3*x^11 + 4*x^10 - 4*x^9 + 3*x^8 - 2*x^7 + 2*x^6 - 3*x^4 + 4*x^3 - 4*x^2 + 3*x - 1, 1.2171505497171526),
(x^12 - x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^7 - x^4 + 2*x^3 - 2*x^2 + x - 1, 1.217150549717162),
(x^12 + x^11 - x^8 - 2*x^7 - 2*x^6 + x^4 - x - 1, 1.2171505497171635),
(x^12 - x^11 - x^8 + 2*x^5 - x^4 - x + 1, 1.2171505497171644),
(x^12 - 2*x^11 + 3*x^10 - 3*x^9 + 2*x^8 - 2*x^7 + x^6 - 2*x^4 + 3*x^3 - 3*x^2 + 2*x - 1, 1.2171505497171662),
(x^12 + x^10 - x^9 - 2*x^7 - x^6 + x^3 - x^2 - 1, 1.2171505497171682),
(x^12 - 3*x^11 + 3*x^10 - x^9 - x^8 + 2*x^7 - 2*x^5 + x^4 - x^3 + 3*x^2 - 3*x + 1, 1.217281181744849),
(x^12 + x^11 + x^10 + x^9 - x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^3 - x^2 - x - 1, 1.2172811817448657),
(x^12 + 3*x^11 + 3*x^10 + x^9 - x^8 - 4*x^7 - 6*x^6 - 4*x^5 - x^4 - x^3 - 3*x^2 - 3*x - 1, 1.2172811817448657),
(x^12 - x^9 - x^8 - x^7 + x^5 + x^4 - x^3 + 1, 1.2172811817448659),
(x^12 + 2*x^11 + 2*x^10 + x^9 - x^8 - 3*x^7 - 4*x^6 - 3*x^5 - x^4 - x^3 - 2*x^2 - 2*x - 1, 1.2172811817448665),
(x^12 - x^11 + x^10 - x^9 - x^8 + x^4 - x^3 + x^2 - x + 1, 1.217281181744867),
(x^12 + x^9 - x^8 - x^7 - x^5 - x^4 - x^3 - 1, 1.2172811817448672),
(x^12 + x^11 - x^10 - x^9 - x^8 - 2*x^7 + 2*x^5 + x^4 - x^3 - x^2 + x + 1, 1.2172811817448692),
(x^12 - 2*x^11 + 2*x^10 - x^9 - x^8 + x^7 - x^5 + x^4 - x^3 + 2*x^2 - 2*x + 1, 1.2172811817448694),
(x^12 - x^11 - x^10 + x^9 - x^8 + 2*x^6 - x^4 - x^3 + x^2 + x - 1, 1.217281181744873),
(x^12 - x^10 + x^9 - 2*x^5 - x^3 - x^2 - 1, 1.2224521787260243),
(x^12 - 2*x^11 + x^10 + x^9 - 2*x^8 + 2*x^7 - 2*x^6 + 2*x^4 - 3*x^3 + 3*x^2 - 2*x + 1, 1.222452178726025),
(x^12 - x^7 - x^6 - x^5 - 1, 1.2226151967494492),
(x^12 + x^11 - x^10 + x^8 - 2*x^7 - 2*x^6 - x^4 + x^2 - x - 1, 1.2226151967494494),
(x^12 - 3*x^11 + 3*x^10 - 3*x^8 + 2*x^7 + 2*x^6 - 4*x^5 + 3*x^4 - 3*x^2 + 3*x - 1, 1.2226151967494496),
(x^12 - x^11 + x^10 - x^8 - 2*x^5 + x^4 - x^2 + x - 1, 1.2226151967494518),
(x^12 - x^11 - x^10 + 2*x^9 - x^8 - 2*x^7 + 2*x^6 - x^4 + 2*x^3 - x^2 - x + 1, 1.2226151967494523),
(x^12 - 2*x^11 + 2*x^10 - 2*x^8 + x^7 + x^6 - 3*x^5 + 2*x^4 - 2*x^2 + 2*x - 1, 1.2226151967494532),
(x^12 - 3*x^11 + 4*x^10 - 4*x^9 + 2*x^8 + 2*x^7 - 4*x^6 + 4*x^5 - 4*x^4 + 4*x^3 - 4*x^2 + 3*x - 1, 1.223819577363795),
(x^12 + x^10 - x^9 - x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 - 1, 1.2238195773638108),
(x^12 - 2*x^11 + 3*x^10 - 3*x^9 + x^8 + x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1, 1.2238195773638112),
(x^12 + x^11 - 2*x^8 - 2*x^7 - x - 1, 1.223819577363812),
(x^12 - x^11 - 2*x^8 + 2*x^7 - x + 1, 1.2238195773638123),
(x^12 - x^11 + 2*x^10 - 2*x^9 - 2*x^6 + 2*x^5 - 2*x^4 + 2*x^3 - 2*x^2 + x - 1, 1.2238195773638127),
(x^12 + x^11 - x^10 - x^9 - x^7 + x^5 - x^3 - x^2 - x - 1, 1.2253771855899727),
(x^12 - x^11 - x^10 + x^9 - x^7 + 2*x^6 - x^5 - x^3 + x^2 - x + 1, 1.225377185589974),
(x^12 - x^11 - x^6 + 2*x^5 - 2*x^4 + x - 1, 1.2270982484786555),
(x^12 + x^11 - x^9 - x^8 - 2*x^7 - x^6 + x^4 + x^3 - x - 1, 1.2300632025017078),
(x^12 - x^9 - x^8 - x^6 + x^4 + x^3 - 1, 1.2301520433430366),
(x^12 - x^8 - x^7 + x^6 - x^5 - x^4 - 1, 1.232573739469347),
(x^12 - 3*x^11 + 3*x^10 - x^9 + x^7 - 4*x^6 + 5*x^5 - 2*x^4 + x^3 - 3*x^2 + 3*x - 1, 1.2333284309598),
(x^12 - x^11 + x^10 - x^9 + x^7 - 2*x^6 + x^5 - 2*x^4 + x^3 - x^2 + x - 1, 1.2333284309598074),
(x^12 + x^11 - x^10 - x^9 + x^7 - 3*x^5 - 2*x^4 + x^3 + x^2 - x - 1, 1.2333284309598074),
(x^12 - x^11 - x^10 + x^9 + x^7 - 2*x^6 - x^5 + 2*x^4 + x^3 - x^2 - x + 1, 1.2333284309598096),
(x^12 - x^9 + x^7 - x^6 - x^5 - 2*x^4 + x^3 - 1, 1.2333284309598112),
(x^12 - 2*x^11 + 2*x^10 - x^9 + x^7 - 3*x^6 + 3*x^5 - 2*x^4 + x^3 - 2*x^2 + 2*x - 1, 1.2333284309598118),
(x^12 - x^11 + x^10 - 2*x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^2 + x - 1, 1.2339551780878162),
(x^12 - x^11 - x^7 + x^6 + x^5 - 2*x^4 + x - 1, 1.2341997914505352),
(x^12 - 2*x^11 + x^9 + x^8 - 2*x^7 + 2*x^6 - x^4 - x^3 + 2*x^2 - 2*x + 1, 1.2365057033914317),
(x^12 - 2*x^11 + 3*x^9 - 3*x^8 + 2*x^6 - 2*x^5 + x^4 + x^3 - 2*x^2 + 2*x - 1, 1.2365057033914924),
(x^12 - 2*x^10 - x^9 + x^8 + 2*x^5 + x^4 - x^3 - 1, 1.2365057033914926),
(x^12 - x^9 - x^8 - 2*x^7 + x^4 + x^3 + 2*x^2 + 1, 1.2365057033914941),
(x^12 + 2*x^11 - 3*x^9 - 3*x^8 - 2*x^7 - 2*x^6 + 3*x^4 + 3*x^3 + 2*x^2 + 2*x + 1, 1.2365057033914961),
(x^12 + x^9 - x^8 - 2*x^5 - x^4 - x^3 - 2*x^2 - 1, 1.236505703391497),
(x^12 + x^11 - x^6 - 2*x^5 - 2*x^4 - 2*x^3 - 2*x^2 - x - 1, 1.2365057033914981),
(x^12 + 2*x^11 - x^9 + x^8 - 2*x^6 - 2*x^5 - 3*x^4 - 3*x^3 - 2*x^2 - 2*x - 1, 1.2365057033915),
(x^12 - x^11 + 2*x^9 - 2*x^8 + x^6 - 2*x^5 - 2*x^2 + x - 1, 1.2365057033915001),
(x^12 - 2*x^10 + x^9 + x^8 - 2*x^7 - x^4 + x^3 + 1, 1.2365057033915003),
(x^12 - x^11 - 2*x^7 + x^6 + 2*x^2 - x + 1, 1.236505703391501),
(x^12 + x^11 - 2*x^9 - 2*x^8 - 2*x^7 - x^6 + 2*x^4 + 2*x^3 + 2*x^2 + x + 1, 1.236505703391505),
(x^12 - 4*x^11 + 8*x^10 - 11*x^9 + 10*x^8 - 4*x^7 - 4*x^6 + 10*x^5 - 12*x^4 + 11*x^3 - 8*x^2 + 4*x - 1, 1.237428905126001),
(x^12 - 3*x^11 + 3*x^10 - x^9 - 2*x^8 + 6*x^7 - 6*x^6 + 2*x^5 - x^3 + 3*x^2 - 3*x + 1, 1.237428905126009),
(x^12 - x^11 + 3*x^10 - 3*x^9 + 2*x^8 - 2*x^7 - 2*x^6 + 2*x^5 - 4*x^4 + 3*x^3 - 3*x^2 + x - 1, 1.2374289051260368),
(x^12 + 3*x^11 + 3*x^10 + x^9 - 2*x^8 - 6*x^7 - 6*x^6 - 2*x^5 - x^3 - 3*x^2 - 3*x - 1, 1.2374289051260376),
(x^12 - x^11 - x^8 + x^7 + x^6 - x^5 - x^4 + x - 1, 1.2374289051260385),
(x^12 - 2*x^11 + 2*x^10 - x^9 - 2*x^8 + 4*x^7 - 4*x^6 + 2*x^5 - x^3 + 2*x^2 - 2*x + 1, 1.237428905126039),
(x^12 - x^11 + 2*x^10 - 2*x^9 + x^8 - x^7 - x^6 + x^5 - 3*x^4 + 2*x^3 - 2*x^2 + x - 1, 1.237428905126039),
(x^12 - x^11 - x^10 + x^9 - 2*x^8 + 2*x^7 + 2*x^6 - 2*x^5 - x^3 + x^2 + x - 1, 1.2374289051260396),
(x^12 + x^11 - x^10 - x^9 - 2*x^8 - 2*x^7 + 2*x^6 + 2*x^5 - x^3 - x^2 + x + 1, 1.2374289051260396),
(x^12 + x^10 - x^8 - x^7 - x^6 - x^5 - x^4 - x^2 - 1, 1.2374289051260399),
(x^12 + 2*x^11 + 2*x^10 + x^9 - 2*x^8 - 4*x^7 - 4*x^6 - 2*x^5 - x^3 - 2*x^2 - 2*x - 1, 1.2374289051260403),
(x^12 - 2*x^11 + 4*x^10 - 5*x^9 + 4*x^8 - 2*x^7 - 2*x^6 + 4*x^5 - 6*x^4 + 5*x^3 - 4*x^2 + 2*x - 1, 1.237428905126041),
(x^12 - x^9 - 2*x^8 + 2*x^5 - x^3 + 1, 1.237428905126041),
(x^12 + x^11 + x^10 + x^9 - 2*x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^3 - x^2 - x - 1, 1.237428905126041),
(x^12 + 2*x^10 - x^9 - 2*x^7 - 2*x^6 - 2*x^4 + x^3 - 2*x^2 - 1, 1.2374289051260412),
(x^12 - 3*x^11 + 6*x^10 - 8*x^9 + 7*x^8 - 3*x^7 - 3*x^6 + 7*x^5 - 9*x^4 + 8*x^3 - 6*x^2 + 3*x - 1, 1.2374289051260419),
(x^12 - x^11 + x^10 - x^9 - 2*x^8 + 2*x^7 - 2*x^6 + 2*x^5 - x^3 + x^2 - x + 1, 1.2374289051260423),
(x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 3*x^8 - x^7 - x^6 + 3*x^5 - 5*x^4 + 4*x^3 - 3*x^2 + 2*x - 1, 1.2374289051260425),
(x^12 - x^11 + x^10 - x^9 - 2*x^4 + x^3 - x^2 + x - 1, 1.2374289051260425),
(x^12 + x^9 - 2*x^8 - 2*x^5 - x^3 - 1, 1.2374289051260425),
(x^12 - 2*x^11 + 2*x^10 - 3*x^9 + 2*x^8 + 2*x^5 - 4*x^4 + 3*x^3 - 2*x^2 + 2*x - 1, 1.2374289051260432),
(x^12 + x^11 + 2*x^10 - x^8 - 3*x^7 - 3*x^6 - x^5 - x^4 - 2*x^2 - x - 1, 1.2374289051260432),
(x^12 - 3*x^11 + 5*x^10 - 7*x^9 + 6*x^8 - 2*x^7 - 2*x^6 + 6*x^5 - 8*x^4 + 7*x^3 - 5*x^2 + 3*x - 1, 1.2374289051260532),
(x^12 - 5*x^11 + 11*x^10 - 15*x^9 + 14*x^8 - 6*x^7 - 6*x^6 + 14*x^5 - 16*x^4 + 15*x^3 - 11*x^2 + 5*x - 1, 1.2374289051260647),
(x^12 - x^9 - x^7 - 2*x^6 - x^5 + x^3 + 2*x^2 + 2*x + 1, 1.2377747004221065),
(x^12 + x^11 - 2*x^10 - 3*x^9 + x^8 + 3*x^7 - x^5 - x^4 - x^3 + x + 1, 1.2381636732853454),
(x^12 + x^11 - x^9 - x^8 - x^7 + x^5 - x^4 - x^3 - 2*x^2 - x - 1, 1.2381636732853454),
(x^12 - x^11 - 2*x^10 + x^9 + 3*x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 + x - 1, 1.2381636732853458),
(x^12 + 3*x^11 + 2*x^10 - 3*x^9 - 5*x^8 - x^7 + 2*x^6 + x^5 - x^4 - 3*x^3 - 4*x^2 - 3*x - 1, 1.2381636732853458),
(x^12 - x^10 + x^8 - x^7 - x^6 + x^5 - x^4 - x^2 - 1, 1.2381636732853485),
(x^12 + 2*x^11 + x^10 - 2*x^9 - 3*x^8 - x^7 + x^6 + x^5 - x^4 - 2*x^3 - 3*x^2 - 2*x - 1, 1.2381636732853485),
(x^12 - x^10 + x^9 - 2*x^7 + x^6 - 2*x^4 - x^3 + x^2 - 1, 1.238721769655264),
(x^12 + x^11 + x^8 - 2*x^6 - 2*x^5 - 3*x^4 - 2*x^3 - 2*x^2 - x - 1, 1.2393378079735542),
(x^12 + x^11 - x^9 - 2*x^8 - x^7 + x^5 - x^3 - x - 1, 1.2395129108220568),
(x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 4*x^8 - 4*x^7 + 3*x^6 - 2*x^5 + 2*x^3 - 3*x^2 + 2*x - 1, 1.2395129108220586),
(x^12 + x^10 - 2*x^9 - 2*x^7 + x^6 - x^2 - 1, 1.2395129108220595),
(x^12 - x^11 + 2*x^10 - 3*x^9 + 2*x^8 - 3*x^7 + 2*x^6 - x^5 + x^3 - 2*x^2 + x - 1, 1.2395129108220597),
(x^12 - x^11 - x^9 + x^7 + x^5 - 2*x^4 + x^3 - x + 1, 1.2395129108220624),
(x^12 - 3*x^11 + 4*x^10 - 5*x^9 + 6*x^8 - 5*x^7 + 4*x^6 - 3*x^5 + 3*x^3 - 4*x^2 + 3*x - 1, 1.2395129108220813),
(x^12 - x^9 - 2*x^7 + x^6 - x^3 + 1, 1.2400741298418678),
(x^12 - x^8 - 2*x^7 + x^4 - 1, 1.240530070122331),
(x^12 + x^11 + x^10 - x^9 - 2*x^8 - 3*x^7 - x^6 - x^5 + x^3 + x^2 + x + 1, 1.2406843840135455),
(x^12 - x^10 - x^7 - x^6 + x^5 + x^2 - 1, 1.2417610053620818),
(x^12 - x^8 - x^7 - x^6 - x^5 + x^4 + 1, 1.241877844116519),
(x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 3*x^8 - 3*x^7 + x^6 - x^5 + x^4 - 2*x^2 + 2*x - 1, 1.2427954473188327),
(x^12 - x^11 - x^10 + x^9 + x^8 - 2*x^7 + x^4 + x^3 - x^2 - x + 1, 1.2444237118811827),
(x^12 + x^11 - x^10 - x^9 + x^8 - 2*x^6 - 2*x^5 - x^4 + x^3 + x^2 - x - 1, 1.2444237118811858),
(x^12 - x^11 - 2*x^10 + 2*x^9 + x^8 - x^7 - x^5 + x^4 + x - 1, 1.2449941133546525),
(x^12 - x^10 - x^9 + x^7 - x^5 - x^3 + x^2 + 1, 1.2449941133546534),
(x^12 - x^11 - x^8 + x^7 - x^5 + x^4 - 2*x^3 + 2*x^2 - x + 1, 1.2449941133546558),
(x^12 + x^11 - x^8 - x^7 - x^5 - x^4 - 2*x^3 - 2*x^2 - x - 1, 1.2449941133546565),
(x^12 + 2*x^11 + x^10 - x^9 - 2*x^8 - x^7 - x^5 - 2*x^4 - 3*x^3 - 3*x^2 - 2*x - 1, 1.2449941133546574),
(x^12 - 2*x^11 + x^10 + x^9 - 2*x^8 + x^7 - x^5 + 2*x^4 - 3*x^3 + 3*x^2 - 2*x + 1, 1.2449941133546576),
(x^12 - x^10 + x^9 - x^7 - x^5 - x^3 - x^2 - 1, 1.2449941133546583),
(x^12 + 3*x^11 + 2*x^10 - 2*x^9 - 3*x^8 - x^7 - x^5 - 3*x^4 - 4*x^3 - 4*x^2 - 3*x - 1, 1.2449941133546596),
(x^12 + x^11 - 2*x^10 - 2*x^9 + x^8 + x^7 - x^5 - x^4 + x + 1, 1.2449941133546605),
(x^12 - 3*x^11 + 2*x^10 + 2*x^9 - 3*x^8 + x^7 - x^5 + 3*x^4 - 4*x^3 + 4*x^2 - 3*x + 1, 1.244994113354706),
(x^12 - x^9 - x^7 - x^5 + x^3 - 1, 1.245328504045331),
(x^12 - 2*x^7 - x^6 - 1, 1.2469917463868065),
(x^12 - x^11 - x^10 + x^9 + x^8 - 3*x^6 + 2*x^5 + x^4 - 3*x^3 + x^2 + x - 1, 1.2470184339955448),
(x^12 - 2*x^11 + x^10 - x^2 + 1, 1.2470478623827934),
(x^12 - x^10 - x^2 - 2*x - 1, 1.2470478623827936),
(x^12 + x^11 - x^7 - x^6 - x^5 - 2*x^4 - 2*x^3 - 2*x^2 - x - 1, 1.2473984225894796),
(x^12 - 2*x^11 + 3*x^10 - 4*x^9 + 4*x^8 - 5*x^7 + 5*x^6 - 3*x^5 + 2*x^3 - 3*x^2 + 2*x - 1, 1.2479760484711473),
(x^12 - x^8 - x^6 - 2*x^5 + x^4 - 1, 1.2482596744131684),
(x^12 - x^11 + x^10 - x^9 + 2*x^8 - 2*x^7 - x^6 - x^3 - x^2 + x - 1, 1.2486277289080179),
(x^12 - 2*x^11 + x^10 + x^8 - 2*x^7 + 2*x^6 - 4*x^5 + 5*x^4 - 4*x^3 + 3*x^2 - 2*x + 1, 1.248652405265881),
(x^12 - x^10 + x^8 - 2*x^5 - x^4 - x^2 - 1, 1.2486524052658834),
(x^12 + x^11 - x^9 - x^8 - x^6 - 2*x^5 - x^4 + x^3 - x - 1, 1.2487739150039698),
(x^12 - x^10 + x^9 + x^8 - x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - 1, 1.2493144823143894),
(x^12 - x^11 + x^10 - 2*x^9 + x^8 - 2*x^7 + 2*x^6 - 2*x^5 + 3*x^4 - 2*x^3 + x^2 - x + 1, 1.2498515888642558),
(x^12 + x^11 + x^10 - x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^2 - x - 1, 1.2498515888642558),
(x^12 - 3*x^11 + 5*x^10 - 6*x^9 + 5*x^8 - 4*x^7 + 4*x^6 - 4*x^5 + 3*x^4 - 2*x^3 + x^2 - x + 1, 1.2498515888642563),
(x^12 - x^11 + x^10 - x^8 - x^4 - x^2 - x - 1, 1.249851588864259),
(x^12 - x^10 - x^7 - x^5 + x^2 + 1, 1.250654100616717),
(x^12 + 2*x^11 + x^10 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 4*x^3 - 3*x^2 - 2*x - 1, 1.2506541006167196),
(x^12 - 2*x^10 - 2*x^9 + 3*x^8 + 2*x^7 - 3*x^6 - 2*x^5 + 3*x^4 + 2*x^3 - 2*x^2 + 1, 1.2520730858822462),
(x^12 - x^10 - x^9 + 2*x^7 - x^6 - 2*x^5 + x^3 + x^2 - 1, 1.2525891634349218),
(x^12 - x^11 + x^8 - x^7 - x^6 - x^5 + x^4 + x - 1, 1.2526364965764816),
(x^12 - x^11 + x^9 - x^8 - x^7 + x^6 + x^5 - x^4 - x^3 - x - 1, 1.2533245605734835),
(x^12 + x^11 - x^10 - x^9 - 2*x^5 - 2*x^4 + x^3 + x^2 - x - 1, 1.2535474081280769),
(x^12 - x^11 - x^10 + x^9 - 2*x^5 + 2*x^4 + x^3 - x^2 - x + 1, 1.2535474081280786),
(x^12 + x^11 - 3*x^7 - 3*x^6 - x^5 + x + 1, 1.2536304368456626),
(x^12 - x^11 + x^10 - x^9 - x^7 - x^6 + x^5 + x^3 - x^2 + x - 1, 1.2540018246575528),
(x^12 - x^10 + x^9 + x^8 - 3*x^7 - x^6 + x^5 - x^4 - x^3 + x^2 - 1, 1.2540113528372716),
(x^12 - x^11 - x^8 + 2*x^6 - x^4 - x - 1, 1.2540640003365655),
(x^12 + x^9 - x^7 - x^6 - x^5 - 2*x^4 - x^3 - 2*x^2 - 1, 1.2551953496186983),
(x^12 + x^11 - x^10 - x^9 - x^7 - 2*x^6 - x^5 + x^3 + x^2 + x + 1, 1.255195349618699),
(x^12 + 2*x^11 + 2*x^10 + x^9 - x^7 - 3*x^6 - 5*x^5 - 6*x^4 - 5*x^3 - 4*x^2 - 2*x - 1, 1.2551953496186994),
(x^12 + 3*x^11 + 3*x^10 + x^9 - x^7 - 4*x^6 - 7*x^5 - 8*x^4 - 7*x^3 - 5*x^2 - 3*x - 1, 1.2551953496186994),
(x^12 - x^11 - x^10 + x^9 - x^7 + x^5 + x^3 - x^2 + x - 1, 1.2551953496186998),
(x^12 + x^11 + x^10 + x^9 - x^7 - 2*x^6 - 3*x^5 - 4*x^4 - 3*x^3 - 3*x^2 - x - 1, 1.2551953496187023),
(x^12 + x^11 - x^9 - x^8 - x^7 - 2*x^6 - x^5 + x^4 + x^3 - x - 1, 1.2558136703447271),
(x^12 - x^11 - x^9 + x^8 - x^7 + x^5 + x^4 - x^3 - x + 1, 1.25581367034473),
(x^12 - x^10 - x^6 - 2*x^5 + 2*x^3 + x^2 - 1, 1.2558528462225136),
(x^12 - x^8 - x^7 - x^6 + x^5 - x^4 - 1, 1.256005329442428),
(x^12 + 2*x^11 + x^10 - 2*x^9 - 3*x^8 - x^7 - x^5 - x^4 - x^2 - 2*x - 1, 1.2562965091159424),
(x^12 - x^10 - 2*x^9 + x^8 + x^7 - x^5 + x^4 - x^2 + 1, 1.2562965091159444),
(x^12 - x^10 + x^8 - 2*x^7 - x^6 + x^4 - x^2 + 1, 1.256998185563202),
(x^12 - x^9 - x^8 - x^4 + x^3 - 1, 1.2573488355182125),
(x^12 - 2*x^8 - x^7 + x^6 - x^5 + 1, 1.2581784051484295),
(x^12 - x^11 + x^10 - x^9 + x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.2595418996582586),
(x^12 - x^11 + x^10 - x^9 - 2*x^7 + x^6 + x^3 - x^2 + x - 1, 1.2607856610457406),
(x^12 + x^10 - 2*x^9 - x^8 - 2*x^7 - x^6 + 2*x^5 + x^4 + 2*x^3 - x^2 - 1, 1.2614741887386525),
(x^12 + x^11 - 2*x^9 - 2*x^8 - x^7 + x^6 + x^5 - x - 1, 1.2616574713757553),
(x^12 - x^8 - 2*x^5 - x^4 - 1, 1.262650279648631),
(x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^5 + 3*x^4 - 2*x^3 + 2*x^2 - 2*x + 1, 1.2626502796486327),
(x^12 + 2*x^11 + x^10 - 2*x^8 - 5*x^7 - 4*x^6 - x^5 + x^2 + 2*x + 1, 1.2631112188121962),
(x^12 - x^11 + x^10 - 2*x^8 + x^7 - x^6 - x^5 + x^2 - x + 1, 1.2631112188121998),
(x^12 - x^10 - 2*x^8 - x^7 + 2*x^6 + x^5 + x^2 - 1, 1.263111218812202),
(x^12 - 2*x^11 + x^10 - 2*x^8 + 3*x^7 - x^5 + x^2 - 2*x + 1, 1.2631112188122022),
(x^12 + x^11 + x^10 - 2*x^8 - 3*x^7 - 3*x^6 - x^5 + x^2 + x + 1, 1.2631112188122022),
(x^12 + x^10 - 2*x^8 - x^7 - 2*x^6 - x^5 + x^2 + 1, 1.263111218812204),
(x^12 + x^11 - 2*x^9 - x^8 - 2*x^5 - x^4 + x + 1, 1.263236296657127),
(x^12 + x^9 - 2*x^8 - x^7 - 3*x^5 + x^3 + 1, 1.2636945650463565),
(x^12 - x^11 + x^10 - 2*x^9 + x^8 - x^6 - x^4 + 2*x^3 - x^2 + x - 1, 1.2639212733128775),
(x^12 - 2*x^11 + 2*x^10 - x^9 - x^8 + 3*x^7 - 2*x^6 - x^5 + x^4 - x^3 - 1, 1.264033278765612),
(x^12 - x^11 + x^10 - x^9 + x^8 - 2*x^7 - x^4 + x^3 - x^2 + x - 1, 1.2642719423546873),
(x^12 + x^10 - x^9 - x^8 - 2*x^7 - x^6 - x^4 + x^3 + x^2 + 1, 1.264514110621497),
(x^12 + x^11 - 2*x^7 - 2*x^6 - 2*x^4 - 2*x^3 - x - 1, 1.26455394859932),
(x^12 - x^11 - 2*x^7 + 2*x^6 - 2*x^4 + 2*x^3 - x + 1, 1.2645539485993242),
(x^12 + x^11 - 3*x^9 - 2*x^8 + 2*x^6 + x^3 - x - 1, 1.2646426627834177),
(x^12 - x^11 - x^9 + x^8 - x^6 + x^4 + x^3 - x - 1, 1.264885053671951),
(x^12 + x^11 - x^9 - 2*x^8 - x^7 - x^5 + x^3 - x - 1, 1.2650841293082782),
(x^12 - x^11 - x^9 + x^7 - x^5 + 2*x^4 - x^3 - x + 1, 1.2650841293082795),
(x^12 - 2*x^10 - x^9 + 2*x^8 - x^6 + 2*x^5 - x^3 - 1, 1.265992324494349),
(x^12 - 2*x^11 + x^10 + x^9 - x^8 - x^6 + 2*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1, 1.2662396265223361),
(x^12 - x^11 - x^10 + x^9 + x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.2667209749143484),
(x^12 - x^11 - x^10 + 2*x^9 - x^8 - x^7 + x^6 - x^5 + x^4 - x^2 + x - 1, 1.266912517317132),
(x^12 - x^10 - x^8 - x^7 + x^6 + x^5 - x^4 + x^2 - 1, 1.2669808343956634),
(x^12 - x^11 - 2*x^8 + 3*x^7 - x^6 - x^5 + 2*x^4 - 2*x^3 + x - 1, 1.2669993188941833),
(x^12 - x^10 + x^9 - 3*x^7 + x^6 + x^5 - 2*x^4 - x^3 + x^2 - 1, 1.2675892278965684),
(x^12 - x^11 + x^7 - 3*x^6 + x^5 + x - 1, 1.2680565397536405),
(x^12 - x^11 + x^9 - 2*x^8 + x^7 - x^6 - x^5 + 2*x^4 - x^3 + x - 1, 1.2681654032856606),
(x^12 - x^11 + x^10 - x^9 - x^8 + x^7 - 3*x^6 + 3*x^5 - x^4 + x^3 - x^2 + x - 1, 1.2683717551756313),
(x^12 - x^10 - x^7 + x^5 - 2*x^4 + x^2 - 1, 1.268487994162529),
(x^12 - x^11 - x^10 + 2*x^9 - x^7 - x^5 - x^2 + x - 1, 1.2689059410193957),
(x^12 - x^11 + x^10 - x^9 - x^6 + x^3 - x^2 - x - 1, 1.269371335380608),
(x^12 - x^11 + x^9 - x^7 - x^6 + x^5 - 2*x^4 + x^3 - 2*x^2 + x - 1, 1.269882171549021),
(x^12 - x^11 + x^10 - x^8 + x^7 - 3*x^6 + x^5 - x^4 - x^2 + x - 1, 1.2701382308024767),
(x^12 - x^8 - 2*x^6 + x^4 - 2*x^3 - 1, 1.2705939367082013),
(x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 - 2*x^6 + 4*x^5 - 3*x^4 + 2*x^2 - 2*x + 1, 1.270593936708203),
(x^12 + x^10 - 2*x^9 + x^8 - 2*x^7 + x^6 - 2*x^5 - x^4 - x^2 - 1, 1.2718760651858498),
(x^12 - 8*x^11 + 27*x^10 - 48*x^9 + 41*x^8 + 8*x^7 - 70*x^6 + 104*x^5 - 97*x^4 + 64*x^3 - 29*x^2 + 8*x - 1, 1.272019651375559),
(x^12 + x^11 - 2*x^10 - x^9 - 3*x^7 + x^5 + 3*x^3 + 2*x^2 + x + 1, 1.272019657971005),
(x^12 - 2*x^10 - 2*x^8 + 4*x^6 + 2*x^4 - 2*x^2 - 1, 1.2720196594307056),
(x^12 - x^10 - 2*x^8 - x^6 + 2*x^4 + 3*x^2 + 1, 1.272019660325574),
(x^12 + x^11 - x^10 - x^9 - 2*x^8 - 3*x^7 - x^6 + x^5 + 2*x^4 + 3*x^3 + 3*x^2 + x + 1, 1.272019663134068),
(x^12 - x^11 - x^10 + x^9 - 2*x^8 + 3*x^7 - x^6 - x^5 + 2*x^4 - 3*x^3 + 3*x^2 - x + 1, 1.2720196641474293),
(x^12 + 2*x^11 - 2*x^9 - 4*x^8 - 6*x^7 - 2*x^6 + 2*x^5 + 4*x^4 + 6*x^3 + 4*x^2 + 2*x + 1, 1.2720196648954905),
(x^12 - x^11 + x^9 - 4*x^8 + 3*x^7 - 2*x^6 - x^5 + 4*x^4 - 3*x^3 + 4*x^2 - x + 1, 1.2720196713810683),
(x^12 + 2*x^11 + x^10 - 2*x^9 - 6*x^8 - 6*x^7 - 3*x^6 + 2*x^5 + 6*x^4 + 6*x^3 + 5*x^2 + 2*x + 1, 1.2720196730631466),
(x^12 - 2*x^11 - 2*x^10 + 6*x^9 - 2*x^8 - 2*x^7 + 4*x^6 - 6*x^5 + 2*x^4 + 2*x^3 - 2*x^2 + 2*x - 1, 1.2720196918940654),
(x^12 - 3*x^11 + 2*x^10 + 3*x^9 - 8*x^8 + 9*x^7 - 4*x^6 - 3*x^5 + 8*x^4 - 9*x^3 + 6*x^2 - 3*x + 1, 1.2720197254618224),
(x^12 - 4*x^11 + 4*x^10 + 4*x^9 - 12*x^8 + 12*x^7 - 6*x^6 - 4*x^5 + 12*x^4 - 12*x^3 + 8*x^2 - 4*x + 1, 1.2720200707042646),
(x^12 - x^11 + x^9 - 2*x^7 - x^3 + x - 1, 1.272171499220264),
(x^12 - x^11 + x^8 - 2*x^7 - x^6 + 2*x^5 - x^4 + x - 1, 1.2722509444112455),
(x^12 + x^11 + 2*x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - 3*x^4 - 2*x^3 - 2*x^2 - x - 1, 1.2730548368172276),
(x^12 + x^11 - x^9 - x^8 - x^7 + x^6 - x^5 - 3*x^4 - x^3 - x - 1, 1.2734211820907637),
(x^12 + x^11 + x^10 - x^9 - 2*x^8 - 2*x^7 - x^6 - x^3 - x^2 - x - 1, 1.2737258144938992),
(x^12 - x^9 - x^7 - x^6 + x^5 - x^3 - 1, 1.2739092675666115),
(x^12 - 2*x^11 + x^10 - x^8 + x^7 + x^5 - x^4 - x^2 + 2*x - 1, 1.2740495198502455),
(x^12 - x^10 - x^8 - x^7 + x^5 + x^4 - x^2 + 1, 1.2740495198502462),
(x^12 + x^10 - x^8 - x^7 - 2*x^6 - x^5 - x^4 - x^2 - 1, 1.2740495198502466),
(x^12 - x^11 + x^10 - x^8 - x^6 - x^4 - x^2 + x - 1, 1.2740495198502475),
(x^12 + 2*x^11 + x^10 - x^8 - 3*x^7 - 4*x^6 - 3*x^5 - x^4 - x^2 - 2*x - 1, 1.274049519850248),
(x^12 + x^11 + x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - x^4 - x^2 - x - 1, 1.2740495198502486),
(x^12 - x^8 - x^7 + x^6 + x^5 - x^4 - 2*x^3 - 2*x^2 - 2*x - 1, 1.2743457940121865),
(x^12 - x^9 - x^8 - x^7 + x^5 + x^4 - x^3 - 1, 1.2744455074652803),
(x^12 - 2*x^11 + x^10 - x^9 + x^8 + x^7 - 2*x^6 + 3*x^5 - 3*x^4 + x^3 - x^2 + 2*x - 1, 1.2747320028963445),
(x^12 + x^11 + x^10 - x^9 - 2*x^8 - 2*x^7 - 2*x^6 + x^3 - x^2 - x - 1, 1.2747320028963496),
(x^12 + x^10 - x^9 - x^8 - x^7 - 2*x^6 + x^5 - x^4 + x^3 - x^2 - 1, 1.2747320028963507),
(x^12 - x^10 - x^9 - x^8 + x^7 + x^5 + x^4 - x^3 - x^2 + 1, 1.274732002896351),
(x^12 - x^11 + x^10 - x^9 - 2*x^6 + 2*x^5 - 2*x^4 + x^3 - x^2 + x - 1, 1.2747320028963522),
(x^12 + 2*x^11 + x^10 - x^9 - 3*x^8 - 3*x^7 - 2*x^6 - x^5 + x^4 + x^3 - x^2 - 2*x - 1, 1.2747320028963538),
(x^12 - x^7 - 2*x^6 - x^5 - 1, 1.2748145998618212),
(x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 3*x^7 + 2*x^6 - x^5 + 2*x^4 - 2*x^3 + 2*x^2 - 2*x + 1, 1.2748145998618232),
(x^12 - 2*x^11 + 2*x^10 - 2*x^9 + x^8 + x^7 - 2*x^6 + x^5 - x^4 + 1, 1.2750786271025785),
(x^12 - x^8 + x^7 - x^5 - x^4 - 2*x^3 - 2*x^2 - 2*x - 1, 1.2750786271025818),
(x^12 - x^8 - x^7 - x^6 - x^5 + x^4 - 1, 1.2753786576636972),
(x^12 - x^10 + x^8 - x^7 - x^5 - x^4 - x^2 - 1, 1.27538669072562),
(x^12 - 2*x^11 + x^10 + x^8 - 3*x^7 + 4*x^6 - 5*x^5 + 5*x^4 - 4*x^3 + 3*x^2 - 2*x + 1, 1.275386690725623),
(x^12 - x^11 - x^10 + 2*x^9 - 2*x^7 + 2*x^6 - 2*x^5 - x^2 + x - 1, 1.2754254114359713),
(x^12 + x^11 - x^9 - x^8 - 2*x^7 - 2*x^6 + x^4 + x^3 - x - 1, 1.2756154391778474),
(x^12 - x^11 - x^9 + x^8 - 2*x^7 + 2*x^6 + x^4 - x^3 - x + 1, 1.2756154391778478),
(x^12 + x^11 - x^10 - 2*x^9 + x^7 - x^5 - 2*x^4 + x^2 - x - 1, 1.2756697653683013),
(x^12 - x^11 - x^10 + 2*x^8 - x^7 - x^5 + 2*x^3 - x^2 - x + 1, 1.2756697653683027),
(x^12 + x^11 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 + x^4 + 2*x^3 - x - 1, 1.2759721492294875),
(x^12 - x^11 + x^10 - x^9 - x^8 - x^6 + 2*x^5 - x^4 + x^3 - x^2 + x - 1, 1.2759863718119344),
(x^12 - x^11 + x^10 - x^9 - x^6 - 2*x^5 + 2*x^4 - x^3 + x^2 - x + 1, 1.276008155417501),
(x^12 + x^11 - 2*x^9 - x^8 - x^7 - x^6 - x^5 + x^4 + x + 1, 1.2761494952463848),
(x^12 - 2*x^7 - 2*x^6 + 1, 1.2762376411827725),
(x^12 - x^10 - x^7 + x^6 - x^5 - x^2 - 1, 1.2762405082142507),
(x^12 - x^9 - x^7 - 2*x^6 + x^5 + x^3 - 1, 1.2762889346849815),
(x^12 - x^11 - x^6 + x - 1, 1.2763513300268818)],

13: [(x^13 - x^7 - x^6 - 1, 1.145506428012675)],

14: [(x^14 + x^13 - x^8 - 2*x^7 - x^6 - x - 1, 1.1455064280126723),
(x^14 - x^13 - x^8 + x^6 - x + 1, 1.145506428012673),
(x^14 - x^13 - x^9 + x^8 - x^6 + x^5 - x + 1, 1.1478833949961844),
(x^14 + x^13 - x^9 - x^8 - x^6 - x^5 - x - 1, 1.1478833949961875),
(x^14 - x^13 - x^10 + x^9 - x^5 + x^4 - x + 1, 1.1531044139368083),
(x^14 + x^13 - x^10 - x^9 - x^5 - x^4 - x - 1, 1.1531044139368103),
(x^14 - 3*x^13 + 4*x^12 - 4*x^11 + 4*x^10 - 4*x^9 + 2*x^8 + 2*x^7 - 4*x^6 + 4*x^5 - 4*x^4 + 4*x^3 - 4*x^2 + 3*x - 1, 1.1543552117091167),
(x^14 - 2*x^13 + 3*x^12 - 3*x^11 + 3*x^10 - 3*x^9 + x^8 + x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 2*x - 1, 1.1543552117091518),
(x^14 + x^13 - 2*x^8 - 2*x^7 - x - 1, 1.154355211709154),
(x^14 - x^13 - 2*x^8 + 2*x^7 - x + 1, 1.1543552117091542),
(x^14 - x^13 + 2*x^12 - 2*x^11 + 2*x^10 - 2*x^9 - 2*x^6 + 2*x^5 - 2*x^4 + 2*x^3 - 2*x^2 + x - 1, 1.154355211709156),
(x^14 + x^12 - x^11 + x^10 - x^9 - x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 - 1, 1.154355211709157),
(x^14 - x^11 - x^10 - x^7 + x^4 + x^3 + 1, 1.158551893096862),
(x^14 - x^11 - x^10 + x^7 - x^4 + x^3 - 1, 1.1585518930968695),
(x^14 - x^13 - x^11 + x^10 - x^4 + x^3 - x + 1, 1.1624234880569666),
(x^14 + x^13 - x^11 - x^10 - x^4 - x^3 - x - 1, 1.1624234880569673),
(x^14 - x^10 - x^9 - x^7 + x^5 + x^4 - 1, 1.1649884603711251),
(x^14 - x^9 - x^8 - x^6 + x^5 - 1, 1.1696915642911438),
(x^14 - 2*x^8 - x^7 + 1, 1.1703977848482694),
(x^14 - 6*x^13 + 14*x^12 - 15*x^11 + 6*x^10 - x^8 + 12*x^7 - 29*x^6 + 28*x^5 - 6*x^4 - 13*x^3 + 14*x^2 - 6*x + 1, 1.1739849969294294),
(x^14 + x^12 - 2*x^11 - 2*x^9 - x^8 - x^6 + 2*x^5 + 2*x^3 + x^2 + 1, 1.1739850137001069),
(x^14 - 3*x^13 + 3*x^12 - x^11 - x^9 + 2*x^8 - 2*x^6 + x^5 - x^3 + 3*x^2 - 3*x + 1, 1.1742908848526108),
(x^14 + 2*x^13 + 2*x^12 + x^11 - x^9 - 3*x^8 - 4*x^7 - 3*x^6 - x^5 - x^3 - 2*x^2 - 2*x - 1, 1.17429088485262),
(x^14 - 2*x^13 + 2*x^12 - x^11 - x^9 + x^8 - x^6 + x^5 - x^3 + 2*x^2 - 2*x + 1, 1.174290884852623),
(x^14 - x^11 - x^9 - x^8 + x^6 + x^5 - x^3 + 1, 1.1742908848526232),
(x^14 + 3*x^13 + 3*x^12 + x^11 - x^9 - 4*x^8 - 6*x^7 - 4*x^6 - x^5 - x^3 - 3*x^2 - 3*x - 1, 1.1742908848526232),
(x^14 + x^11 - x^9 - x^8 - x^6 - x^5 - x^3 - 1, 1.1742908848526235),
(x^14 + x^13 - x^12 - x^11 - x^9 - 2*x^8 + 2*x^6 + x^5 - x^3 - x^2 + x + 1, 1.1742908848526246),
(x^14 - x^13 - x^12 + x^11 - x^9 + 2*x^7 - x^5 - x^3 + x^2 + x - 1, 1.174290884852625),
(x^14 - x^13 + x^12 - x^11 - x^9 + x^5 - x^3 + x^2 - x + 1, 1.1742908848526272),
(x^14 + x^13 + x^12 + x^11 - x^9 - 2*x^8 - 2*x^7 - 2*x^6 - x^5 - x^3 - x^2 - x - 1, 1.1742908848526286),
(x^14 - 2*x^13 + 2*x^12 - 2*x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^5 - 2*x^4 + 2*x^2 - 2*x + 1, 1.1743131589071245),
(x^14 - 2*x^6 - 2*x^3 - 1, 1.1743131589071292),
(x^14 - 2*x^9 - 1, 1.1745178431664562),
(x^14 + x^13 - x^9 - 2*x^8 - 2*x^7 + x^5 - x - 1, 1.1745309682797738),
(x^14 - x^13 - x^9 + 2*x^6 - x^5 - x + 1, 1.1745309682797742),
(x^14 - x^9 - x^8 - x^7 + x^6 - x^5 + 1, 1.175969998971801),
(x^14 + x^13 - x^11 - x^10 - x^7 - 2*x^6 + x^4 + x^3 - x - 1, 1.177954173108655),
(x^14 - x^13 - 2*x^9 + 2*x^8 - x + 1, 1.1781079809007178),
(x^14 + x^13 - 2*x^9 - 2*x^8 - x - 1, 1.1781079809007196),
(x^14 - 3*x^13 + 3*x^12 - x^11 - x^10 + 3*x^9 - 3*x^8 + 3*x^6 - 3*x^5 + x^4 - x^3 + 3*x^2 - 3*x + 1, 1.1783750148550596),
(x^14 - x^13 + x^12 - x^11 - x^10 + x^9 - x^8 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1, 1.1783750148550858),
(x^14 - x^13 - x^12 + x^11 - x^10 + x^9 + x^8 - 2*x^7 + x^6 + x^5 - x^4 - x^3 + x^2 + x - 1, 1.178375014855086),
(x^14 + 3*x^13 + 3*x^12 + x^11 - x^10 - 3*x^9 - 3*x^8 - 2*x^7 - 3*x^6 - 3*x^5 - x^4 - x^3 - 3*x^2 - 3*x - 1, 1.1783750148550867),
(x^14 + x^13 - x^12 - x^11 - x^10 - x^9 + x^8 - x^6 + x^5 + x^4 - x^3 - x^2 + x + 1, 1.1783750148550876),
(x^14 - 2*x^13 + 2*x^12 - x^11 - x^10 + 2*x^9 - 2*x^8 + 2*x^6 - 2*x^5 + x^4 - x^3 + 2*x^2 - 2*x + 1, 1.1783750148550889),
(x^14 + x^11 - x^10 - 2*x^7 - x^4 - x^3 - 1, 1.178375014855089),
(x^14 + 2*x^13 + 2*x^12 + x^11 - x^10 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 - 2*x^5 - x^4 - x^3 - 2*x^2 - 2*x - 1, 1.1783750148550893),
(x^14 + x^13 + x^12 + x^11 - x^10 - x^9 - x^8 - 2*x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - x - 1, 1.1783750148550896),
(x^14 - x^11 - x^10 + x^4 - x^3 + 1, 1.1783750148550909),
(x^14 - 2*x^13 + x^12 - x^11 + 2*x^10 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^6 - 2*x^4 + x^3 + x^2 - 2*x + 1, 1.1785964146669963),
(x^14 - x^12 - x^11 + 2*x^5 - x^3 + x^2 - 1, 1.1785964146670187),
(x^14 + 3*x^13 + 5*x^12 + 5*x^11 + 3*x^10 - 3*x^8 - 6*x^7 - 9*x^6 - 10*x^5 - 9*x^4 - 7*x^3 - 5*x^2 - 3*x - 1, 1.178596414667019),
(x^14 + 4*x^13 + 7*x^12 + 7*x^11 + 4*x^10 - 4*x^8 - 8*x^7 - 12*x^6 - 14*x^5 - 12*x^4 - 9*x^3 - 7*x^2 - 4*x - 1, 1.178596414667019),
(x^14 + 2*x^13 + 3*x^12 + 3*x^11 + 2*x^10 - 2*x^8 - 4*x^7 - 6*x^6 - 6*x^5 - 6*x^4 - 5*x^3 - 3*x^2 - 2*x - 1, 1.1785964146670191),
(x^14 + x^13 + x^12 + x^11 + x^10 - x^8 - 2*x^7 - 3*x^6 - 2*x^5 - 3*x^4 - 3*x^3 - x^2 - x - 1, 1.1785964146670194),
(x^14 - x^13 + x^12 - x^11 + x^10 - 2*x^9 + x^8 - 2*x^7 + x^6 - x^4 + x^3 + x^2 - x + 1, 1.17859641466702),
(x^14 + x^13 + x^12 - x^11 - x^10 - 2*x^9 - x^8 - 2*x^7 - x^6 + x^4 + x^3 + x^2 + x + 1, 1.1785964146670205),
(x^14 + x^12 - x^11 - 2*x^9 - 2*x^7 + x^3 + x^2 + 1, 1.178596414667021),
(x^14 + 2*x^13 + x^12 - x^11 - 2*x^10 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 + 2*x^4 + x^3 + x^2 + 2*x + 1, 1.1785964146670231),
(x^14 - x^13 - x^12 + x^11 - x^3 + x^2 - x + 1, 1.179102948240262),
(x^14 + x^13 - x^12 - x^11 - x^3 - x^2 - x - 1, 1.179102948240263),
(x^14 + 2*x^13 + x^12 - 2*x^9 - 4*x^8 - 2*x^7 - x^2 - 2*x - 1, 1.1822598571242349),
(x^14 - 2*x^13 + x^12 - 2*x^9 + 4*x^8 - 2*x^7 - x^2 + 2*x - 1, 1.1822598571242364),
(x^14 + x^13 + x^12 - 2*x^9 - 2*x^8 - 2*x^7 - x^2 - x - 1, 1.1822598571242378),
(x^14 + x^12 - 2*x^9 - 2*x^7 - x^2 - 1, 1.18225985712424),
(x^14 - x^13 + x^12 - 2*x^9 + 2*x^8 - 2*x^7 - x^2 + x - 1, 1.1822598571242402),
(x^14 - x^12 - 2*x^9 + 2*x^7 - x^2 + 1, 1.1822598571242404),
(x^14 - x^13 + x^12 - x^11 - x^7 + x^3 - x^2 + x - 1, 1.1835340233164475),
(x^14 - x^10 - x^9 + x^5 - x^4 - 1, 1.183699208761436),
(x^14 - x^12 - x^11 + x^9 + x^8 - 2*x^7 - x^6 + x^5 + 2*x^4 - x^3 - x^2 + 1, 1.1846044861349116),
(x^14 + 2*x^13 + x^12 - x^11 - 2*x^10 - x^9 + x^8 - 3*x^6 - 3*x^5 + x^3 - x^2 - 2*x - 1, 1.184604486134912),
(x^14 + x^13 + x^12 - 2*x^10 - 2*x^9 - 2*x^8 - x^7 + x^2 + x + 1, 1.1861349405170576),
(x^14 - x^13 + x^12 - x^11 - x^8 + x^7 - x^6 + x^3 - x^2 + x - 1, 1.1869401852498893),
(x^14 - 3*x^13 + 3*x^12 - x^11 - 2*x^9 + 6*x^8 - 6*x^7 + 2*x^6 - x^3 + 3*x^2 - 3*x + 1, 1.1870930029027726),
(x^14 + x^13 - x^12 - x^11 - 2*x^9 - 2*x^8 + 2*x^7 + 2*x^6 - x^3 - x^2 + x + 1, 1.1870930029028146),
(x^14 - x^13 + x^12 - x^11 - 2*x^9 + 2*x^8 - 2*x^7 + 2*x^6 - x^3 + x^2 - x + 1, 1.1870930029028155),
(x^14 - 2*x^13 + 2*x^12 - x^11 - 2*x^9 + 4*x^8 - 4*x^7 + 2*x^6 - x^3 + 2*x^2 - 2*x + 1, 1.1870930029028162),
(x^14 + x^13 + x^12 + x^11 - 2*x^9 - 2*x^8 - 2*x^7 - 2*x^6 - x^3 - x^2 - x - 1, 1.1870930029028162),
(x^14 - x^11 - 2*x^9 + 2*x^6 - x^3 + 1, 1.1870930029028164),
(x^14 + x^11 - 2*x^9 - 2*x^6 - x^3 - 1, 1.1870930029028164),
(x^14 + 2*x^13 + 2*x^12 + x^11 - 2*x^9 - 4*x^8 - 4*x^7 - 2*x^6 - x^3 - 2*x^2 - 2*x - 1, 1.1870930029028164),
(x^14 + 3*x^13 + 3*x^12 + x^11 - 2*x^9 - 6*x^8 - 6*x^7 - 2*x^6 - x^3 - 3*x^2 - 3*x - 1, 1.187093002902817),
(x^14 - x^13 - x^12 + x^11 - 2*x^9 + 2*x^8 + 2*x^7 - 2*x^6 - x^3 + x^2 + x - 1, 1.1870930029028188),
(x^14 - x^8 - x^7 - x^6 - 1, 1.1874990526977027)],

15: [(x^15 - x^8 - x^7 - 1, 1.124883709386441)],

16: [(x^16 - x^15 + x^14 - x^13 + x^12 - x^11 + x^10 - x^9 - x^8 + x^7 - x^6 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.114330918847059),
(x^16 + x^15 - x^9 - 2*x^8 - x^7 - x - 1, 1.124883709386437),
(x^16 - x^15 - x^9 + x^7 - x + 1, 1.1248837093864392),
(x^16 - x^11 - x^10 - x^8 + x^6 + x^5 - 1, 1.1267797071279628),
(x^16 - 6*x^15 + 14*x^14 - 14*x^13 - x^12 + 20*x^11 - 28*x^10 + 20*x^9 - 2*x^8 - 8*x^7 + 8*x^5 + x^4 - 14*x^3 + 14*x^2 - 6*x + 1, 1.1278384884112818),
(x^16 - 2*x^12 - x^8 + 2*x^4 + 1, 1.1278384900599838),
(x^16 - 8*x^15 + 28*x^14 - 56*x^13 + 69*x^12 - 48*x^11 + 48*x^9 - 70*x^8 + 64*x^7 - 56*x^6 + 64*x^5 - 71*x^4 + 56*x^3 - 28*x^2 + 8*x - 1, 1.1278386794353443),
(x^16 - 2*x^9 - 1, 1.1288530614231511),
(x^16 - x^12 - x^9 - x^7 + x^4 + 1, 1.129687457751344),
(x^16 + 2*x^15 + 2*x^14 + 2*x^13 + x^12 - x^9 - 2*x^8 - 3*x^7 - 4*x^6 - 4*x^5 - 3*x^4 - 2*x^3 - 2*x^2 - 2*x - 1, 1.1296874577513463),
(x^16 - x^15 - 2*x^9 + 2*x^8 - x + 1, 1.1313634795322054),
(x^16 + x^15 - 2*x^9 - 2*x^8 - x - 1, 1.1313634795322058),
(x^16 + x^15 - x^12 - x^11 - x^5 - x^4 - x - 1, 1.1350871349472627),
(x^16 - x^15 - x^12 + x^11 - x^5 + x^4 - x + 1, 1.1350871349472638),
(x^16 - x^15 + x^14 - x^13 + x^12 - x^11 - x^8 + x^5 - x^4 + x^3 - x^2 + x - 1, 1.137728226611204),
(x^16 - x^15 + x^13 - x^12 + x^10 - x^9 - x^8 + x^7 - x^6 + x^4 - x^3 + x - 1, 1.1425912667449034)],

17: [(x^17 - x^9 - x^8 - 1, 1.1093819335598598)],

19: [(x^19 - x^10 - x^9 - 1, 1.097304131450219)]},

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