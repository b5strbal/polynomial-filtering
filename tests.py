"""
Tests for the polynomial filtering code. 
"""

from data import results, upper_bounds
from polynomial_filtering import find_polynomial_candidates
from timed_tests import name_to_bool


def test(or_type, degree):
    """
    Check the correctness of the results in data.sage for the specified case. 
    
    The polynomials filtering algorithm is run and the results are compared with previously computed results in data.sage.

    INPUT:

    - ``or_type`` -- either "preserving", "reversing" or "nonor"
    - ``degree`` -- the degree of the polynomials

    EXAMPLE:

        sage: from tests import test
        sage: test('preserving',4)
        Test (type: preserving, degree: 4) is successful.

    """
    res = find_polynomial_candidates(degree, upper_bounds[or_type][degree], *name_to_bool(or_type))
    res2 = results[or_type][degree]
    def fail():
        raise ValueError("Test (type: {0}, degree: {1}) has failed.\n{2}\n{3}\n").format(or_type, degree, res, res2)

    if len(res) != len(res2):
        fail(or_type, degree)
    else:
        for i in range(len(res)):
            pair1 = res[i]
            pair2 = res2[i]
            if pair1[0] != pair2[0] or abs(pair1[1]-pair2[1]) > 0.0001:
                fail(or_type, degree)
    print "Test (type: {0}, degree: {1}) is successful.".format(or_type, degree)

def test_all(fast=True):
    """
    Runs all (or most) of the tests.

    INPUT: 

    - ``fast`` (default=True) -- if True, then only the simpler tests are run that are completed in a few seconds. If False, all tests are run, which can take a very long time.

    EXAMPLES:

        sage: from tests import test_all
        sage: test_all()
        Test (type: nonor, degree: 3) is successful.
        Test (type: nonor, degree: 4) is successful.
        Test (type: nonor, degree: 5) is successful.
        Test (type: nonor, degree: 6) is successful.
        Test (type: nonor, degree: 7) is successful.
        Test (type: nonor, degree: 8) is successful.
        Test (type: preserving, degree: 4) is successful.
        Test (type: preserving, degree: 6) is successful.
        Test (type: preserving, degree: 8) is successful.
        Test (type: reversing, degree: 8) is successful.
        Test (type: reversing, degree: 4) is successful.
        Test (type: reversing, degree: 6) is successful.

    """
    for or_type in ['nonor', 'preserving', 'reversing']:
        for degree in results[or_type].keys():
            if fast and degree > 8:
                continue
            test(or_type, degree)

