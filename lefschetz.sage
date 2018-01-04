import itertools
import numpy as np
from collections import namedtuple
from data import results

x = polygen(ZZ, 'x')


def partitions_compact(n):
    """
    Return the list of partitions of an integer in a compact form.

    Compact means that instead of describing the partition of 4 as (1, 1, 1, 1),
    it is written as (1, 4), meaning that 1 occurs with multiplicity 4.

    INPUT:

    - ``n`` - a positive integer

    EXAMPLE:

        sage: partitions_compact(4)
        [[(4, 1)], [(1, 1), (3, 1)], [(2, 2)], [(1, 2), (2, 1)], [(1, 4)]]

    """
    temp = Partitions(n).list()
    final = []
    for partition in temp:
        s = set(partition)
        final.append([(k,list(partition).count(k)) for k in s])
    return final

def strata(genus):
    """
    Return the possible strata for the genus g orientable surface that is a lift
    from a nonorientable surface by an unbranched cover.

    The strata are partitions of 4``genus``-4, with 2 added to each element of the partition (the number of prongs is 2 plus the order of the zero). 
    
    The partitions have to have certain symmetries: every number should happen an even number of times, since the involution deck transformation creates a pairing of the singularities.

    INPUT: 
    
    - ``genus`` -- the genus of the surface.

    EXAMPLE:

        sage: strata(5)
        [[(10, 2)], [(4, 2), (8, 2)], [(6, 4)], [(4, 4), (6, 2)], [(4, 8)]]

    """
    return [[(2*item[0]+2,2*item[1]) for item in partition] for partition in partitions_compact(genus-1)]
    temp = Partitions(genus-1).list()


def is_symmetric(partition):
    """
    Decide if a partition is symmetric.

    Symmetric means that every element of the partition is either even or if odd, it occurs an even number of times.

    Lifts of orbits have to satisfy these symmetries, since an orbit of odd length lifts to one orbit of twice the length. The lift of an orbit of even length can have one one or two components, but in either case, the orbit lengths are even.

    INPUT:

    - ``partition`` - a partition (a list) with elements in the compact form (number, multiplicity)

    EXAMPLE:

        sage: load('lefschetz.sage')
        sage: is_symmetric([(3, 4), (2,3)])
        True
        sage: is_symmetric([(3, 4), (3,3)])
        False

    """
    return all(orbit%2 == 0 or mult%2 ==0 for (orbit, mult) in partition)

def orbit_patterns(set_size):
    """
    Return the possible orbit pattern on a set of given size.

    The orbit pattern is simply a partition which is symmetric in the sense of the method ``is_symmetric``.

    INPUT:

    - ``set_size`` -- a positive integer

    EXAMPLE:

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
    Takes a partition an decorate its elements with certain markers in all possible ways.

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
    """
    The recursive method doing the job of ``append_options``.
    """
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
    Return the list of possible orbit combinations of a number of singularities with specified number of prongs.

    INPUT:

    - ``num_prongs`` -- the number of prongs of a singularity
    - ``multiplicity`` -- the multiplicity of singularities with said number of prongs

    EXAMPLE:

    If we have two 6-pronged singularities, either both are fixed or they are interchanged. In each of the cases, on the first return, either all prongs are fixed, or the action of an order 3 rotation (order 6 and 2 are not possible, since the unstable foliation goes to the unstable foliation):

        sage: orbit_combinations_fixed_num_prongs(6, 2)

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
    """
    Return a generator to all the possible orbit combinations in a stratum. (How many orbits are and of what length and what the order of the first return map is for each orbit.)

    INPUT:

    - ``stratum`` -- a list of tuples (num_prongs, multiplicity), encoding a stratum.

    EXAMPLE:

        sage: list(orbit_combinations_in_stratum([(4, 2), (6, 2)]))

        [[Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=2, num_orbits=1, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=1),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=2, num_orbits=1, order_after_first_return=3)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=1)],
        [Orbit(num_prongs=4, orbit_length=1, num_orbits=2, order_after_first_return=2),
        Orbit(num_prongs=6, orbit_length=1, num_orbits=2, order_after_first_return=3)]]

    """
    for item in itertools.product(*[orbit_combinations_fixed_num_prongs(*x) for x in
    stratum]):
        yield list(itertools.chain(*item))

def lefschetz_contribution(orbit, power, lambda_pos=True):
    """
    INPUT:

    - ``orbit`` -- an Orbit
    - ``power`` -- the power of the mapping class to be considered
    - ``lambda_pos`` (default=True) -- whether the stretch factor is positive or negative. Negative means that both (transversely orientable) invariant foliations are reversed.

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


def max_num_singularities(genus):
    """
    Return the maximum number of singularities on the orientable surface with specified genus.

    EXAMPLE:

        sage: load('lefschetz.sage')
        sage: max_num_singularities(3)
        4

    """
    return 2*genus-2

def lefschetz_numbers(poly, max_power=10, is_orientable=True):
    """
    Compute the Lefschetz numbers of powers of a homeomorphism of a surface.

    INPUT:

    - ``poly`` -- the characteristic polynomial of the action on homology
    - ``max_power`` -- the Lefschetz numbers are computed from power 1 to ``max_power`` and organized into a list
    - ``is_orientable`` -- whether the surface is orientable

    EXAMPLE:

        sage: load('lefschetz.sage')
        sage: lefschetz_numbers(x^4-x^3-x^2-x+1)
        [1, -1, -5, -5, -14, -25, -41, -77, -131, -226]
        sage: lefschetz_numbers(x^3-x^2-x-1,is_orientable=False)
        [0, -2, -6, -10, -20, -38, -70, -130, -240, -442]

    """
    mat = companion_matrix(poly)
    base = 2 if is_orientable else 1
    return [base-(mat**k).trace() for k in range(1, max_power+1)]

def reciprocalize(poly):
    """
    Multiply a polynomial by its reciprocal to make it reciprocal.

    EXAMPLE:

        sage: reciprocalize(x^3-x^2-x-1)
        x^6 - x^4 - 4*x^3 - x^2 + 1

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
    """
    Decide which orbit combinations are compatible with a polynomial.

    INPUT:

    - ``poly`` -- a polynomial, the charpoly of the action on homology on a nonorientable surface. So this is polynomials is not reciprocal.
    - ``orbit_combinations`` -- a list of orbit combinations to be checked. The format is the same as the output of ``orbit_combinations_in_stratum``.
    - ``max_power`` -- the maximum power until which consistency is checked. 

    OUTPUT:

    The list of 2-element lists of the form [compatible orbit combination, True/False], where the second element is True is lambda is positive and False otherwise.

    """
    rec_poly = reciprocalize(poly)
    candidate_strata = []

    print "-------------------------------------"
    print "Testing ", poly

    lef_nums = lefschetz_numbers(rec_poly,max_power,True)
    max_lefschetz = max_num_singularities(poly.degree())
    for k in lef_nums:
        if k > max_lefschetz:
            print "-------------------------------------"
            return []

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
            candidate_strata.append([orbits, lambda_pos])
    print "-------------------------------------"
    return candidate_strata


def print_summary(poly, orbits, lambda_pos):
    """
    Print a polynomial, orbit patterns and whether the stretch factor is positive or negative.
    """
    print "Polynomial:", poly
    print "Orbit patterns:", orbits
    print "The stretch factor is: ", "POSITIVE" if lambda_pos else "NEGATIVE"


def get_possible_polynomials(degree, max_power=50):
    """
    Takes the polynomials in data.py and runs the Lefschetz tests on them.

    Since the Lefschetz tests currently work only for nonorientable surfaces, this method as well.

    INPUT:

    - ``degree`` -- the degree of the polynomials
    - ``max_power`` -- the maximum power until which the consistency of Lefschetz numbers are tested

    OUTPUT:

    a dictionary whose keys are polynomials and the values are lists of possible orbit patterns

    EXAMPLE:

        sage: from lefschetz import get_possible_polynomials
        sage: get_possible_polynomials(5)
        -------------------------------------
        Testing  x^5 - x^3 - x^2 - 1
        Polynomial: (x^5 - x^3 - x^2 - 1) * (x^5 + x^3 + x^2 - 1)
        Orbit patterns: [Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)]
        The stretch factor is:  POSITIVE
        4 x 50 dense matrix over Integer Ring
        -------------------------------------

        {x^5 - x^3 - x^2 - 1: [[[Orbit(num_prongs=10, orbit_length=1, num_orbits=2, order_after_first_return=5)],
        True]]}
        

        sage: get_possible_polynomials(10)
        ...
        {x^10 - x^9 - x^6 + x^4 - x + 1: [[[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=1)], True],
        [[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=2)], True],
        [[Orbit(num_prongs=4, orbit_length=9, num_orbits=2, order_after_first_return=2)], True]],
        x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1: [[[Orbit(num_prongs=8, orbit_length=3, num_orbits=2, order_after_first_return=4)], True]]}

    """
    res = {}
    
    for (poly, factor) in results['nonor'][degree]:
        poly = ZZ[x](poly)
        orbit_combinations = [orbits for stratum in strata(poly.degree())
                              for orbits in orbit_combinations_in_stratum(stratum)]
        candidate_strata = get_compatible_strata(poly, orbit_combinations,
                                                 max_power)
        if len(candidate_strata) > 0:
            res[poly] = candidate_strata
    return res


# -------------------------------------------------------------------------
# The polynomials and possible orbit patterns AFTER the Lefschetz tests. 
# -------------------------------------------------------------------------


polynomials_and_orbits = {
    8: {x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1:
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
                  True]]},
    10: {x^10 - x^9 - x^6 + x^4 - x + 1:
                  [[[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=1)],
                   True],
                  [[Orbit(num_prongs=4, orbit_length=18, num_orbits=1, order_after_first_return=2)],
                   True],
                  [[Orbit(num_prongs=4, orbit_length=9, num_orbits=2, order_after_first_return=2)],
                   True]],
                 x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1:
                 [[[Orbit(num_prongs=8, orbit_length=3, num_orbits=2, order_after_first_return=4)],
                   True]]},
    
    12: {
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
}


def double_check_results(degree, max_power=100, debug=False, 
                        print_tables=False):
    """
    Tests that the polynomials and orbit patterns in ``polynomials_and_orbits``
    are indeed consistent.

    INPUT:

    - ``degree`` -- the degree of the polynomials
    - ``max_power`` -- the maximum power until which the consistency of Lefschetz numbers are tested
    - ``debug`` -- if True, info useful for debugging is printed out
    - ``print_tables`` -- if True, the consistent Lefschetz tables are printed out

    EXAMPLE:

    sage: load('lefschetz.sage')
    sage: double_check_results(10)
    ...
    -------------------------------------------
    x^10 - x^9 - x^6 + x^4 - x + 1 [[Orbit(num_prongs=4, orbit_length=9, num_orbits=2, order_after_first_return=2)], True] SEEMS OKAY.
    -------------------------------------------

    """
    poly_dict = polynomials_and_orbits[degree]
    for poly in poly_dict.keys():
        rec_poly = reciprocalize(poly)
        for data in poly_dict[poly]:
            orbits, lambda_pos = data
            try:
                summary = create_summary(rec_poly, orbits, max_power,
                                         lambda_pos, debug)
                print "-------------------------------------------"
                print poly, data, "SEEMS OKAY."
                print "-------------------------------------------"
                if print_tables:
                    print_summary(rec_poly, orbits, lambda_pos)
            except RuntimeError as ex:
                print "-------------------------------------------"
                print poly, data, "HAS FAILED!!!!!!!"
                print "-------------------------------------------"

