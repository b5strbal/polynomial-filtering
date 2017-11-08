import itertools
x = polygen(ZZ, 'x')

deg_7_candidates = [
    x^7 - x^4 - x^3 - 1
]

deg_8_candidates = [
    x^8 - x^7 + x^6 - x^5 - x^4 + x^3 - x^2 + x - 1, # LEF 32
    x^8 + x^7 - x^5 - 2*x^4 - x^3 - x - 1,
    x^8 - x^7 - x^5 + x^3 - x + 1,
    x^8 + 2*x^7 + x^6 - x^5 - 2*x^4 - 3*x^3 - 3*x^2 - 2*x - 1, # LEF 28
    x^8 - x^6 - x^5 - x^3 + x^2 + 1,  # LEF 24
    x^8 + x^7 - x^6 - x^5 - x^3 - x^2 - x - 1,
    x^8 - x^7 - x^6 + x^5 - x^3 + x^2 - x + 1,
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
    x^8 - x^5 - x^4 - x^3 - 1]

deg_10_candidates = [
    x^10 - x^7 - x^6 - x^5 + x^4 + x^3 - 1,
    x^10 - x^9 + x^8 - x^7 - x^5 + x^3 - x^2 + x - 1,
    x^10 + x^9 - x^6 - 2*x^5 - x^4 - x - 1,
    x^10 - x^9 - x^6 + x^4 - x + 1,
    x^10 - x^9 + x^7 - x^6 - x^5 + x^4 - x^3 + x - 1]

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
    """
    return [[(2*item[0]+2,2*item[1]) for item in partition] for partition in partitions_compact(genus-1)]
    temp = Partitions(genus-1).list()

def orbit_patterns(set_size):
    """
    Return a list whose elements describe the possible orbits of
    ``set_size`` singularities of the same degree. Every element of the
    list is also a list and it describes one scenario. Elements of each list
    are tuples (``orbit_length``, ``num_orbits``). Every orbit length has to be
    even unless it is one. The number of orbits of length 1 has to be even, so
    in that case we return (1, ``num_orbits``/2).
    """
    assert(set_size % 2) == 0
    half_mult = set_size // 2
    temp = [[(2*k,mult) for (k,mult) in partition]
            for partition in partitions_compact(half_mult)]
    # return temp
    final = []
    for partition in temp:
        final.append(partition)
        if partition[0][0] != 2:
            continue
        for i in range(partition[0][1]):
            start = [(1,1*i+1)]  # the true multiplicity is double, but this is
            # convenient, since these come in pairs
            if partition[0][1]-i-1 > 0:
                start.append((2,partition[0][1]-i-1))
            final.append(start + partition[1:])
    return final

def orbit_patterns_with_orders(num_prongs, multiplicity):
    """

    OUTPUT:

    list of list such that each list contains tuples (``orbit_length``,
    ``num_orbits``, ``order``).

    EXAMPLE:

    sage: orbit_patterns_with_orders(6, 2)
    [[(2, 1, 1)], [(2, 1, 3)], [(1, 2, 1)], [(1, 2, 3)]]

    sage: orbit_patterns_with_orders(2, 4)
    [[(4, 1, 1)], [(2, 2, 1)], [(1, 2, 1), (2, 1, 1)], [(1, 4, 1)]]
    """
    result = []
    for pattern in orbit_patterns(multiplicity):
        result.extend(append_options(pattern, (num_prongs//2).divisors()))
    return result

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
    comps = Compositions(num_options+multiplicity, length=num_options).list()
    comps = [[i-1 for i in comp] for comp in comps]
    for comp in comps:
        new_list = list(list_so_far)
        for k in range(num_options):
            mult = comp[k] if data > 1 else 2*comp[k]
            if mult > 0:
                new_list.append((data, mult, options[k]))
        append_options_recursive(current_idx + 1, list_with_mult, options,
                                 new_list, result)

def detailed_orbit_patterns(stratum):
    return list(itertools.product(*[orbit_patterns_with_orders(*x) for x in stratum]))

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

def test(poly):
    rec_poly = reciprocalize(poly)
    lef = lefschetz_numbers(rec_poly,30,True)
    max_lefschetz = max_num_singularities(poly.degree())
    for k in lef:
        if k > max_lefschetz:
            # print k, "is too big Lefschetz number:"
            # print poly
            # break
            return
    # print poly, poly.factor()
    print lef

for poly in deg_8_candidates:
    test(poly)
