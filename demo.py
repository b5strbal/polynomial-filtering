import argparse
import time
from polynomial_filtering import find_polynomial_candidates

upper_bounds = {3: 1.84, 4: 1.52, 5: 1.43, 6: 1.422, 7: 1.2885, 8:
                1.3568, 9: 1.2173, 10: 1.22262, 11: 1.1743, 12:1.2764, 13:
                1.145507, 14: 1.1875, 15:1.1249, 16:1.1426, 17:1.10939,
                18:1.20515, 19:1.097305, 21:1.087629, 23:1.079704}


def timed_test_by_first_trace(g, largest_root_bound):
    trace_bounds = [get_trace_bound(g, largest_root_bound, 1,
                                    False, False, typ) for
                    typ in [LOWER, UPPER]]

    import os
    traces = sorted(range(trace_bounds[0], trace_bounds[1] + 1), key = lambda x: abs(x))
    print traces
    for i in traces:
        # print 'sage polynomial_filtering.pyx ' + str(g) + ' --first_trace ' + str(i)
        os.system('sage polynomial_filtering.pyx ' + str(g) + ' --first_trace ' + str(i) + ' &')



def timed_test(g, largest_root_bound, first_trace=None):
    """
    Tests and times a genus, and prints the results to a file and to the screen
    as well. (For nonorientable surfaces.)
    """
    start = time.time()
    my_str = "\nTesting g=" + str(g) + " with limit on dilatation " +\
                         str(largest_root_bound) + '\n'
    if first_trace is not None:
        my_str += "Assuming that the first trace is " + str(first_trace) + "\n"
    result = find_polynomial_candidates(g,largest_root_bound,False,
                                        False,first_trace)

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
