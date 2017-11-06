from polynomial_filtering import find_polynomial_candidates
import time
import argparse

parser = argparse.ArgumentParser(description='Test polynomials.')
parser.add_argument('degree', nargs = '?', default = 'all', type=int)
parser.add_argument('largest_root_bound', nargs = '?', type=float)
parser.add_argument('--first_trace', type=int)
args = parser.parse_args()

upper_bounds = {3: 1.84, 4: 1.52, 5: 1.43, 6: 1.422, 7: 1.2885, 8:
                1.3568, 9: 1.2173, 10: 1.22262, 11: 1.1743, 13: 1.145507} 

def timed_test(g, largest_root_bound, first_trace=None):
    start = time.time()
    print "\nTesting g=" + str(g) + " with limit on dilatation " +\
                         str(largest_root_bound)
    if first_trace != None:
        print "The first trace is assumed to be " + str(first_trace)
    print find_polynomial_candidates(g,largest_root_bound,False, first_trace)
    print "Elapsed time:" + str(time.time() - start) + " seconds."


    
if args.degree == 'all':
    for g in range(3,12):
        timed_test(g, upper_bounds[g])
        sys.exit()

g = args.degree

largest_root_bound = args.largest_root_bound
if largest_root_bound is None:
    largest_root_bound = upper_bounds[g]

    
timed_test(g, largest_root_bound, args.first_trace)





