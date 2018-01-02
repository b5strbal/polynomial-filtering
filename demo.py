"""
This is script to be run from the command line. The script should be run
in Sage's python, such as 

$ sage demo.py nonor 3
$ sage -python demo.py nonor 3

The script finds all polynomials of the specified degree that can possible occur are the characteristic polynomials of the action on homology of pseudo-Anosov maps (with an orientable invariant foliation) on the closed surface whose first homology has the same dimension as the degree. One can choose to consider orientable or nonorientable surface and on orientable surfaces, one can consider orientation-preserving or -reversing maps.

The results are printed on the screen as well written into a file. The name of the file takes the form [type]_[degree].txt or [type]_[degree]_[first_trace].txt, depending on whether ``first_trace`` is specified or not. See below for the meaning of the three pieces of data. 

USAGE:

sage demo.py [-h] [--largest_root_bound [LARGEST_ROOT_BOUND]]
                    [--separate_by_first_trace] [--first_trace [FIRST_TRACE]]
                    type degree

- ``type`` -- there are three possible choices. 
    "nonor": the surface is nonorientable
    "preserving": the surface is orientable and the map is orientation-preserving
    "reversing": the surface is orientable and the map is orientation-reversing

- ``degree`` -- the degree of polynomials to be searched. When ``type`` is "nonor", this is an integer at least 3. For the other two types, this is an even integer at least 4.

- ``largest_root_bound`` -- polynomials are searched whose Perron root are at most this number. This is optional. If not provided, its value defaults to value or ``upper_bounds`` in `data.sage`. This default value should be close to the minimal stretch factor.

- ``first_trace`` -- also optional. If specified, only polynomials with that first trace are searched. If not specified, all possible first traces are considered.

- ``separate_by_first_trace`` -- if this option is used, then the process of searching for polynomials is split up to multiple jobs depending on the first trace. These separate jobs are invoked using `sage demo.py --first_trace ...`, so as separate processes, they run on separate cores of the processor. Parallelizing computation this way can be really useful on machines with processors with many cores.

EXAMPLES:

Listing degree 3 polynomials on nonorientable surfaces:

>>>>>$ sage demo.py nonor 3
Testing g=3 with limit on dilatation 1.84000000000000
The surface is nonorientable.
(x^3 - x^2 - x - 1, 1.839286755214161?)
Elapsed time:0.0206820964813 seconds.

Specifying the upper bound for the Perron root:

>>>>>$ sage demo.py nonor 3 --largest_root_bound 2
Testing g=3 with limit on dilatation 2.0
The surface is nonorientable.
(x^3 - x^2 - x - 1, 1.839286755214161?)
Elapsed time:0.0174880027771 seconds.

Restricting to examples where the first trace is -1 (no examples found):

>>>>>$ sage demo.py nonor 3 --first_trace -1
Testing g=3 with limit on dilatation 1.84000000000000
The surface is nonorientable.
Assuming that the first trace is -1
Elapsed time:0.00100302696228 seconds.

Parallelizing the computation by the first trace:

>>>>>$ sage demo.py reversing 4 --separate_by_first_trace

This starts 9 processes of the form 

>>>>>$  sage demo.py reversing 4 --first_trace T

where T ranges from -4 to 4.

"""

import argparse
from data import upper_bounds
from timed_tests import name_to_bool, timed_test, timed_test_by_first_trace



parser = argparse.ArgumentParser(description='Test polynomials.')
parser.add_argument('type', type=str)
parser.add_argument('degree', type=int)
parser.add_argument('--largest_root_bound', nargs = '?', type=float)
parser.add_argument('--separate_by_first_trace', action='store_true', default=False)
parser.add_argument('--first_trace', nargs='?', type=int)
args = parser.parse_args()

# print args

is_orientable, is_orientation_rev = name_to_bool(args.type)

if args.degree is not None:
    largest_root_bound = args.largest_root_bound
    if largest_root_bound is None:
        largest_root_bound = upper_bounds[args.type][args.degree]

    if args.separate_by_first_trace:
        timed_test_by_first_trace(args.degree, largest_root_bound,
        is_orientable, is_orientation_rev)
    else:
        timed_test(args.degree, largest_root_bound, is_orientable, is_orientation_rev, args.first_trace)
