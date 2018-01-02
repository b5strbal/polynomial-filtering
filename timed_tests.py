"""
Methods that run the polynomial filtering, measure the running time and print the result on the screen and into files.
"""
import time
from polynomial_filtering import find_polynomial_candidates, get_trace_bound, LOWER, UPPER


def or_type_name(is_orientable, is_orientation_rev):
    """
    Convert the booleans encoding orientations to strings.

    Booleans are more convenient to use in the Cython code, but strings are more intuitive for running tests from the shell.
    """
    if is_orientable:
        if is_orientation_rev:
            return "reversing"
        else:
            return "preserving"
    else:
        return "nonor"

def name_to_bool(or_type_name):
    """
    Convert the strings describing orientations to booleans.
    """
    if or_type_name == "reversing":
        return (True, True)
    elif or_type_name == "preserving":
        return (True, False)
    elif or_type_name == "nonor":
        return (False, False)
    else:
        raise ValueError("The type should be 'reversing', 'preserving' or 'nonor'.")


def timed_test(degree, largest_root_bound, is_orientable, 
                is_orientation_rev, first_trace=None):
    """
    Run and time the polynomial filtering for some surface, and prints the results to a file and to the screen as well.

    EXAMPLES:

        sage: from timed_tests import timed_test
        sage: timed_test(3, 2, False, False)
        Testing degree=3 with limit on dilatation 2
        The surface is nonorientable.
        (x^3 - x^2 - x - 1, 1.839286755214161?)
        Elapsed time:0.0290589332581 seconds.

        sage: timed_test(4, 1.9, True, False)
        Testing degree=4 with limit on dilatation 1.90000000000000
        The surface is orientable. The map is orientation-preserving.
        (x^4 - x^3 - x^2 - x + 1, 1.722083805739043?)
        (x^4 - 2*x^3 + x^2 - 2*x + 1, 1.883203505913526?)
        Elapsed time:0.0429930686951 seconds.

        sage: timed_test(4, 1.9, True, True)
        Testing degree=4 with limit on dilatation 1.90000000000000
        The surface is orientable. The map is orientation-reversing.
        (x^4 - x^3 - x - 1, 1.618033988749895?)
        Elapsed time:0.0113308429718 seconds.

        sage: timed_test(4, 1.9, True, False, 2)
        Testing degree=4 with limit on dilatation 1.90000000000000
        The surface is orientable. The map is orientation-preserving.
        Assuming that the first trace is 2
        (x^4 - 2*x^3 + x^2 - 2*x + 1, 1.883203505913526?)
        Elapsed time:0.0154829025269 seconds.

    """
    start = time.time()
    my_str = "\nTesting degree=" + str(degree) + " with limit on dilatation " +\
                         str(largest_root_bound) + '\n'
    my_str += "The surface is {0}.".format("orientable" if is_orientable else "nonorientable")
    if is_orientable:
        my_str += " The map is orientation-{0}.\n".format("reversing" if is_orientation_rev else "preserving")
    else:
        my_str += "\n"
    if first_trace is not None:
        my_str += "Assuming that the first trace is " + str(first_trace) + "\n"
    result = find_polynomial_candidates(degree,largest_root_bound,is_orientable,
                                        is_orientation_rev,first_trace)

    for item in result:
        my_str += repr(item) + '\n'

    my_str += "Elapsed time:" + str(time.time() - start) + " seconds.\n"
    print my_str

    prefix = or_type_name(is_orientable, is_orientation_rev)

    if first_trace is None:
        myfile = open("{1}_{0}.txt".format(degree, prefix),'w')
    else:
        myfile = open("{2}_{0}_{1}.txt".format(degree, first_trace, prefix),'w')
    myfile.write(my_str)
    myfile.close()


def timed_test_by_first_trace(degree, largest_root_bound, is_orientable, 
                is_orientation_rev):
    """
    Starts several instances of the demo.py script, each testing a different value for the first trace. 

    The results are printed on the screen and to files also.

    EXAMPLE:

        sage: timed_test_by_first_trace(3, 1.84, False, False)

        Testing degree=3 with limit on dilatation 1.84000000000000
        The surface is nonorientable.
        Assuming that the first trace is 0
        Elapsed time:0.0156300067902 seconds.

        Testing degree=3 with limit on dilatation 1.84000000000000
        The surface is nonorientable.
        Assuming that the first trace is 1
        (x^3 - x^2 - x - 1, 1.839286755214161?)
        Elapsed time:0.0665440559387 seconds.

        Testing degree=3 with limit on dilatation 1.84000000000000
        The surface is nonorientable.
        Assuming that the first trace is 2
        Elapsed time:0.0248429775238 seconds.

        Testing degree=3 with limit on dilatation 1.84000000000000
        The surface is nonorientable.
        Assuming that the first trace is 3
        Elapsed time:0.00219202041626 seconds.

    """
    trace_bounds = [get_trace_bound(degree, largest_root_bound, 1,
                                    is_orientable, is_orientation_rev, typ) for
                    typ in [LOWER, UPPER]]

    import os
    traces = sorted(range(trace_bounds[0], trace_bounds[1] + 1), key = lambda x: abs(x))
    for i in traces:
        os.system('sage demo.py ' + or_type_name(is_orientable, is_orientation_rev) + " " + str(degree) + ' --first_trace ' + str(i) + ' &')
