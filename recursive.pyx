import numpy as np
cdef int cUPPER = 0
cdef int cLOWER = 1
cdef long [:] backward_traces = np.zeros(30, dtype=np.int)
cdef long MAXINT = 1000000

cdef do_traces_pass_backwards(long[:] coeffs, long[:] trace_abs_bounds):
    cdef long [:] backward_coeffs = coeffs[-2::-1]

    # coeffs[-1] is the constant term of the polynomial. When it is negative,
    # then all coefficients should be taken with a negative sign.
    backward_traces[0] = -coeffs[-1]*backward_coeffs[0]
    if abs(backward_traces[0]) > trace_abs_bounds[0]:
        return False
    cdef int i, j
    for i in range(1, len(coeffs)-1):
        # using Newton's formula to figure out the next trace
        # E.g. p_4 = -a_1p_3 - a_2p_2 - a_3p_1 - 4a_4
        backward_traces[i] = -coeffs[-1] * (i+1) * backward_coeffs[i]
        for j in range(i):
            backward_traces[i] -= coeffs[-1] * backward_coeffs[j] * backward_traces[i-1-j]
        if abs(backward_traces[i]) > trace_abs_bounds[i]:
            return False
    return True


cpdef find_polynomials_recursive(int num_coeffs,
                                int is_orientable,
                                int is_orientation_rev,
                                long[:] coeffs,
                                long[:] traces,
                                long[:,:]trace_bounds,
                                int next_idx,
                                list good_coeffs,
                                long first_trace):

    """
    INPUT:

    - ``trace_bounds`` -- a list of 2-tuples, each tuple contaning a lower and
    upper bound for the trace. When the recursion starts, the length of this
    list is the number of middle coefficients (in the nonorientable case) and
    about half of the middle coeffs in the orientable case. As we get deeper in
    the recursion, we keep only the part of the tail that we haven't used yet.

    - ``largest_root_bound`` -- the upper bound for the largest root.

    - ``is_orientable`` -- whether the surface is orientable or not

    - ``next_idx`` -- The index of the coefficient to be figured out next. When 0, it means the first coefficient after the main coefficient (which is always 1). 

    """
    cdef int j
    if next_idx == num_coeffs -1:
        if not is_orientable:
            # Choosing the constant coeff to be +1 or -1 if the surface is nonorientable.
            for j in range(-1,2,2):  # [-1,1],making sure it compiles to fast C code
                coeffs[next_idx] = j
                find_polynomials_recursive(num_coeffs,
                                        is_orientable, 
                                        is_orientation_rev,
                                        coeffs, traces,
                                        trace_bounds, next_idx+1,
                                        good_coeffs, first_trace)
            return
        elif is_orientation_rev and num_coeffs % 2 == 0:
            # If the map is orientation reversing and the degree is divisible by four, then the middle coefficient has to vanish.abs
            coeffs[next_idx] = 0
            find_polynomials_recursive(num_coeffs,
                                    is_orientable, 
                                    is_orientation_rev,
                                    coeffs, traces,
                                    trace_bounds, next_idx+1,
                                    good_coeffs, first_trace)
            return

    # There are no more coefficients to determine.
    if next_idx == num_coeffs:
        # For a nonorientable surface, we check if the traces check out backwards as well. If not, we don't append the polynomial to the list of good polynomials. Otherwise we do.
        if not is_orientable and \
           not do_traces_pass_backwards(coeffs, trace_bounds[cUPPER]):
            return
        good_coeffs.append(list(coeffs))
        return

    cdef long trace_upper_bound = trace_bounds[cUPPER, next_idx]
    cdef long trace_lower_bound = trace_bounds[cLOWER, next_idx]

    # Trying to determine the next coefficient from the previous coefficients
    # and traces using Newton's formula, e.g.:
    # 4a_4 = -a_3p_1 - a_2p_2 - a_1p_3 - p_4
    # Here the coefficients are the a_i and the traces of the powers are p_1,
    # p_2, etc.
    coeffs[next_idx] = 0


    for j in range(next_idx):
        # The right hand side of the formula above, without p_4.
        coeffs[next_idx] -= traces[j] * coeffs[next_idx - 1 - j]

    # print 'b'
    # calculating the largest number less than the upper bound for the
    # trace (current_trace_bounds[1]) such that next_coeff minus this number is
    # divisible by n + 1
    traces[next_idx] = trace_upper_bound - (trace_upper_bound -
                                           coeffs[next_idx]) % (next_idx + 1)

    # Completing Newton's formula above to determine the next coefficient
    coeffs[next_idx] -= traces[next_idx]
    coeffs[next_idx] = coeffs[next_idx] / (next_idx + 1)


    while trace_lower_bound <= traces[next_idx] <= trace_upper_bound:
        if next_idx != 0 or first_trace == MAXINT or traces[next_idx] == first_trace:
            if not is_orientable and \
               2 * next_idx >= num_coeffs-1 and \
               coeffs[next_idx] % 2 != coeffs[num_coeffs-next_idx-2] % 2:
                # In the nonorientable case, if at least half of the coefficients
                # have been determined, then we make sure that the polynomial is
                # reciprocal mod 2.
                pass
            else:
                find_polynomials_recursive(num_coeffs,
                                           is_orientable, 
                                           is_orientation_rev,
                                           coeffs, traces,
                                           trace_bounds, next_idx+1,
                                           good_coeffs, first_trace)
        # Choosing a smaller trace and updating the coefficient
        traces[next_idx] -= next_idx + 1
        coeffs[next_idx] += 1