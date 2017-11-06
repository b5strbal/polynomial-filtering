import os
tb = (-10, 13)
traces = sorted(range(tb[0], tb[1] + 1), key = lambda x: abs(x))
for i in traces:
    os.system('sage demo.sage 13 1.145507 --first_trace ' + str(i))
