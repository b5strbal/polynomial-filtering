import os
tb = (-9, 12)
traces = sorted(range(tb[0], tb[1] + 1), key = lambda x: abs(x))
for i in traces:
    os.system('sage demo.sage 12 1.2 --first_trace ' + str(i))
