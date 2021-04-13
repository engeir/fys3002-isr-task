import time

import isr
import isr2

t0 = time.perf_counter()
isr.isr()
t1 = time.perf_counter()
isr2.isr()
t2 = time.perf_counter()
print(t1 - t0)
print(t2 - t1)
