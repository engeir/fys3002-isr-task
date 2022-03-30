import time

import isr_task.isr as isr
import isr_task.isr2 as isr2

t0 = time.perf_counter()
isr2.isr()
t1 = time.perf_counter()
print(t1 - t0)
isr.isr()
t2 = time.perf_counter()
print(t2 - t1)
