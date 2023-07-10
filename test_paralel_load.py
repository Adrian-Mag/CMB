import numpy as np
import time
from concurrent.futures import ThreadPoolExecutor

def fun(x, y):
    return x + y

N = 1000
X = np.linspace(0, 10, N)
Y = np.linspace(0, 10, N)
result = np.zeros((N, N))

# Single-threaded code
start_time_single = time.time()

for index_x, x in enumerate(X):
    for index_y, y in enumerate(Y):
        result[index_y, index_x] = fun(x, y)

end_time_single = time.time()
execution_time_single = end_time_single - start_time_single

# Multithreaded code
start_time_multi = time.time()

def compute(index_x, x):
    for index_y, y in enumerate(Y):
        result[index_y, index_x] = fun(x, y)

num_threads = 2
executor = ThreadPoolExecutor(max_workers=num_threads)
futures = [executor.submit(compute, index_x, x) for index_x, x in enumerate(X)]
for future in futures:
    future.result()

end_time_multi = time.time()
execution_time_multi = end_time_multi - start_time_multi

# Print the execution times
print("Single-threaded: ", execution_time_single, "seconds")
print("Multithreaded: ", execution_time_multi, "seconds")
print("Speedup: ", execution_time_single / execution_time_multi)
