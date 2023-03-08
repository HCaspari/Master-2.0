import numpy as np

list = []
print(list)
list.append(0)
print(list)
a = np.array()
list.append([1, 2, 3])
print(list)
list[0] = [2, 3, 4]
print(list)

array = np.array(list)
print(array)