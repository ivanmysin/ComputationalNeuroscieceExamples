import numpy as np
import matplotlib.pyplot as plt

u = [np.zeros(100, dtype=np.float64), ]

for i in range(1, 1001):
    u.append(u[-1] + np.random.normal(0, 1, size=100))

print(np.std(u[-1]))
print(np.sqrt(1000))

u = np.stack(u)

plt.plot(u)
plt.show()