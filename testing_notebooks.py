from matplotlib.pyplot import * as plt
%matplotlib inline

x = linspace(0, 3*pi, 500)
plt.plot(x, sin(x**2))
plt.show()

