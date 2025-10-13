
import matplotlib.pyplot as plt
import numpy as np

x1,y1 = np.loadtxt('Set_cardinality.csv', unpack=True, delimiter = ',')
# x2,y2 = np.loadtxt('Ibm_Overall_NewEncryption_time_samples.csv', unpack=True, delimiter = ',')
# x3,y3 = np.loadtxt('Our_Overall_NewEncryption_time_samples.csv', unpack=True, delimiter = ',')

fig, ax1 = plt.subplots(1)#, 2)
plt.bar(x1,y1)#,  label='frequency vs key switching steps' ,marker='s')#, linestyle='--', linewidth=2)#, markersize=5)

# ax1.plot(x1, y1,'tab:green',label='Zero')
# ax1.plot(x1, y1, label='frequency vs key switching steps' ,marker='s', linestyle='--', linewidth=2, markersize=5)
# ax1.plot(x2, y2, label='Singh et al.'  ,marker='d', linestyle=':', linewidth=2, markersize=5) #color='red'  ,
# ax1.plot(x3, y3, label='Proposed Scheme' ,marker='o', linestyle='--', linewidth=2, markersize=5)
# ax1.plot(x2, y2,'tab:red',label='One')
# ax1.set_linestyle(self, ls)

ax1.grid(True)
# ax1.legend(loc='upper left')

# ax1.set_title('Zhu Scheme')

# ax2.plot(x1, y1,'tab:green',label='Zero')
# ax2.plot(x2, y2,'tab:red',label='One')
# ax2.grid(True)
# ax2.legend(loc='upper right')
# ax2.set_title('Our Scheme')

# plt.subplots_adjust(left=None, bottom=None, right=None, top=None , wspace=None, hspace=None)

ax1.set(xlabel='set_size', ylabel='Frequency')
# ax2.set(xlabel='Dimension', ylabel='Time')
plt.savefig('Set_size_vs_frequency.jpg',dpi=300)
plt.show()




















# x = np.random.randint(low=1, high=11, size=50)
# y = x + np.random.randint(1, 5, size=x.size)
# data = np.column_stack((x, y))

# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,figsize=(8, 4))

# ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
# ax1.set_title('Scatter: $x$ versus $y$')
# ax1.set_xlabel('$x$')
# ax1.set_ylabel('$y$')

# ax2.hist(data, bins=np.arange(data.min(), data.max()),label=('x', 'y'))
# ax2.legend(loc=(0.65, 0.8))
# ax2.set_title('Frequencies of $x$ and $y$')
# ax2.yaxis.tick_right()

# Some example data to display
# x = np.linspace(0, 2 * np.pi, 400)
# y = np.sin(x ** 2)

# axs[0, 0].plot(x, y)
# axs[0, 0].set_title('Axis [0,0]')

# axs[0, 1].plot(x, y, 'tab:orange')
# axs[0, 1].set_title('Axis [0,1]')

# axs[1, 0].plot(x, -y, 'tab:green')
# axs[1, 0].set_title('Axis [1,0]')

# axs[1, 1].plot(x, -y, 'tab:red')
# axs[1, 1].set_title('Axis [1,1]')

# for ax in axs.flat:

# Hide x labels and tick labels for top plots and y ticks for right plots.
# for ax in axs.flat:
# ax.label_outer()
