import matplotlib.pyplot as plt
from scipy import signal
lst=[]
fname = 'denisovaTSS15kb.txt'
with open(fname) as f:
    next(f)
    for l in f:
        lst.append(float(l.rstrip().split('\t')[2]))
lst=lst[300:len(lst)]
        
# x , y = signal.periodogram(lst)
# ding = 1/x[1:1010]
# plt.plot(ding, y[1:1010])
# plt.show()
x , y = signal.periodogram(lst)
x , y = signal.periodogram(lst,detrend='linear')
# x , y = signal.periodogram(lst,detrend='constant')  # default
# x , y = signal.periodogram(lst,detrend=False)
print (1/x[y.argmax()], y.max())
plt.plot(1/x[1:], y[1:])
#plt.ylim([-1000, 100000])
plt.show()

#(183.72727272727272, 45.889402503661501)
