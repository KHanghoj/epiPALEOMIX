import re
import numpy as np
import matplotlib.pyplot as plt

def unpack(t, gccont, obs, ref):
    return gccont, obs, ref


def readdata(f):
    g, o, r = [],[],[]
    with open(f, 'r') as f_in:
        #next(f_in)  # do not need the header
        gen = (unpack(*re.split(r'\s', line.rstrip('\n')))
                           for line in f_in)
    
        for gccont, obs, ref in gen:
            g.append(float(gccont))
            o.append(float(obs))
            r.append(float(ref))
    return np.asarray(g), np.asarray(o), np.asarray(r)

# def running_mean(x, N):
#     cumsum = np.cumsum(np.insert(x, 0, 0))
#     return (cumsum[N:] - cumsum[:-N]) / N
#def running_mean(x,N):
#    return np.convolve(x, np.ones((N,))/N, mode='same')
def running_mean(x,N):
    score=[]
    for idx in xrange(0,len(x),N):
        score.append(sum(x[idx:idx+N])/N)
    return score
gccontent, obs, ref = readdata('BAM1_GCcorrect_5_finescale')
obs = obs/max(obs)
ref = ref/max(ref)
zeropad = [-0.1,-0.05,1,1.05]
span = 0.2
winsize, tabsize = len(obs)-1, 3

tab_means =  np.asarray(running_mean(obs, tabsize))
ori_point = np.asarray(xrange(0, winsize, tabsize))/float(winsize)
rates = tab_means

pd.DataFrame(ori_point, rates)

# plt.figure()
# plt.plot(gccontent, obs)
# plt.plot(gccontent, running_mean(obs, 3), color='green')
# plt.show()
