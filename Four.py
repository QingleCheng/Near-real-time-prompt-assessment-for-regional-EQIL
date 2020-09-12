import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import signal
import matplotlib.font_manager as fm
sampRat = 100
T = 6
times12 = fm.FontProperties(fname='C:/Windows/Fonts/times.ttf', size=12)
times14 = fm.FontProperties(fname='C:/Windows/Fonts/times.ttf', size=14)
arial16 = fm.FontProperties(fname='C:/Windows/Fonts/arial.ttf', size=16)
def elastic_spectrum_fft(type, t, ag, n, zeta, dT, Tmax):
    '''
    Compute the elastic response spectrum using FFT method.

    :param type: str
        - SA:   绝对加速度反应谱
        - SV:   相对速度反应谱
        - SD:   相对位移反应谱
        - PSA:  拟加速度反应谱
        - PSV:  拟速度反应谱

    :param t: ndarray
        time

    :param ag: ndarray
        ground motion acceleration

    :param n: int
        number of points of FFT, better be 2^m and larger than sample length

    :param zeta: float
        damping ratio

    :param dT: float
        interval of period of the response spectrum

    :param Tmax: float
        maximum period considered of the response spectrum, should be integer times of dT

    :return:
        (Tn, Sa)
    '''

    sp = np.fft.fft(-ag, n)
    dt = t[1] - t[0]
    freq = np.fft.fftfreq(n, d=dt)
    cf = 2 * np.pi * freq

    Tn = np.linspace(dT, Tmax, int(Tmax / dT))
    cfn = 2 * np.pi / Tn

    # add initial point
    Tn = np.append(np.array([0]), Tn)
    if type in ['SA', 'PSA']:
        ag_max = np.max(ag)
        S = np.array([ag_max])
    elif type in ['SV', 'SD', 'PSV']:
        S = np.array([0])
    else:
        S = np.array([0])
        print('Error: undefined spectrum type!')
        exit(0)

    for cfn1 in cfn:
        H = 1 / (-cf ** 2 + (1j) * 2 * zeta * cfn1 * cf + cfn1 ** 2)
        U = sp * H
        u = np.fft.ifft(U, n)  # u.real is relative displacement time history

        if type == 'SA':
            rd = u.real
            rv = np.gradient(rd, dt)
            ra = np.gradient(rv, dt)
            ag_zero = np.zeros(ra.shape)
            ag_zero[:ag.size] = ag
            aa = ra + ag_zero
            # aa = ra[:ag.size] + ag
            SA1 = np.max(np.abs(aa))
            S = np.append(S, SA1)
        if type == 'SV':
            rd = u.real
            rv = np.gradient(rd, dt)
            SV1 = np.max(np.abs(rv))
            S = np.append(S, SV1)
            pass
        if type == 'SD':
            SD1 = np.max(np.abs(u.real))
            S = np.append(S, SD1)
        if type == 'PSA':
            SD1 = np.max(np.abs(u.real))
            PSA1 = (cfn1 ** 2) * SD1
            S = np.append(S, PSA1)
        if type == 'PSV':
            SD1 = np.max(np.abs(u.real))
            PSV1 = cfn1 * SD1
            S = np.append(S, PSV1)

    return (Tn, S)


def mkdir(path):
    import os
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path) 
        return True
    else:
        return False
parent_path = os.path.abspath(os.path.dirname(os.getcwd()))
temppath=parent_path.split("\\")
staitonname1=temppath[-2]
temppath1=staitonname1.split("_")

mkdir(".//mo1groundmotion//")
with open("input.txt", 'r', encoding='gbk') as f:
    filename= []
    i = 0
    x=[]
    y=[]
    for line in f.readlines():
        if i >= 1:   #跳0行
            temp = line.split('	')
            temp[-1] = temp[-1].strip('\n')#列表最后一个元素 去换行符
            filename.append(temp[0])
            x.append(temp[1])
            y.append(temp[2])
        i += 1

for i in range(len(filename)):
    print(i)
    filename_temp=".//groundmotion//"+filename[i]+".txt"
    a = np.loadtxt(filename_temp, delimiter="\t", skiprows=2)
    period1, sa1 = elastic_spectrum_fft('SA', a[:, 0], a[:, 1], 32768*4, 0.05, 0.01, 10)
    filename_temp1 = ".//mo1groundmotion//" + filename[i] +".txt" 
    f1 = open(filename_temp1, 'w')
    f1.write(x[i]+"\t"+y[i]+"\n")
    f1.write(str(period1.size)+"\n")
    for j in range(period1.size):
        f1.write(str(period1[j]) + "\t"+str(sa1[j])+"\n")
    f1.close()

print("Finish computing the response spectra.")
runInterstation="InteStation.exe 4 2"
