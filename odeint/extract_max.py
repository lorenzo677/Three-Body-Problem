#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# MOON
# def relative_max(arr, t, beta):
#     max_values = []
#     time = []
#     for i in range(1, len(arr) - 1):
#         threshold = 0.0025 * np.exp(-2 * beta * t[i-1])-0.0004/(i*i)
#         if arr[i] > arr[i - 1] and arr[i] > arr[i + 1] and arr[i] >= threshold:
#             max_values.append(arr[i])
#             time.append(t[i])
#     time=np.array(time)
#     return time, max_values

# MERCURY
# def relative_max(arr, t, beta):
#     max_values = []
#     time = []
#     for i in range(1, len(arr) - 1):
#         threshold = 0.2 * np.exp(-3 * beta * t[i-1])-0.06/(i*i)
#         if arr[i] > arr[i - 1] and arr[i] > arr[i + 1] and arr[i] >= threshold:
#             max_values.append(arr[i])
#             time.append(t[i])
#     time=np.array(time)
#     return time, max_values

# def func(x, p1,p2):
#   return p1*np.exp(p2*x)

# EUROPA
def relative_max(arr, t, beta):
    max_values = []
    time = []
    for i in range(1, len(arr) - 1):
        threshold = 0.0025 * np.exp(-4*beta * t[i-1])-0.005/i
        if arr[i] > arr[i - 1] and arr[i] > arr[i + 1] and arr[i] >= threshold:
            max_values.append(arr[i])
            time.append(t[i])
    time=np.array(time)
    return time, max_values

def func(x, p1,p2):
  return p1*np.exp(p2*x)
# MOON
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.04.csv")
# file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.004.csv")
# file2 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.0004.csv")

# MERCURY
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.04.csv")
# file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.008.csv")
# file2 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.004.csv")

# EUROPA
file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta0.04.csv")
file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta0.004.csv")
#%%
time, max = relative_max(file.distance, file.t, 0.04)
time1, max1 = relative_max(file1.distance, file1.t, 0.08)
# time2, max2 = relative_max(file2.distance, file2.t, 0.004)

#%%

# Here you give the initial parameters for p0 which Python then iterates over
# to find the best fit
# popt, pcov = curve_fit(func,time[:900],max[:900],p0=(0.0025,-2*0.04))
# popt1, pcov1 = curve_fit(func,time1,max1,p0=(0.0025,-2*0.004))
# popt2, pcov2 = curve_fit(func,time2,max2,p0=(0.0025,-2*0.0004))

# popt, pcov = curve_fit(func,time[:900],max[:900],p0=(0.0025,-3*0.04))
popt, pcov = curve_fit(func,time,max,p0=(0.0025,-3*0.04))
popt1, pcov1 = curve_fit(func,time1,max1,p0=(0.0025,-3*0.004))
# popt2, pcov2 = curve_fit(func,time2,max2,p0=(0.2,-3*0.004))

# residuals = max[:900]- func(time[:900], *popt)
# ss_res = np.sum(residuals**2)
# ss_tot = np.sum((max[:900]-np.mean(max[:900]))**2)
# r_squared = 1 - (ss_res / ss_tot)
residuals = max- func(time, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((max-np.mean(max))**2)
r_squared = 1 - (ss_res / ss_tot)

residuals1 = max1- func(time1, *popt1)
ss_res1 = np.sum(residuals1**2)
ss_tot1 = np.sum((max1-np.mean(max1))**2)
r_squared1 = 1 - (ss_res1 / ss_tot1)

# residuals2 = max2- func(time2, *popt2)
# ss_res2 = np.sum(residuals2**2)
# ss_tot2 = np.sum((max2-np.mean(max2))**2)
# r_squared2 = 1 - (ss_res2 / ss_tot2)
#%%
fig = plt.figure(figsize=(8, 6))

plt.style.use('seaborn')
# plt.scatter(time[:8250], max[:8250], label=r'$\beta = 0.04$')
# plt.scatter(time1[:800], max1[:800], label=r'$\beta = 0.004$')
# # plt.scatter(time2[:600], max2[:600], label=r'$\beta = 0.0004$')

# plt.plot(time[:8250], popt[0]*np.exp(popt[1]*time)[:8250], color='black')
# plt.plot(time1[:800], popt1[0]*np.exp(popt1[1]*time1)[:800], color='black')
# plt.plot(time2[:600], popt2[0]*np.exp(popt2[1]*time2)[:600], color='black')
plt.scatter(time, max, label=r'$\beta = 0.04$')
plt.scatter(time1, max1, label=r'$\beta = 0.004$')

plt.plot(time, popt[0]*np.exp(popt[1]*time), color='black')
plt.plot(time1, popt1[0]*np.exp(popt1[1]*time1), color='black')
plt.xlabel('Time [Years]', fontsize=18)
plt.ylabel('Distance [A.U.]', fontsize=18)
plt.legend()
perr = np.sqrt(np.diag(pcov))

print(f'BEST PARAMETERS:\t d_0\t\t\t E_d_0\t\t\t alpha\t\t\t E_alpha\t\t R^2')
print(f'beta=0.04\t {popt[0]}\t {np.sqrt(np.diag(pcov))[0]}\t {popt[1]}\t {np.sqrt(np.diag(pcov))[1]}\t {r_squared}')
print(f'beta=0.004\t {popt1[0]}\t {np.sqrt(np.diag(pcov1))[0]}\t {popt1[1]}\t {np.sqrt(np.diag(pcov1))[1]}\t {r_squared1}')
# print(f'beta=0.0004\t {popt2[0]}\t {np.sqrt(np.diag(pcov2))[0]}\t {popt2[1]}\t {np.sqrt(np.diag(pcov2))[1]}\t {r_squared2}')

fig.subplots_adjust(top = 0.949, right=0.98)
# plt.savefig('images/max_fit_distance_europa.pdf', dpi=400)
plt.show()


# %%
