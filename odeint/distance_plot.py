import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta_0.4_final.csv")
file = pd.read_csv("/Users/lorenzo/Downloads/final_beta10000.csv") 
file = pd.read_csv("/Users/lorenzo/Downloads/final_beta1e8.csv") 
# EUROPA
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta0.04.csv")
# file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta0.004.csv")

# MERCURY
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.04.csv")
# file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.008.csv")
# file2 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_mercury_beta0.004.csv")

# MOON
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.04.csv")
# file1 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.004.csv")
# file2 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_beta0.0004.csv")

fig = plt.figure(figsize=(8, 6))
plt.style.use('seaborn')
# plt.plot(file2.t, file2.distance,  color=sns.color_palette()[1], label= r'$\beta=0.004$')
# plt.plot(file1.t, file1.distance,  color=sns.color_palette()[3], label= r'$\beta=0.008$' )
plt.plot(file.t, file.distance,  color=sns.color_palette()[2], label= r'$\beta=10^8$')

# for MOON
# plt.plot(file.t, 0.0025*np.exp(-0.004*3*file.t), color='black', linewidth=1)
# plt.plot(file.t, 0.0025*np.exp(-0.0004*2*file.t), color='black', linewidth=1)
# plt.plot(file.t, 0.0025*np.exp(-0.04*3*file.t), color='black', linewidth=1)

# for MERCURY
# plt.plot(file1.t, 0.22*np.exp(-0.004*3*file1.t), color='black', linewidth=1)
# plt.plot(file.t, 0.22*np.exp(-0.04*3*file.t), color='black', linewidth=1)
# plt.plot(file2.t, 0.22*np.exp(-0.008*3*file2.t), color='black', linewidth=1)

plt.title('Distance', fontsize=18)
plt.xlabel('Time [Years]', fontsize=18)
plt.ylabel('Distance [A.U.]', fontsize=18)
# plt.xlim(0, 10)
# plt.ylim(-0.0001,0.00265 )
plt.legend(frameon=True, facecolor='white', bbox_to_anchor=(0.8, 0.08))
fig.subplots_adjust(top = 0.949, right=0.98)
plt.savefig('/Users/lorenzo/Desktop/Three-Body-Problem/odeint/images/final_distance_beta1e8.pdf', dpi=400)
plt.show()