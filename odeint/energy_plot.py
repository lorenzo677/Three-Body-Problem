#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%
# file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output_europa_beta_0.4_final.csv")
file = pd.read_csv("/Users/lorenzo/Downloads/final_beta10000.csv")
file = pd.read_csv("/Users/lorenzo/Downloads/final_beta1e8.csv")

#%%
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
plt.style.use('seaborn')
plt.plot(file.t, file.energy, color=sns.color_palette()[1],)
plt.title('Total Energy', fontsize = 18)
plt.xlabel('Time [Years]', fontsize = 18)
# plt.grid(True)
plt.ylabel(r'Energy [$M_\oplus AU^2y^{-2}$]', fontsize = 18)
# plt.xlim(-1, 51)
# plt.ylim( -0.00089, -0.00085 )
ax.tick_params(axis='both', labelsize=16)
fig.subplots_adjust(top = 0.933, right=0.99, left=0.151)
plt.savefig('/Users/lorenzo/Desktop/Three-Body-Problem/odeint/images/final_energy_BETA1e8.pdf', dpi=400)
plt.show()
# %%
