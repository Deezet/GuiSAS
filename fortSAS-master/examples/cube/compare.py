from modelexp.data import XyeData
import matplotlib.pyplot as plt

old = XyeData()
old.loadFromFile('./oldCubeData.xye')
q_old, I_old, sI_old = old.getData()

new = XyeData()
new.loadFromFile('./newCubeData.xye')
q_new, I_new, sI_new = new.getData()



fig, ax = plt.subplots()
ax.errorbar(q_old, I_old, sI_old, label='old')
ax.errorbar(q_new, I_new, sI_new, label='new')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.show()