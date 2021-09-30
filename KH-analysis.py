import matplotlib.pyplot as plt
import h5py

# Read in the data
f = h5py.File('analysis_tasks/analysis_tasks_s44/analysis_tasks_s44_p0.h5','r')
y = f['/scales/y/1.0'][:]
t = f['scales']['sim_time'][:]
s_ave = f['tasks']['s profile'][:]
f.close()

s_ave = s_ave[:,0,:] # remove length-one x dimension

fig, axis = plt.subplots(figsize=(10,5))

for i in range(0,21,5):
  plt.plot(s_ave[i,:],y,label='t=%4.2f' %t[i])
  
plt.ylim([-0.5,0.5])
plt.xlim([0,1])
plt.xlabel(r'$\frac{\int \ s dx}{L_x}$',fontsize=24)
plt.ylabel(r'$y$',fontsize=24)
plt.legend(loc='lower right').draw_frame(False)
plt.draw()
plt.waitforbuttonpress()

