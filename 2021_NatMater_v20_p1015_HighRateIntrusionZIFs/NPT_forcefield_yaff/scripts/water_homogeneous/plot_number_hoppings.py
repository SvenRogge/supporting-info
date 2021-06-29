import numpy as np
import matplotlib.pyplot as pt
from molmod.units import pascal
from mpl_toolkits.axes_grid1 import make_axes_locatable

n_hops = np.loadtxt('hopping_events_long.csv',delimiter=',',skiprows=1)

press = np.arange(0,110,10)
loading = n_hops[:,0]
n_hops = n_hops[:,1:]

print(loading)
print(press)

pt.clf()
cmap=pt.cm.get_cmap('hot')

for i in range(len(press)):
    pt.plot(loading[1:], n_hops[1:,i]/loading[1:], 'o',linewidth=0, alpha=0.5, color=cmap(i*0.5/len(press)), label='%i MPa' %press[i])

pt.plot(loading[1:], n_hops[1:,:].mean(axis=1)/loading[1:], 'k', label='mean')
pt.xlim([0,80])
pt.xlabel('Water loading [molecules / unit cell]')
pt.ylabel('Number of hopping events per water molecule')
legend = pt.legend(loc='upper left',fancybox=True, ncol = 3)
pt.savefig('hoppings.pdf', bbox_inches='tight', file='pdf')
pt.savefig('hoppings.svg', bbox_inches='tight', file='svg')

pt.clf()

fig, ax = pt.subplots()
im = ax.imshow(n_hops.T, cmap="RdYlGn_r")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

cbar = ax.figure.colorbar(im, cax=cax)
cbar.ax.set_ylabel('Number of hopping events', rotation=-90, va="bottom",color='0.3')
cbar.ax.spines['bottom'].set_color('0.3')
cbar.ax.spines['top'].set_color('0.3')
cbar.ax.spines['right'].set_color('0.3')
cbar.ax.spines['left'].set_color('0.3')
cbytick_obj = pt.getp(cbar.ax.axes, 'yticklabels')                #tricky
pt.setp(cbytick_obj, color='0.3')


ax.set_xticks(np.arange(len(loading)))
ax.set_yticks(np.arange(len(press)))

ax.set_xticklabels(np.arange(4,80,4))
ax.set_xlim([0,len(loading)-2])
ax.set_yticklabels(press)

ax.set_xlabel('Number of adsorbed water molecules')
ax.set_ylabel('Pressure [MPa]')
ax.invert_yaxis()

ax.tick_params(
    axis='both',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    direction='in',
    length=10,
    colors='0.3') # labels along the left edge are off

ax.spines['bottom'].set_color('0.3')
ax.spines['top'].set_color('0.3')
ax.spines['right'].set_color('0.3')
ax.spines['left'].set_color('0.3')
ax.yaxis.label.set_color('0.3')
ax.xaxis.label.set_color('0.3')

fig.tight_layout()
pt.savefig('hopping_abs_heatmap_long.pdf', bbox_inches='tight', file='pdf')
pt.savefig('hopping_abs_heatmap_long.png', bbox_inches='tight', file='png')

pt.clf()

fig, ax = pt.subplots()
im = ax.imshow((n_hops[1:,:]/loading[1:, None]).T, cmap="RdYlGn_r")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)

cbar = ax.figure.colorbar(im, cax=cax)
cbar.ax.set_ylabel('Number of hopping events per water molecule', rotation=-90, va="bottom",color='0.3')
cbar.ax.spines['bottom'].set_color('0.3')
cbar.ax.spines['top'].set_color('0.3')
cbar.ax.spines['right'].set_color('0.3')
cbar.ax.spines['left'].set_color('0.3')
cbytick_obj = pt.getp(cbar.ax.axes, 'yticklabels')                #tricky
pt.setp(cbytick_obj, color='0.3')


ax.set_xticks(np.arange(len(loading)))
ax.set_yticks(np.arange(len(press)))

ax.set_xticklabels(np.arange(4,80,4))
ax.set_xlim([0,len(loading)-2])
ax.set_yticklabels(press)

ax.set_xlabel('Number of adsorbed water molecules')
ax.set_ylabel('Pressure [MPa]')
ax.invert_yaxis()

ax.tick_params(
    axis='both',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    direction='in',
    length=10,
    colors='0.3') # labels along the left edge are off

ax.spines['bottom'].set_color('0.3')
ax.spines['top'].set_color('0.3')
ax.spines['right'].set_color('0.3')
ax.spines['left'].set_color('0.3')
ax.yaxis.label.set_color('0.3')
ax.xaxis.label.set_color('0.3')

fig.tight_layout()
pt.savefig('hopping_rel_heatmap_long.pdf', bbox_inches='tight', file='pdf')
pt.savefig('hopping_rel_heatmap_long.png', bbox_inches='tight', file='png')
