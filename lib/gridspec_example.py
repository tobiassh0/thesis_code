def annotate_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def outside_ticks(fig):
	for i, ax in enumerate(fig.axes):
		ax.tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)

fig = plt.figure()

fig.suptitle("GirdSpec w/ different subplotpars")

## 0% & 11%
gs1 = GridSpec(2, 5, left=0.05, right=0.49, wspace=0., hspace=0.05)# right spacing creates space for gs2 
# 0%
ax1 = fig.add_subplot(gs1[0, :3])
ax2 = fig.add_subplot(gs1[0, 3:],sharey=ax1)
# 11%
ax3 = fig.add_subplot(gs1[1, :3],sharex=ax1,sharey=ax1)
ax4 = fig.add_subplot(gs1[1, 3:],sharex=ax2,sharey=ax1)

# 1% & 50%
gs2 = GridSpec(2, 5, left=0.51, right=0.98, wspace=0., hspace=0.05)
# 1%
ax5 = fig.add_subplot(gs2[0, :3],sharex=ax1,sharey=ax1)
ax6 = fig.add_subplot(gs2[0, 3:],sharex=ax2,sharey=ax1)
# 50%
ax7 = fig.add_subplot(gs2[1, :3],sharex=ax1,sharey=ax1)
ax8 = fig.add_subplot(gs2[1, 3:],sharex=ax2,sharey=ax1)


annotate_axes(fig)
outside_ticks(fig)
plt.savefig('gridspec.png')
