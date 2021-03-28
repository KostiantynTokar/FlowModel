import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import argparse
import tqdm
import os

parser = argparse.ArgumentParser(description = "Plots scalar field, obstacle and vortex")
parser.add_argument("-i", "--input", type=str, default="input.txt", help="input file name")
parser.add_argument("-o", "--output", type=str, default="output.gif", help="output file name")
parser.add_argument("-t", "--title", type=str, default="", help="figure title")
#group = parser.add_mutually_exclusive_group()
parser.add_argument("--dont_show", action="store_true", help="don't show the plot")
parser.add_argument("--dont_save", action="store_true", help="don't save the plot")
parser.add_argument("--save_steps", type=int, nargs="*", default=[], help="save figure at particular steps")
parser.add_argument("-s", "--speed", type=float, default=1, help="speed coefficient (1 default)")
parser.add_argument("--fps", type=int, default=15, help="frames per second of animation")
parser.add_argument("--dont_repeat", action="store_true", help="don't repeat animation")
parser.add_argument("-e", "--extent", nargs=4, type=float, metavar=('x0', 'y0', 'x1', 'y1'), help="ploating area on rectangle (x0,y0), (x1,y1)")
parser.add_argument("--no_colorbar", action="store_true", help="don't show the colorbar")
parser.add_argument("--sep", type=str, default=";", help="separator in the input file")
args = parser.parse_args()

input_fname = args.input
output_fname = args.output
show = not args.dont_show
save = not args.dont_save
repeat = not args.dont_repeat
title = args.title

speed = args.speed
fps = args.fps
show_colorbar = not args.no_colorbar
sep = args.sep

save_steps = args.save_steps
def fig_step_name(step):
    return os.path.splitext(output_fname)[0] + "_" + str(step) + ".png"

n = 0
x0, y0, x1, y1 = 0, 0, 0, 0
rows = 0
cols = 0
taus = []
matrices = []
n_parts = 0
parts_x = []
parts_y = []
plume_x = []
plume_y = []
vmax = 1
vmin = -1
for k,line in enumerate(open(input_fname)):
    if k == 0:
        n = int(line.split(sep)[0])
    elif k == 1:
        splitted = line.split(sep)[:4]
        x0 = float(splitted[0])
        y0 = float(splitted[1])
        x1 = float(splitted[2])
        y1 = float(splitted[3])
    elif k == 2:
        splited = line.split(sep)[:2]
        rows = int(splited[0])
        cols = int(splited[1])
    elif k == 3:
        splitted = line.split(sep)[:n]
        taus.extend([float(tau_str) for tau_str in splitted])
    elif 4 <= k < n + 4:
        splitted = line.split(sep)[0:rows*cols]
        flatten = np.array([float(val_str) for val_str in splitted])
        #if k > 4 + 10 * n / 100:
        vmin = min(vmin, np.amin(flatten))
        vmax = max(vmax, np.amax(flatten))
        matrices.append(np.reshape(flatten, (rows, cols)))
    elif k == n + 4:
        n_parts = int(line.split(sep)[0])
    elif n + 5 <= k < n + n_parts + 5:
        splitted = line.split(sep)[:-1]
        parts_x.append([float(x) for x in splitted[0::2]])
        parts_y.append([float(y) for y in splitted[1::2]])
    elif n + n_parts + 5 <= k < 2 * n + n_parts + 5:
        splitted = line.split(sep)[0:-1]
        plume_x.append([float(x) for x in splitted[0::2]])
        plume_y.append([float(y) for y in splitted[1::2]])
    else:
        break

if args.extent:
    extent = args.extent
else:
    extent = [x0, y0, x1, y1]

ms_per_frame = 1000 / fps
remainder = 0
frames = []
for k,tau in enumerate(taus):
    remainder += 1000 * tau / speed
    num_frames = int(remainder // ms_per_frame)
    frames.extend([k] * num_frames)
    remainder -= num_frames * ms_per_frame

#fig, ax = plt.subplots(figsize=(10,6), dpi=103)
fig = plt.figure(figsize=(10,6), dpi=103)
gs = fig.add_gridspec(1,1)
#ax_author = fig.add_subplot(gs[1,0])
ax = fig.add_subplot(gs[0,0])
fig.suptitle(title)
fig.text(0.05, 0.05, "Токар К.С.\nII курс маг., ОМ\n2019", fontsize = 10, 
        verticalalignment = 'bottom', horizontalalignment = "left",
        bbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5))

class AnimationProvider:
    def __init__(self, matrices, x0, y0, x1, y1, parts_x, parts_y, plume_x, plume_y, vmin, vmax, extent, save_steps, show_colorbar):
        self.matrices = matrices
        self.plume_x = plume_x
        self.plume_y = plume_y
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.parts_x = parts_x
        self.parts_y = parts_y
        self.vmin = vmin
        self.vmax = vmax
        self.extent = extent
        self.save_steps = save_steps
        self.show_colorbar = show_colorbar
        self.last_saved_step_ind = -1
        
        self.colorbar = None
        self.last_frame = -1
        xi0 = -vmin / (vmax - vmin)
        cdict = {"red" : ((0,0,0),
                          (xi0,1,1),
                          (1,1,1)),
                 "green" : ((0,0,0),
                            (xi0,1,1),
                            (1,0,0)),
                 "blue" : ((0,1,1),
                           (xi0,1,1),
                           (1,0,0))}
        self.cmap = colors.LinearSegmentedColormap("b_r", cdict)
        #self.cmap = colors.LinearSegmentedColormap.from_list("b_r", ["b", "w", "r"], 100)

    def init_func(self):
        pass

    def update_func(self, i, fig, ax):
        if self.last_frame == i:
            return
        last_saved_ind = self.last_saved_step_ind
        for save_ind in self.save_steps[last_saved_ind + 1:]:
            if save_ind >= i:
                break
            self._redraw(save_ind, fig, ax)
            fig.savefig(fig_step_name(save_ind))
            self.last_saved_step_ind += 1
        self.last_frame = i
        self._redraw(i, fig, ax)
        if(len(self.save_steps) > self.last_saved_step_ind + 1 and self.save_steps[self.last_saved_step_ind + 1] == i):
            fig.savefig(fig_step_name(save_ind))
            self.last_saved_step_ind += 1

    def _redraw(self, i, fig, ax):
        ax.cla()
        ax.set_xlim(extent[0], extent[2])
        ax.set_ylim(extent[1], extent[3])
        im = ax.imshow(self.matrices[i], origin='lower', extent=(self.x0,self.x1,self.y0,self.y1), cmap = self.cmap, vmin = self.vmin, vmax = self.vmax)
        for part_x, part_y in zip(self.parts_x, self.parts_y):
            ax.plot(part_x, part_y, color='black', lw=3)
        ax.scatter(self.plume_x[i], self.plume_y[i], c='g', s=0.1)
        ax.set_title("step {}".format(i))
        if not self.show_colorbar:
            return
        if self.colorbar is not None:
            self.colorbar.update_normal(im)
            self.colorbar.set_ticks([self.vmin, *self.colorbar.get_ticks(), self.vmax])
        else:
            self.colorbar = fig.colorbar(im)
            self.colorbar.set_ticks([self.vmin, *self.colorbar.get_ticks(), self.vmax])

pr = AnimationProvider(matrices, x0, y0, x1, y1, parts_x, parts_y, plume_x, plume_y, vmin, vmax, extent, save_steps, show_colorbar)
anim = animation.FuncAnimation(fig, pr.update_func, fargs = (fig, ax), init_func = pr.init_func, frames=tqdm.tqdm(frames), interval=ms_per_frame, repeat = repeat)
if pr.last_saved_step_ind < len(save_steps) - 1:
    pr.update_func(save_steps[-1], fig, ax)
if save:
    anim.save(output_fname, writer = animation.PillowWriter(fps=fps))
if show:
    plt.show()

