import itertools

import matplotlib.pyplot as plt
import matplotlib.gridspec as grid_spec
import numpy as np


def ridge_plot(data, *args, xlabel=False, ylabel=False, labels=False, figname=None, y_scale=1., **kwargs):
    """Plot data in a ridge plot with fixed width and fixed height per ridge.

    Args:
        data (list): a list of n tuples/lists of length 2: (x, y)-pairs; list of n np.ndarrays: (y)
        xlabel (bool or str, optional): x-label placed at the bottom, nothing if False. Defaults to False.
        ylabel (bool or str, optional): y-label placed at all ridge y-axis, nothing if False. Defaults to False.
        labels (bool or list, optional): list of str with the labels of the ridges; must be the same length as `data`. Defaults to False.
        figname (None or str, optional): first arg in plt.figure(); useful for tracking figure-object. Defaults to None.
        y_scale (float, optional): scale of y axis relative to the default. Defaults to 1.

    *args:
        'dots': add the 'o' str to the plt.<plt_type>
        'slalomaxis': numbers on the y axis change between left and right to prevent overlap
        'x_lim_S': limit the x axis based on the smallest x ticks insted of the largest (default)
        'grid': turn on grid
        'squeeze': set spacing between subplots to -.5. This turns ON 'slalomaxis' and OFF 'grid'.

    **kwargs:
        plt_type (str, optional): plt class (loglog, plot, semilogx etc.) Defaults to plot.
        xlim (list or tuple, optional): min and max value along x axis (len: 2)
        ylim (list or tuple, optional): min and max value along y axis (len: 2)
        color (str, optional): override the color loop with a constant color
        lt (str, optional): override the linetype
    """
    if 'plt_type' in kwargs.keys():
        plt_type = kwargs['plt_type']
    else:
        plt_type = 'plot'
    fsize = (4, y_scale * len(data))
    gs = grid_spec.GridSpec(len(data), 1)
    if figname is not None:
        fig = plt.figure(figname, figsize=fsize)
    else:
        fig = plt.figure(figsize=fsize)
    ax_objs = []
    l2 = []
    ls = ['-', '--']
    ls = itertools.cycle(ls)
    c = ['r', 'g', 'b', 'magenta', 'darkorange',
         'chartreuse', 'firebrick', 'yellow', 'royalblue']
    c = itertools.cycle(c)
    if 'xlim' in kwargs.keys():
        x_min, x_max = kwargs['xlim']
    # TODO: only given T_max and dt (optional), calculate time / x axis
    # elif len([a for a in args if not isinstance(args, str)]) > 0:
    #     if len([a for a in args if not isinstance(args, str)]) == 1:
    #         T = [a for a in args if not isinstance(args, str)][0]
    #         dt = 1
    #     elif len([a for a in args if not isinstance(args, str)]) == 2:
    #         T = [a for a in args if not isinstance(args, str)][0]
    #         dt = [a for a in args if not isinstance(args, str)][1]
    elif len(data[0]) != 2:
        print('reset')
        x_min, x_max = 0, len(data[0])
    elif 'x_lim_S' in args:
        x_min, x_max = x_limit([d[0] for d in data], plt_type, False)
    else:
        x_min, x_max = x_limit([d[0] for d in data], plt_type)

    if ylabel:
        ax = fig.add_subplot(111, frame_on=False)
        ax.tick_params(labelcolor='w', axis='both', which='both', zorder=-1,  # labelleft=False,
                       labelbottom=False, top=False, bottom=False, left=False, right=False)

    # Loop through data
    y_min = np.inf
    y_max = - np.inf
    for i, s in enumerate(data):
        col = next(c)
        lnst = next(ls)
        ax_objs.append(fig.add_subplot(gs[i:i + 1, 0:]))
        if i == 0:
            spines = ["bottom"]
        elif i == len(data) - 1:
            spines = ["top"]
        else:
            spines = ["top", "bottom"]
        y_min = s[1].min() if s[1].min() < y_min else y_min
        y_max = s[1].max() if s[1].max() > y_max else y_max

        # Plot data
        p_func = getattr(ax_objs[-1], plt_type)
        line_type = '-o' if 'dots' in args else '-'
        line_type = kwargs['lt'] if 'lt' in kwargs.keys() else line_type
        clr = kwargs['color'] if 'color' in kwargs.keys() else col
        if len(s) == 2:
            l = p_func(s[0], s[1], line_type, color=clr, markersize=1.5)[0]
        else:
            l = p_func(s, line_type, color=clr, markersize=1.5)[0]

        # Append in line-list to create legend
        l2.append(l)
        ax_objs[-1].patch.set_alpha(0)
        # Scale all subplots to the same x axis
        plt.xlim([x_min, x_max])
        if 'ylim' in kwargs.keys():
            plt.ylim(kwargs['ylim'])

        # The length of data is greater than one, fix the plot according to the input args and kwargs.
        if 'blank' in args:
            spine = ['top', 'bottom', 'left', 'right']
            for sp in spine:
                ax_objs[-1].spines[sp].set_visible(False)
            plt.tick_params(axis='both', which='both', bottom=False, left=False,
                            top=False, right=False, labelbottom=False, labelleft=False)
        else:
            if len(data) != 1:
                if 'squeeze' in args:
                    if i % 2:
                        ax_objs[-1].tick_params(axis='y', which='both', left=False,
                                                labelleft=False, labelright=True)
                        ax_objs[-1].spines['left'].set_color('k')
                    else:
                        ax_objs[-1].tick_params(axis='y', which='both', right=False,
                                                labelleft=True, labelright=False)
                        ax_objs[-1].spines['right'].set_color('k')
                elif 'slalomaxis' in args:
                    if i % 2:
                        ax_objs[-1].tick_params(axis='y', which='both',
                                        labelleft=False, labelright=True)
                for sp in spines:
                    ax_objs[-1].spines[sp].set_visible(False)
                if 'squeeze' not in args:
                    ax_objs[-1].spines['left'].set_color(col)
                    ax_objs[-1].spines['right'].set_color(col)
                ax_objs[-1].tick_params(axis='y', which='both', colors=col)
                ax_objs[-1].yaxis.label.set_color(col)
            if ('grid' in args and 'squeeze' not in args) or ('grid' in args and len(data) == 1):
                plt.grid(True, which="major", ls="-", alpha=0.2)
            elif 'grid' in args and 'squeeze' in args:
                plt.minorticks_off()
                alpha = .2 if i in (0, len(data) - 1) else .1
                plt.grid(True, axis='y', which="major", ls=lnst, alpha=0.2)
                plt.grid(True, axis='x', which="major", ls='-', alpha=alpha)
            if i == len(data) - 1:
                if xlabel:
                    plt.xlabel(xlabel)
                if len(data) != 1:
                    plt.tick_params(axis='x', which='both', top=False)
            elif i == 0:
                plt.tick_params(axis='x', which='both',
                                bottom=False, labelbottom=False)  # , labeltop=True
            else:
                plt.tick_params(axis='x', which='both', bottom=False,
                                top=False, labelbottom=False)

    if ylabel:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt_type = 'log' if plt_type in ['semilogy', 'loglog'] else plt_type
        if plt_type != 'plot':
            ax.set_yscale(plt_type)
        ax.set_ylabel(ylabel)
        y_min = 1e-3 if plt_type == 'log' and y_min <= 0 else y_min
        ylim = kwargs['ylim'] if 'ylim' in kwargs.keys() else [y_min, y_max]
        ax.set_ylim(ylim)
        # fig.text(0.01, 0.5, ylabel, ha='left', va='center', rotation='vertical')

    if labels:
        if len(labels) == len(data):
            l_d = len(data)
            c_max = 4
            n_row = int(np.ceil(l_d / c_max))
            n_col = 1
            while l_d > n_col * n_row:
                n_col += 1
            fig.legend(l2, labels, loc='lower center',  bbox_to_anchor=(.5, 1.),
                    bbox_transform=ax_objs[0].transAxes, ncol=n_col)
        else:
            print('Length of labels and data was not equal.')
    if 'squeeze' in args:
        gs.update(hspace=-0.5)
    else:
        gs.update(hspace=0.)


def x_limit(data, plt_type, maxx=True):
    t_min = data[0]
    x_max = data[0][-1]
    for t in data[1:]:
        t_0, t_max = np.min(t), np.max(t)
        if maxx:
            t_min = t if t_0 < t_min[0] else t_min
            # t_max = t if t_1 > t_max[-1] else t_max
            x_max = t_max if t_max > x_max else x_max
        else:
            t_min = t if t[0] > t_min[0] else t_min
            # t_max = t if t[-1] < t_max[-1] else t_max
            x_max = t_max if t_max < x_max else x_max
    diff = .05 * (x_max - t_min[0])
    # x_max = t_max[-1] + diff
    x_max += diff
    if plt_type in ['loglog', 'semilogx']:
        x_min = .8 * t_min[t_min > 0][0] if t_min[0] < diff else t_min[0] - diff
        # if x_min < 0:
        #     x_min = 1e-10
    else:
        x_min = t_min[0] - diff
    return x_min, x_max
