#       drv_pylab.py
#       
#       Copyright 2009 charley <charley@hosts-137-205-164-145.phys.warwick.ac.uk>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

# This file contain a number of drivers classes to matplotlib.

from matplotlib.pyplot import figure
from matplotlib import cm
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
from matplotlib import rcParams

from numpy import linspace
from numpy import logspace
from numpy import sign
from numpy import average
from numpy import dtype
from numpy import real
from numpy import imag
from numpy import abs
from numpy import angle

from numpy import log10
from numpy import unwrap
from numpy import arctan2
from numpy import imag
from numpy import real

from scipy.signal import freqz

from misc import unitPrefix, prettyunit


def plotcomplexpolar(metaAry, size = (10, 7.5), dpi = 75, grid = True, legend = 0, fontsize = 15):
    """
    metaArray function to do a simple 1D plot of complex array as magnitude and phase angle.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    """
    
    if legend is None:
        legend = 0
    
    axis = metaAry['range']    
    mag = abs(metaAry.data)
    pha = angle(metaAry.data, deg=True)
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    my0 = min(mag)
    my1 = max(mag)
    py0 = min(pha)
    py1 = max(pha)
    xunit = axis['unit'][0]
    myunit = metaAry['unit']
    pyunit = 'Degree'
    
    # Leave 10% margin in the y axis
    mmean = average((my0, my1))
    mreach = abs(my0-my1) / 2 / 0.9
    my0 = sign(my0-mmean) * mreach + mmean
    my1 = sign(my1-mmean) * mreach + mmean
    
    pmean = average((py0, py1))
    preach = abs(py0-py1) / 2 / 0.9
    py0 = sign(py0-pmean) * preach + pmean
    py1 = sign(py1-pmean) * preach + pmean
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    myunit, my0, my1, myscale = prettyunit(myunit, my0, my1)
    pyunit, py0, py1, pyscale = prettyunit(pyunit, py0, py1)
    
    if myscale != 1:
        mag = mag * myscale
        
    if pyscale != 1:
        pha = pha.copy() * pyscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    mylabl = lbl_repr(metaAry['label'], myunit, "Magnitude")
    pylabl = lbl_repr(metaAry['label'], pyunit, "Phase angle")
    
    title = metaAry['name']
    
    fig = figure(figsize=size, dpi = dpi)
    host = SubplotHost(fig, 111)
    
    fig.add_subplot(host)
    par = host.twinx()
    
    if axis['log'][0] is False:
        x = linspace(x0, x1, len(metaAry))
    else:
        raise NotImplemented, "Log axis is not yet implemented."
    
    host.plot(x, mag, 'b-', label=lbl_repr(axis['label'][0], '', "Magnitude"))
    par.plot(x, pha, 'r--', label=lbl_repr(axis['label'][0], '', "Phase"))
    
    host.grid(grid)
    
    host.set_xlabel(xlabl, fontsize=fontsize)
    host.set_ylabel(mylabl, fontsize=fontsize)
    par.set_ylabel(pylabl, fontsize=fontsize)
    
    host.set_xlim([x0, x1])
    host.set_ylim([my0, my1])
    par.set_ylim([py0, py1])
    
    if fontsize is not None:
        host.set_title(title, fontsize=int(fontsize*1.3))
    else:
        host.set_title(title)
        
    if legend >= 0:
        host.legend(loc=legend)
    
    return fig, host, par
    


def plotcomplex(metaAry, size = (10, 7.5), dpi = 75, grid = True, legend = 0, fontsize = 15):
    """
    metaArray function to do a simple 1D plot of complex array as real and imaginary parts.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    """
    
    if legend is None:
        legend = 0
    
    axis = metaAry['range']
    rdata = metaAry.data.real
    idata = metaAry.data.imag
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    ry0 = min(rdata)
    ry1 = max(rdata)
    iy0 = min(idata)
    iy1 = max(idata)
    xunit = axis['unit'][0]
    ryunit = metaAry['unit']
    iyunit = metaAry['unit']
    
    # Leave 10% margin in the y axis
    rmean = average((ry0, ry1))
    rreach = abs(ry0-ry1) / 2 / 0.9
    ry0 = sign(ry0-rmean) * rreach + rmean
    ry1 = sign(ry1-rmean) * rreach + rmean
    
    imean = average((iy0, iy1))
    ireach = abs(iy0-iy1) / 2 / 0.9
    iy0 = sign(iy0-imean) * ireach + imean
    iy1 = sign(iy1-imean) * ireach + imean
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    ryunit, ry0, ry1, ryscale = prettyunit(ryunit, ry0, ry1)
    iyunit, iy0, iy1, iyscale = prettyunit(iyunit, iy0, iy1)
    
    if ryscale != 1:
        rdata = rdata.copy() * ryscale
        
    if iyscale != 1:
        idata = idata.copy() * iyscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    rylabl = lbl_repr(metaAry['label'], ryunit, "Real part")
    iylabl = lbl_repr(metaAry['label'], iyunit, "Imaginary part")
    
    title = metaAry['name']
    
    fig = figure(figsize=size, dpi = dpi)
    host = SubplotHost(fig, 111)
    
    fig.add_subplot(host)
    par = host.twinx()
    
    if axis['log'][0] is False:
        x = linspace(x0, x1, len(metaAry))
    else:
        raise NotImplemented, "Log axis is not yet implemented."
    
    host.plot(x, rdata, 'b-', label=lbl_repr(axis['label'][0], '', "Real"))
    par.plot(x, idata, 'r--', label=lbl_repr(axis['label'][0], '', "Imaginary"))
    
    host.grid(grid)
    
    host.set_xlabel(xlabl, fontsize=fontsize)
    host.set_ylabel(rylabl, fontsize=fontsize)
    par.set_ylabel(iylabl, fontsize=fontsize)
    
    host.set_xlim([x0, x1])
    host.set_ylim([ry0, ry1])
    par.set_ylim([iy0, iy1])
    
    if fontsize is not None:
        host.set_title(title, fontsize=int(fontsize*1.3))
    else:
        host.set_title(title)
        
    if legend >= 0:
        host.legend(loc=legend)
    
    return fig, host, par
    

def plot1d(metaAry, size = (10, 7.5), dpi = 75, grid = True, legend = None, fontsize = 15,\
            fig = None, ax = None, label = None):
    """
    metaArray function to do a simple 1D plot.
    
    legend:
        'best'  0
        'upper right'   1
        'upper left'    2
        'lower left'    3
        'lower right'   4
        'right'         5
        'center left'   6
        'center right'  7
        'lower center'  8
        'upper center'  9
        'center'        10
    
    label   Label for the legend display, default to metaAry['range']['label'][0]
    
    """
    
    if metaAry.dtype is dtype('complex'):
        return plotcomplex(metaAry, size = size, dpi = dpi, grid = grid, legend = legend, fontsize = fontsize)
        
    if legend is None:
        legend = -1
        
    axis = metaAry['range']
    data = metaAry.data
    
    # Load the plotting ranges and units
    x0 = axis['begin'][0]
    x1 = axis['end'][0]
    y0 = min(metaAry.data)
    y1 = max(metaAry.data)
    xunit = axis['unit'][0]
    yunit = metaAry['unit']
    
    # Leave 10% margin in the y axis
    mean = average((y0, y1))
    reach = abs(y0-y1) / 2 / 0.9
    y0 = sign(y0-mean) * reach + mean
    y1 = sign(y1-mean) * reach + mean
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    yunit, y0, y1, yscale = prettyunit(yunit, y0, y1)
    
    if yscale != 1:
        data = data.copy() * yscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    ylabl = lbl_repr(metaAry['label'], yunit)
    
    title = metaAry['name']
    
    # check if object is 1D metaArray object
    if fig is None:
        fig = figure(figsize=size, dpi = dpi)
    
    if ax is None:
        ax = fig.add_subplot(111)
    else:
        x00, x01 = ax.get_xlim()
        y00, y01 = ax.get_ylim()
        
        x0 = min((x0, x00))
        y0 = min((y0, y00))
        x1 = max((x1, x01))
        y1 = max((y1, y01))
    
    if axis['log'][0] is False:
        x = linspace(x0, x1, len(metaAry))
    else:
        raise NotImplemented
    
    if label is None:
        label = axis['label'][0]
    
    ax.plot(x, data, label=label)
    
    ax.grid(grid)
    
    ax.set_xlabel(xlabl, fontsize=fontsize)
    ax.set_ylabel(ylabl, fontsize=fontsize)
    
    ax.set_xlim([x0, x1])
    ax.set_ylim([y0, y1])
    
    if fontsize is not None:
        ax.set_title(title, fontsize=int(fontsize*1.3))
    else:
        ax.set_title(title)
        
    if legend >= 0:
        ax.legend(loc=legend)
    
    return fig, ax
    


def plot2d(metaAry, size = (10, 7.5), dpi = 75, fontsize = 15, cmap = None, \
            nticks = 5, aspect_ratio = 1.0, corient = 'vertical', cformat = None, 
            vmin = None, vmax = None):
    """
    metaArray function to do a simple 2D plot.
    
    cmap              Colour map, default is pyplot.cm.spectral
    nticks          Number of ticks in the colour bar
    aspect_ratio    Aspect ratio of the plot {float|'ij'|'xy'}
                        float:  Fixed aspect ratio by the given number
                        'ij':   Same aspect ratio as ij space
                        'xy':   Same aspect ratio as xy space 
    corient         Colorbar orientation ('vertical'|'horizontal')
    cformat         Colorbar format [ None | format string | Formatter object ]
    vmin            Minimum value for the colour scale
    vmax            Maximum value for the coloir scale
    """
    
    if cmap is None:
        cmap = cm.spectral
    
    if corient is not 'horizontalt':
        corient = 'vertical'
    
    axis = metaAry['range']
    data = metaAry.data
    
    x0 = axis['begin'][0] 
    x1 = axis['end'][0]
    y0 = axis['begin'][1]
    y1 = axis['end'][1]
    if vmin is None:
        v0 = metaAry.data.min()
    else:
        v0 = vmin
    
    if vmax is None:
        v1 = metaAry.data.max()
    else:
        v1 = vmax
    
    xunit = axis['unit'][0]
    yunit = axis['unit'][1]
    vunit = metaAry['unit']
    
    # Apply unit prefix if unit is defined
    xunit, x0, x1, xscale = prettyunit(xunit, x0, x1)
    yunit, y0, y1, yscale = prettyunit(yunit, y0, y1)
    vunit, v0, v1, vscale = prettyunit(vunit, v0, v1)
    
    if vscale != 1:
        data = data.copy() * vscale
        
    xlabl = lbl_repr(axis['label'][0], xunit)
    ylabl = lbl_repr(axis['label'][1], yunit)
    vlabl = lbl_repr(metaAry['label'], vunit)
    
    ticks = linspace(v0, v1, nticks)
    ticks_lbl = []
    
    for i in range(nticks):
        ticks_lbl.append("%(val)0.4g" % {'val':ticks[i]})
        
    # Aspect ration of the plot
    if aspect_ratio == 'ij':
        ratio = data.shape
        ratio = float(ratio[1]) / ratio[0]
    elif aspect_ratio == 'xy':
        ratio = float(y1 - y0) / float(x1 - x0)
    else:
        try:
            ratio = float(aspect_ratio)
        except:
            print "*** Warning! Unrecognisable aspect ratio spec. Using the default instead."
            ratio = 1.0
    
    ratio /= float(y1 - y0) / float(x1 - x0)
    ratio = abs(ratio)
    
    # Make plot with vertical (default) colorbar
    fig = figure(figsize=size, dpi = dpi)
    ax = fig.add_subplot(111)
    
    extent = (x0, x1, y0, y1)
    cax = ax.imshow(data.transpose()[::-1], cmap=cmap, extent=extent, interpolation = 'bicubic', vmin = v0, vmax = v1, aspect=ratio)
    cbar = fig.colorbar(cax, ticks=ticks, orientation=corient, format=cformat)
    
    # ax.set_size(fontsize)
    ax.set_xlabel(xlabl, fontsize=fontsize)     #   Label font size
    ax.set_ylabel(ylabl, fontsize=fontsize)
    rcParams.update({'font.size': fontsize})    #   Value font size
    
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar.ax.set_yticklabels(ticks_lbl)
    cbar.set_label(vlabl, fontsize=fontsize)
    
    if fontsize is not None:
        ax.set_title(metaAry['name'], fontsize=int(fontsize*1.3))
    else:
        ax.set_title(metaAry['name'])
    
    return fig, ax
    

def lbl_repr(label = None, unit = None, string = None):
    """
    Format axis label and unit into a nice looking string
    
    String: Additional string between label and unit.
    """
    lbl = ''
    try:
        # Ignore label if it is not a string, for it can be None also
        lbl += label
    except TypeError:
        pass
        
    try:
        # Append the additional arguement if exist
        lbl += ' [' + string + ']'
    except TypeError:
        pass
    
    try:
        if unit == '':
            pass                    # Unit less quantities
        else:
            lbl += ' (' + unit + ')'
    except TypeError:
        # Most likely unit is not defined, i.e. not a string.
        lbl += ' (Arb.)'
        
    return lbl
    






