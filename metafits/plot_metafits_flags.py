#!/usr/bin/env python

import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

parser = argparse.ArgumentParser(description=
"""Reads MWA metafits.fits files and plots the positions of the tiles.
Tiles without any flagged dipoles are shown as green squares,
which tiles with one or more flagged dipoles are shown in red. 
Plots are created for both XX and YY polarizations""")

parser.add_argument('-i', '--image', metavar='', required=True,
                    help='Path to input metafits file')
args = parser.parse_args()


def plot_metafits_flags(img):
    """Plots positions of MWA tiles with flagged dipoles

    Reads MWA metafits.fits files and plots the positions of the tiles.
    Tiles without any flagged dipoles are shown as green squares,
    which tiles with one or more flagged dipoles are shown in red.
    Plots are created for both XX and YY polarizations
        Args:
            img(str, required):   Path to the matafits file to be plotted
    """

    # Opens Header Data Unit
    hdu = fits.open(img)
    # Extracts N,E array positions from hdu[1]
    north = hdu[1].data['North']
    east = hdu[1].data['East']

    # In metafits files, delays of "32" correspond to flagged dipoles
    # If a tile has one or more flagged dipoles, I consider the tile to be flagged
    delays = hdu[1].data['Delays']
    flags = []
    for i in range(delays.shape[0]):
        if 32 in delays[i]:
            f = 1
        else:
            f = 0
        flags.append(f)
    flag_color = ['#ff4a17' if x == 1 else '#5d5d5d' for x in flags]

    # Plot XX, YY Polarization MWA arrays
    plt.style.use('dark_background')

    # [::2] slices to show only YY polarization
    # [1::2] slices to show only XX polarization

    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(11, 7))

    # XX Pol
    ax1.scatter(east[1::2], north[1::2], marker='s',
                s=36, c=flag_color[1::2], alpha=0.5)
    ax1.set_title('XX Flagging - MWA Array Configuraton')
    ax1.set_aspect('equal')
    ax1.set_xlabel('East ($m$)')
    ax1.set_ylabel('North ($m$)')

    # YY Pol
    ax2.scatter(east[::2], north[::2], marker='s',
                s=36, c=flag_color[::2], alpha=0.5)
    ax2.set_title('YY Flagging - MWA Array Configuraton')
    ax2.set_aspect('equal')
    ax2.set_xlabel('East ($m$)')

    # Legends
    unflagged = mlines.Line2D([], [], color='#5d5d5d', marker='s', alpha=0.5,
                              markersize=6, linestyle="None", label='Unflagged Tile')
    flagged = mlines.Line2D([], [], color='#ff4a17', marker='s', alpha=0.5,
                            markersize=6, linestyle="None", label='Flagged Tile')
    ax1.legend(handles=[flagged, unflagged])
    ax2.legend(handles=[flagged, unflagged])

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_metafits_flags(args.image)
