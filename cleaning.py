# import necessary packages
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np
import pandas as pd
import math
import streamlit as st

# create single images of filters
def plot_single(data, vmin = 0, vmax = 1, grid = None, wcs_header = None):
    """Plot image data.
    Parameters
        data (2D array) : data to plot
        vmin (float) : minimum normalization value, default set to 0
        vmax (float) : maximum normalization value, default set to 1
        grid (boolean) : optionally add grid
        wcs_header (astropy WCS object) : optionally plot in terms of RA/DEC coordinates
    Returns
        fig, ax : figure and axis onto which data is plotted
    """
    # plot
    if wcs_header == None: # without wcs coords
        fig, ax = plt.subplots(figsize = (10, 10))

        if grid: # add grid
            ax.set(xlim = (0, np.shape(data)[1]), ylim = (0, np.shape(data)[0]))
            ax.locator_params(axis = "x", nbins = 5)
            ax.locator_params(axis = "y", nbins = 5)
            ax.tick_params(which = "major", axis = "both", size = 8, labelsize = 14, color = "black")
            ax.grid(which = "major", color = "gray", alpha = 0.75, linestyle = "solid")
        
    else: # with wcs coords
        fig, ax = plt.subplots(figsize = (10, 10), subplot_kw = {'projection': wcs_header})
        
        # axis labels 
        ra = ax.coords["ra"]
        dec = ax.coords["dec"]
        ra.set_ticks(number = 10)
        dec.set_ticks(number = 10)
        ra.set_major_formatter('hh:mm:ss.s')
        dec.set_major_formatter('dd:mm:ss.s')
        ax.tick_params(which = "major", axis = "both", size = 8, labelsize = 14, color = "black")

        ra.set_axislabel("Right Ascension [hms]", size = 15)
        dec.set_axislabel("Declination [degrees]", size = 15)

        if grid: # add grid
            ax.coords.grid(color = "gray", alpha = 0.75, linestyle = "solid")
    
    # plot data 
    ax.imshow(data, vmin = vmin, vmax = vmax, cmap = "Greys", origin = "lower")
    plt.tight_layout()
    
    return fig, ax

# rescale and normalize
def normalize(data, vmin, vmax):
    """
    Normalize data given a minimum and maximum value to values between 0 and 1.
    Parameters
    ----------
        data (2D array) : data to rescale
        vmin (float) : minimum data value for normalization
        vmax (float) : maximum data value for normalization
    Returns
    -------
        data_normalize (2D array) : normalized array
    """
    # clip data below vmin and above vmax
    data_clip = data.clip(min = vmin, max = vmax)
    
    # normalize to values between 0 and 1
    data_normalize = (data_clip - vmin) / (vmax - vmin)
    
    return data_normalize

# show individual filter plots on streamlit
def show_streamlit(filter_data, filter, key_names, wcs_header = None, grid = None):
    """Show plots on streamlit with adjustable min/max values and returns normalized data.
    Parameters
    ----------
        filter_data (2D array) : data to plot
        filter (str) : filter name
        key_names (str) : strings to create unique keys for streamlit input numbers
        wcs_header (astropy WCS object) : optionally plot in terms of RA/DEC coordinates
        grid (boolean) : optionally add grid
    Output
    ------
        norm_data (2D array) : 
    """
    # find upper/lower limits for data and initial vmax/vmin values
    not_nan_data = filter_data[~ pd.isnull(filter_data)]

    lower, upper = float(math.floor(np.percentile(not_nan_data, 0.1))), float(math.ceil(np.percentile(not_nan_data, 99.9)))
    lower_init, upper_init = np.percentile(not_nan_data, 1), np.percentile(not_nan_data, 90)

    # show on streamlit
    st.write(f"Adjust bounds for {filter} filter.")

    col1, col2 = st.columns(2)
    with col1:
        fmin = st.number_input(f"Lower (> {lower}):", step = 0.1, min_value = lower, max_value = upper, 
                               value = lower_init, key = key_names[0])
    with col2:
        fmax = st.number_input(f"Upper (< {upper}):", step = 0.1, min_value = lower, max_value = upper, 
                               value = upper_init, key = key_names[1])
    
    assert fmin < fmax, "Minimum value must be less than maximum value."

    # normalize data based on streamlit input and plot
    norm_data = normalize(filter_data, fmin, fmax)

    fig, ax = plot_single(data = norm_data, grid = grid, wcs_header = wcs_header)
    st.pyplot(fig)

    return norm_data