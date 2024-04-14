# import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import streamlit as st

# create single images of filters
def plot_single(data, vmin, vmax, figsize = (10, 10), wcs_header = None):
    """Plot image data.
    Parameters
        data (2D array) : data to plot
        figsize (tuple) : figure size, default = (15, 13)
        cmap (str) : matplotlib colormap to use
        scale (float) : scale
        wcs : either None or WCS astropy object from fits header
        grid (boolean) : add optional grid (set to False for no grid)
        **kwargs : to be passed to ax.imshow()
    Returns
        fig, ax : figure and axis onto which data is plotted
    """
    # plot
    if wcs_header == None: # without wcs coords
        fig, ax = plt.subplots(figsize = figsize)

        ax.set_xlabel("x", size = 15)
        ax.set_ylabel("y", size = 15)
        
    else: # with wcs coords
        fig, ax = plt.subplots(figsize = figsize, subplot_kw = {'projection': wcs_header})
        ax.coords.grid(color = 'gray', alpha = 0.75, linestyle = 'solid')

        ### force multiple tick marsk
        
        ax.set_xlabel("Right Ascension [hms]", size = 15)
        ax.set_ylabel("Declination [degrees]", size = 15)
        
    ax.imshow(data, vmin = vmin, vmax = vmax, cmap = "Greys", origin = "lower")
    
    return fig, ax

# rescale and normalize
def normalize(data, vmin, vmax):
    """
    Normalize data given a minimum and maximum value to values between 0 and 1.
    Parameters
    ----------
        data (2D array) : data to rescale
        vmin (float) : minimum data value
        vmax (float) : maximum data value
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
def show_streamlit(filter_data, filter, key_names, wcs_header = None):
    """Show plots on streamlit with adjustable min/max values."""
    # streamlit
    not_nan_data = filter_data[~ pd.isnull(filter_data)]

    lower, upper = float(math.floor(np.percentile(not_nan_data, 0.1))), float(math.ceil(np.percentile(not_nan_data, 99.9)))
    lower_init, upper_init = np.percentile(not_nan_data, 1), np.percentile(not_nan_data, 90)

    st.write(f"Adjust bounds for {filter} filter.")

    col1, col2 = st.columns(2)
    with col1:
        fmin = st.number_input(f"Lower (> {lower}):", step = 0.1, min_value = lower, max_value = upper, 
                               value = lower_init, key = key_names[0])
    with col2:
        fmax = st.number_input(f"Upper (< {upper}):", step = 0.1, min_value = lower, max_value = upper, 
                               value = upper_init, key = key_names[1])
    
    assert fmin < fmax, "Minimum value must be less than maximum value."

    norm_data = normalize(filter_data, fmin, fmax)

    # values = st.slider('Select a range of values', lower, upper, [lower_init, upper_init], step = 0.1)
    fig, ax = plot_single(norm_data, vmin = 0, vmax = 1, wcs_header = wcs_header)
    st.pyplot(fig)

    return norm_data