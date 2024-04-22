# import necessary packages
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import numpy as np

# get image for instrument throughputs
def get_filter_throughput(telescope_name, instrument_name, filter_names):
    """Shows image(s) on streamlit from filter_throughput folder based on telescope, instrument, and filters.
    Parameters
    ----------
        telescope_name (str) : telescope
        instrument_name (str) : instrument
        filter_names (list) : list of filters
    Output
    ------
        instrument (str) : instrument name
        image_name (str) : path to corresponding filter throughput image
        caption (str) : caption for filter throughput image
    """

    # filter throughput -- HST/WFC3/UVIS
    wfc3_uvis_filters = [["F200LP", "F300X", "F350LP", "F475X", "F600LP", "F850LP"],
                         ["F218W", "F225W", "F275W", "F336W", "F390W", "F438W",
                          "F475W", "F555W", "F606W", "F625W", "F775W", "F814W"],
                         ["F390M", "F410M", "FQ422M", "F467M", "F547", "F621M", "F689M", "F763", "F845M"],
                         ["FQ232N", "FQ243N", "F280N", "F343N", "F373N", "F373N", "FQ378N", "FQ387N", "F395N"],
                         ["FQ436N", "FQ437N", "F469N", "F487N", "FQ492N", "F502N", "FQ508N"],
                         ["FQ575N", "FQ619N", "F631N", "FQ634N", "F645N", "F656N", "F657", "F658N"],
                         ["F665N", "FQ672N", "F673N", "FQ674N", "F680N", "FQ727N"],
                         ["FQ750N", "FQ889N", "FQ906N", "FQ924N", "FQ937N", "F953N"]]
    wfc3_uvis = pd.DataFrame(wfc3_uvis_filters)

    # filter throughput -- HST/WFC3/IR
    wfc3_ir_filters = [["F105W", "F125W", "F160W", "F110W", "F140W"],
                       ["F098M", "F127M", "F139M", "F153M"],
                       ["F126N", "F128N", "F130N", "F132N", "F164N", "F167N"]]
    wfc3_ir = pd.DataFrame(wfc3_ir_filters)

    # filter throughput -- HST/ACS/WFC
    acs_wfc_filters = [["F435W", "F475W", "F555W", "F606W", "F625W", "F775W", "F814W", "F850LP"],
                       ["F502N", "F550M", "F658N", "F660N", "F892N"]]
    acs_wfc = pd.DataFrame(acs_wfc_filters)

    # filter throughput -- HST/ACS/HRC
    acs_hrc_filters = [["F220W", "F250W", "F330W", "F435W", "F475W", "F555W",
                        "F606W", "F625W", "F775W", "F814W", "F850LP"],
                        ["F344N", "F502N", "F550M", "F658N", "F660N", "F892N"]]
    acs_hrc = pd.DataFrame(acs_hrc_filters)

    if telescope_name == "JWST":
        instrument = instrument_name.split("/")[0]

        # filter throughput -- JWST
        if instrument == "NIRCAM":
            image_name = "filter_throughput/nircam_filters.png"
            caption = "JWST/NIRCAM filter throughputs."
        elif instrument == "MIRI":
            image_name = "filter_throughput/miri_filters.png"
            caption = "JWST/MIRI filter throughputs."
        elif instrument == "NIRISS":
            image_name = "filter_throughput/niriss_filters.png"
            caption = "JWST/MIRI filter throughputs."

    else:
        instrument = instrument_name
        if instrument == "WFC3/UVIS":
            row_indexes = np.where(np.isin(wfc3_uvis.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/wfc3_uvis_filters" + str(ii) + ".png"
                caption = "HST/WFC3/UVIS filter throughputs."

        if instrument == "WFC3/IR":
            row_indexes = np.where(np.isin(wfc3_ir.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/wfc3_ir_filters" + str(ii) + ".png"
                caption = "HST/WFC3/IR filter throughputs."
        
        if instrument == "ACS/WFC":
            row_indexes = np.where(np.isin(acs_wfc.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/acs_wfc_filters" + str(ii) + ".png"
                caption = "HST/ACS/WFC filter throughputs."
        
        if instrument == "ACS/HRC":
            row_indexes = np.where(np.isin(acs_hrc.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/acs_hrc_filters" + str(ii) + ".png"
                caption = "HST/ACS/HRC filter throughputs."

        if instrument == "ACS/SBC":
            image_name = "filter_throughput/acs_sbc_filters.png"
            caption = "HST/ACS/SBC filter throughputs."
    
    return instrument, image_name, caption

# plot rgb image
def plot_rgb(rgb_data, telescope, instrument, filter_list = None, title = "", filter_colors = None):
    """Plot RGB data.
    Parameters
    ----------
        rgb_data (array) : astropy make_lupton_rgb data
        telescope (str) : telescope name
        instrument (str) : instrument name
        filter_list (list of str) : 
        title (str) : optionally add telescope and plot title
        filter_colors (list of str) : optionally add filters and colors used
    Output
    ------
        fig, ax : figure and axis onto which data is plotted
    """
    # plot image
    fig, ax = plt.subplots(figsize = (10, 10), edgecolor = "black")
    ax.imshow(rgb_data, vmin = 0, vmax = 1, cmap = None, origin = "lower")
    
    # customizations -- get rid of axis labels, plot border, white background
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.gcf().set_facecolor('black')
    ax.patch.set_edgecolor('white')
    ax.patch.set_linewidth(2)

    if title != "": # optionally add telescope name and object name above plot
        plt.subplots_adjust(top = 0.9)
        plt.text(0.01, 1.02, telescope, horizontalalignment = "left", verticalalignment = "bottom", fontsize = 20, c = "white", transform = plt.gca().transAxes)
        plt.text(0.99, 1.02, title, horizontalalignment = "right", verticalalignment = "bottom", fontsize = 20, c = "white", transform = plt.gca().transAxes)

    if filter_colors != None: # optionally add instrument and color-coded filters below plot
        plt.subplots_adjust(bottom = 0.1)

        # add instrument label
        plt.text(0.01, -0.02, instrument + " filters", horizontalalignment = "left", verticalalignment = "top", fontsize = 12, c = "white", transform = plt.gca().transAxes)
        plt.text(0.18, -0.02, '|', horizontalalignment = "left", verticalalignment = "top", fontsize = 12, c = "white", transform = plt.gca().transAxes)

        # add filter label
        for ff in range(len(filter_list)):
            filt = filter_list[ff]
            filt_color = filter_colors[ff]
            xloc = 0.2 + (ff * 0.13)

            plt.text(xloc, -0.02, filt, horizontalalignment = "left", verticalalignment = "top", c = filt_color, fontsize = 12, transform = plt.gca().transAxes)

    return fig, ax
