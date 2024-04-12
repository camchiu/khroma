# import necessary packages
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import numpy as np

# method to create and download a full RGB or RGB + CMY astronomical image

# plot rgb image
def plot_rgb(rgb_data, color, figsize = (10, 10), wcs_header = None):
    if wcs_header == None: # without wcs coords
        fig, ax = plt.subplots(figsize = figsize)

        ax.set_xlabel("x", size = 15)
        ax.set_ylabel("y", size = 15)
    
    else: # with wcs coords
        fig, ax = plt.subplots(figsize = figsize, subplot_kw = {'projection': wcs_header})
        ax.coords.grid(color = 'gray', alpha = 0.75, linestyle = 'solid')
        
        ax.set_xlabel("Right Ascension [hms]", size = 15)
        ax.set_ylabel("Declination [degrees]", size = 15)
        
    ax.imshow(rgb_data, vmin = 0, vmax = 1, cmap = color, origin = "lower")

    return fig, ax

def get_filter_throughput(telescope_name, instrument_name, filter_names):
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
            st.image("filter_throughput/nircam_filters.png", caption = "JWST/NIRCAM filter throughput.")
        elif instrument == "MIRI":
            st.image("filter_throughput/miri_filters.png", caption = "JWST/MIRI filter throughput.")
        elif instrument == "NIRISS":
            st.image("filter_throughput/niriss_filters.png", caption = "JWST/MIRI filter throughput.")

    else:
        instrument = instrument_name
        if instrument == "WFC3/UVIS":
            row_indexes = np.where(np.isin(wfc3_uvis.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/wfc3_uvis_filters" + str(ii) + ".png"
                st.image(image_name, caption = "HST/WFC3/UVIS filter throughput.")

        if instrument == "WFC3/IR":
            row_indexes = np.where(np.isin(wfc3_ir.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/wfc3_ir_filters" + str(ii) + ".png"
                st.image(image_name, caption = "HST/WFC3/IR filter throughput.")
        
        if instrument == "ACS/WFC":
            row_indexes = np.where(np.isin(acs_wfc.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/acs_wfc_filters" + str(ii) + ".png"
                st.image(image_name, caption = "HST/ACS/WFC filter throughput.")
        
        if instrument == "ACS/HRC":
            row_indexes = np.where(np.isin(acs_hrc.values, filter_names))[0] + 1 # get indexed values
            for ii in np.unique(row_indexes):
                image_name = "filter_throughput/acs_hrc_filters" + str(ii) + ".png"
                st.image(image_name, caption = "HST/ACS/HRC filter throughput.")

        if instrument == "ACS/SBC":
            st.image("filter_throughput/acs_sbc_filters.png", caption = "HST/ACS/SBC filter throughput.")
    
    return instrument