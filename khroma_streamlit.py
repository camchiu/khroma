# import necessary packages
import numpy as np
import streamlit as st
import re

from astropy.visualization import make_lupton_rgb
from astroquery.mast import Observations
from astropy.wcs import WCS
import matplotlib.pyplot as plt

import query
import cleaning
import imaging

# initial website configs
st.set_page_config(page_title = "Khroma", layout = "wide")

st.title("Khroma")

st.write("need at least 3 filter images.")

# keeping button clicked on
if 'clicked' not in st.session_state: # initialize the key in session state
    st.session_state.clicked = {1: False, 2: False, 3: False}

def clicked(button): # function to update the value in session state
    st.session_state.clicked[button] = True

#################### SECTION 1 : querying images ####################
col1, col2 = st.columns(spec = [0.34, 0.66])

### GALEX - FUV, NUV (2 filters); HLSP - unwell fits files; HUT ???; OPO - public outreach (see actual STScI image!); PS1 - unwell fits files

with col1:
    # choose telescope
    telescope_choice = st.selectbox("Choose a space telescope:",
                                    ("Hubble Legacy Archive (HLA)", "Hubble Space Telescope (HST)", "James Webb Space Telescope (JWST)"),
                                    index = None,
                                    placeholder = "Select telescope ...")

    # input object to image
    obj_name = st.text_input("Input an astronomical object:", placeholder = "e.g. NGC-4303")

    # set initial filter list to none
    filter_list = None
    
    # streamlit
    st.button("Check instrument availability.", on_click = clicked, args = [1]) # button 1 with callback function

    if st.session_state.clicked[1]: # button 1 -- conditional based on value in session state, not the output
        telescope_name = re.findall(r'\((.*?)\)', telescope_choice)[0] # find string in parenthesis to use for query

        # check instrument availability
        observation_table = Observations.query_criteria(obs_collection = telescope_name, target_name = obj_name, 
                                                        intentType = "science", dataproduct_type = "image")

        instruments = list(np.unique(observation_table["instrument_name"])) 
        if len(instruments) == 0:
            st.write(f"There are no observations of {obj_name} available from {telescope_name}. Please try a different telescope or input a different object.")
        else:
            instrument_name = st.radio("Please select an available instrument:", options = instruments)

        # choose filters
        mask = (observation_table["instrument_name"] == instrument_name) & (observation_table["calib_level"] > 1)
        observation_table = observation_table[mask]
        total_filter_list = list(np.unique(observation_table["filters"]))

        filter_list = st.multiselect(label = f"The following {len(total_filter_list)} filters are available. Select which ones to plot:",
                                     options = total_filter_list)

        # download filter data
        for filter_name in filter_list:
            im = query.search_all_images(filter = filter_name, telescope = telescope_name, obj = obj_name)
            vars()["file_" + filter_name] = query.download_image(im, filter_name, telescope_name, obj_name, 0)
            st.write("    " + filter_name)

with col2:
    ### sky map
    st.write(" sky map tbd ... working on it ... ")

st.markdown("""---""")

#################### SECTION 2 : cleaning images ####################
col1, col2, col3 = st.columns(spec = [0.34, 0.33, 0.33])
### check that all filters are from the same observing run ?? / incorporate observing run somehow ?? -- toggle through different available images ??? -- FLIP THROUGH EACH OF THE IMAGES WHEN DOWNLOADED !!!

with col1:
    st.button("Show individual filter images.", on_click = clicked, args = [2]) # button 2 with callback function

    st.write(" ... text ...")

    hms = st.checkbox("Use RA/Dec.")

    grid = st.checkbox("Add grid lines.")

    reproj = st.checkbox("Align images.")
    # if aligning, choose which filter to align against
    if reproj:
        cola, colb, colc = st.columns(spec = [0.1, 0.8, 0.2])
        with colb:
            align = st.selectbox(label = "Select a filter to align against ...", options = filter_list, index = 0)

if st.session_state.clicked[2]: # button 2
    # open filter data
    if reproj:
        match_header, vars()[align + "_data"] = query.load_fits(vars()["file_" + align])

        for filt in filter_list:
            vars()[filt + "_data"] = query.reproject(vars()["file_" + filt], match_header, filt)

    else:
        for filt in filter_list:
            vars()[filt + "_header"], vars()[filt + "_data"] = query.load_fits(vars()["file_" + filt])

    # plot single filter images
    with col2:
        col2_index = np.arange(0, len(filter_list), 2)
        for ii in col2_index:
            filt = filter_list[ii]
            filter_data = vars()[filt + "_data"]
            key_names = [filt + "_key1", filt + "_key2"]
       
            if reproj:
                obj_wcs = WCS(match_header)
            elif hms:
                obj_wcs = WCS(vars()[filt + "_header"])
            else:
                obj_wcs = None

            vars()["norm_" + filt + "_data"] = cleaning.show_streamlit(filter_data, filt, key_names, wcs_header = obj_wcs, grid = grid)

    with col3:
        col3_index = np.arange(1, len(filter_list), 2)
        for ii in col3_index:
            ff = filter_list[ii]
            filter_data = vars()[ff + "_data"]
            key_names = [ff + "_key1", ff + "_key2"]

            if reproj:
                obj_wcs = WCS(match_header)
            elif hms:
                obj_wcs = WCS(vars()[ff + "_header"])
            else:
                obj_wcs = None

            vars()["norm_" + ff + "_data"] = cleaning.show_streamlit(filter_data, ff, key_names, wcs_header = obj_wcs, grid = grid)

st.markdown("""---""")

#################### SECTION 3 : RGB imaging ####################
col1, col2 = st.columns(spec = [0.34, 0.66])

with col1:
    st.button("Show full RGB image.", on_click = clicked, args = [3]) # button 3

    st.write("explanatory text ....")

    # if there are enough filters, show user input for creation of RGB image
    if (telescope_choice != None) & (filter_list != None):
        if (len(filter_list) >= 3):
            instrument = imaging.get_filter_throughput(telescope_name, instrument_name, filter_list)

            # assign rgb values
            r_filter = st.selectbox(label = "Select the filter for red ...", options = filter_list, index = 0)
            g_filter = st.selectbox(label = "Select the filter for green ...", options = filter_list, index = 1)
            b_filter = st.selectbox(label = "Select the filter for blue ...", options = filter_list, index = 2)

            st.write("explanatory text ....")

            # get stretch, Q-factor
            s_value = st.slider(label = "Stretch factor", min_value = 0.01, max_value = 2., step = 0.01, value = 0.5)
            q_value = st.slider(label = "Q-factor", min_value = 0.1, max_value = 20., step = 0.01, value = 5.)

if st.session_state.clicked[3]: # button 3
    with col2:
        # get normalized data
        r = vars()["norm_" + r_filter + "_data"]
        g = vars()["norm_" + g_filter + "_data"]
        b = vars()["norm_" + b_filter + "_data"]

        # create rgb image
        rgb = make_lupton_rgb(r, g, b, stretch = s_value, Q = q_value)
        fig, ax = imaging.plot_rgb(rgb, color = "plasma", wcs_header = obj_wcs)
        st.pyplot(fig)

        # download image
        image_name = obj_name + "_" + instrument + ".png"
        plt.savefig(image_name, dpi = 600)
        st.download_button(label = "Download image.", data = open(image_name, "rb").read(), file_name = image_name, mime = "image/png")

### more information about object + observing run (PI)










# expander = st.expander("See explanation")
# expander.write(\"\"\"
#     The chart above shows some numbers I picked for you.
#     I rolled actual dice for these, so they're *guaranteed* to
#     be random.
# \"\"\")
# expander.image("https://static.streamlit.io/examples/dice.jpg")