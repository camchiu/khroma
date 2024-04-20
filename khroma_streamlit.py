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

st.write("Astronomical images are one of the ways we can visually come to understand the grandeur of the Universe. The creation of these images takes a lot of scientific and artistic work. Using this page, you can create your own astronomical images using raw filter data from the Mikulski Archive for Space Telescopes (MAST). Start by choosing a space telescope and astronomical object below.")

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

            # get filters
            mask = (observation_table["instrument_name"] == instrument_name) & (observation_table["calib_level"] > 1)
            observation_table = observation_table[mask]
            total_filter_list = list(np.unique(observation_table["filters"]))

            # choose filters
            if len(total_filter_list) == 1:
                st.write(f"Only 1 filter image of {obj_name} is available from {telescope_name}/{instrument_name}, which is not enough to create an RGB image. Please try a different instrument or input a different object.")
                filter_list = None
            
            elif len(total_filter_list) == 2:
                st.write(f"Only 2 filter images of {obj_name} are available from {telescope_name}/{instrument_name}, which is not enough to create an RGB image. Please try a different instrument or input a different object.")
                filter_list = None

            else:
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
    st.write("Click below to see the individual filter images that you have chosen. Remember that to create a RGB (red-green-blue) color image, you will need at least 3 different filters to assign to each of the three different colors.")
    st.button("Show individual filter images.", on_click = clicked, args = [2]) # button 2 with callback function

    st.write("Use the toggles above each individual filter image to adjust the minimum and maximum to your liking. For example, if the background of the image looks too dark, consider increasing the lower bound.")

    st.write("")

    st.write("Grid lines can help orient the viewer.")
    grid = st.checkbox("Add grid lines.")

    st.write("Sometimes it can be helpful for astronomers to use sky coordinates, such as right ascension and declination, to see the spatial spread of images across the sky. These coordinates work like longitude and latitude lines, except projected on the sky.")

    hms = st.checkbox("Use RA/Dec.")

    st.write("If the individual filter images were taken during different observing runs, they might be slightly misaligned. You can use the button below to reproject all the images against a reference filter image.")

    reproj = st.checkbox("Align images.")
    # if aligning, choose which filter to align against
    if reproj:
        cola, colb, colc = st.columns(spec = [0.1, 0.8, 0.2])
        with colb:
            align = st.selectbox(label = "Reference filter:", options = filter_list, index = 0)

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
    st.write("Now we are ready to layer all of the individual filter data and create our RGB image!")

    st.button("Show full RGB image.", on_click = clicked, args = [3]) # button 3

    # if there are enough filters, show user input for creation of RGB image
    if (telescope_choice != None) & (filter_list != None):
        if (len(filter_list) >= 3):
            st.write("Assign each color a filter below. This can be based off the filter passband ranges shown above, or it can be whatever you decide!")

            instrument, image_name, caption_label = imaging.get_filter_throughput(telescope_name, instrument_name, filter_list)

            expander_filters = st.expander(f"Scientists often match filters to colors based on what wavelengths of light they let it. The image below shows the throughput for each of the filters on {instrument_name}.")
            # expander.write("The chart above shows some numbers I picked for you. I rolled actual dice for these, so they're *guaranteed to be random")
            expander_filters.image(image_name, caption = caption_label)

            expander_colors = st.expander(f"Expander about creating CMY images.")
            expander_colors.write("text + CMY image?")

            # assign rgb values
            r_filter = st.multiselect(label = "Select the filter for red:", options = filter_list)

            # red filter
            r_num = len(r_filter)
            if r_num == 1:
                r = vars()["norm_" + r_filter[0] + "_data"]
            
            elif r_num > 1:
                columns_text = st.columns(spec = [0.1, 0.9])
                with columns_text[1]:
                    st.write("Choose the weights for each red filter:")

                rspec_list = [0.1] + [0.9/r_num] * r_num
                columns = st.columns(spec = rspec_list)
                rweight = 0

                for rr in range(r_num):
                    with columns[rr+1]:
                        vars()["rweight_" + r_filter[rr]] = st.number_input(label = r_filter[rr], value = 1/r_num, key = "r_" + r_filter[rr])
                        rweight += vars()["rweight_" + r_filter[rr]]
                
                for rr in range(r_num):
                    rf = r_filter[rr]
                    normalized_weight = vars()["rweight_" + rf] / rweight
                    if rr == 0:
                        r = normalized_weight * vars()["norm_" + rf + "_data"]
                    else:
                        r += normalized_weight * vars()["norm_" + rf + "_data"]

            # green filter
            g_filter = st.multiselect(label = "Select the filter for green:", options = filter_list)

            g_num = len(g_filter)
            if g_num == 1:
                g = vars()["norm_" + g_filter[0] + "_data"]

            if g_num > 1:
                columns_text = st.columns(spec = [0.1, 0.9])
                with columns_text[1]:
                    st.write("Choose the weights for each green filter:")

                gspec_list = [0.1] + [0.9/g_num] * g_num
                columns = st.columns(spec = gspec_list)
                gweight = 0

                for gg in range(g_num):
                    with columns[gg+1]:
                        vars()["gweight_" + g_filter[gg]] = st.number_input(label = g_filter[gg], value = 1/g_num, key = "g_" + g_filter[gg])
                        gweight += vars()["gweight_" + g_filter[gg]]
                
                for gg in range(r_num):
                    gf = g_filter[gg]
                    normalized_weight = vars()["gweight_" + gf] / gweight
                    if gg == 0:
                        g = normalized_weight * vars()["norm_" + gf + "_data"]
                    else:
                        g += normalized_weight * vars()["norm_" + gf + "_data"]

            # blue filter
            b_filter = st.multiselect(label = "Select the filter for blue:", options = filter_list)

            b_num = len(b_filter)
            if b_num == 1:
                b = vars()["norm_" + b_filter[0] + "_data"]

            if b_num > 1:
                columns_text = st.columns(spec = [0.1, 0.9])
                with columns_text[1]:
                    st.write("Choose the weights for each blue filter:")

                bspec_list = [0.1] + [0.9/b_num] * b_num
                columns = st.columns(spec = bspec_list)
                bweight = 0

                for bb in range(b_num):
                    with columns[bb+1]:
                        vars()["bweight_" + b_filter[bb]] = st.number_input(label = b_filter[bb], value = 1/b_num, key = "b_" + b_filter[bb])
                        bweight += vars()["bweight_" + b_filter[bb]]
                
                for bb in range(b_num):
                    bf = b_filter[bb]
                    normalized_weight = vars()["bweight_" + bf] / bweight
                    if bb == 0:
                        b = normalized_weight * vars()["norm_" + bf + "_data"]
                    else:
                        b += normalized_weight * vars()["norm_" + bf + "_data"]

            st.write("You can also adjust the stretch factor and Q-factor of the image.")

            # get stretch, Q-factor
            s_value = st.slider(label = "Stretch factor", min_value = 0.01, max_value = 2., step = 0.01, value = 0.5)
            q_value = st.slider(label = "Q-factor", min_value = 0.1, max_value = 20., step = 0.01, value = 5.)

            st.write("Don't forget to download your finished product at the end!")

if st.session_state.clicked[3]: # button 3
    with col2:
        # create rgb image
        try:
            rgb = make_lupton_rgb(r, g, b, stretch = s_value, Q = q_value)
            fig, ax = imaging.plot_rgb(rgb)
            st.pyplot(fig)
        except NameError:
            st.error("Make sure that you have assigned a filter to each color.")
        except ValueError:
            st.error("Make sure that the images are aligned. You can use the 'Align images.' checkbox above to do so manually. If that does not work, there might be a problem with the fits header.")

        # download image
        image_name = obj_name + "_" + instrument + ".png"
        plt.savefig(image_name, dpi = 600)
        st.download_button(label = "Download image.", data = open(image_name, "rb").read(), file_name = image_name, mime = "image/png")

### more information about object + observing run (PI)