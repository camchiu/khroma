# import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import re

import astropy.units as u
from astropy.visualization import make_lupton_rgb
from astroquery.mast import Observations
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning

import query
import cleaning
import imaging

# initial website configs
st.set_page_config(page_title = "Khroma", layout = "wide")

st.title("Khroma")

# introduction
st.write("Astronomical images are one of the ways we can visually come to understand the grandeur of the Universe. \
         The creation of these images involves a combination of scientific and artistic work. Following the steps below, \
         you can create your own color astronomical image from three or more different filter images using \
         James Webb Space Telescope (JWST) data pulled from the [Mikulski Archive for Space Telescopes (MAST)](%s). \
         Start by choosing an astronomical object to image below." % "https://archive.stsci.edu/")

# keeping button clicked on for whole session
if 'clicked' not in st.session_state: # initialize the key in session state
    st.session_state.clicked = {1: False, 2: False, 3: False}

def clicked(button): # function to update the value in session state
    st.session_state.clicked[button] = True

#################### SECTION 1 : querying images ####################
col1, col2 = st.columns(spec = [0.34, 0.66])

with col2:
    # search for objs based on ra/dec
    raw_coords = st.text_input("Search for object names based on RA/DEC [deg]:", placeholder = "e.g. 150, -30")
    if raw_coords != "":
        obj_list = query.search_coords(raw_coords)

        with st.container(height = 200):
            st.write(obj_list)

    st.write("You can also get inspiration from [published JWST images](%s)." % "https://webbtelescope.org/images")

with col1:
    # choose telescope
    #telescope_choice = st.selectbox("Choose a space telescope:",
    #                                ("Hubble Legacy Archive (HLA)", "Hubble Space Telescope (HST)", "James Webb Space Telescope (JWST)"),
    #                                index = None,
    #                                placeholder = "Select telescope ...")
    telescope_choice = "James Webb Space Telescope (JWST)"

    # input object to image
    obj_name = st.text_input("Input an astronomical object:", placeholder = "e.g. NGC-3132")

    # set initial filter list to none
    filter_list = []
    
    # look for instruments
    st.button("Click here to check what instruments are available.", on_click = clicked, args = [1]) # button 1

    if st.session_state.clicked[1]: # button 1 -- conditional based on value in session state, not the output
        telescope_full = re.sub(r'\([^)]*\)', '', telescope_choice) # get full 
        telescope_name = re.findall(r'\((.*?)\)', telescope_choice)[0] # get string in parenthesis to use for query

        # check instrument availability
        observation_table = Observations.query_criteria(obs_collection = telescope_name, target_name = obj_name, 
                                                        intentType = "science", dataproduct_type = "image")                       

        instruments = list(np.unique(observation_table["instrument_name"])) 
        if len(instruments) == 0:
            st.write(f"There are no observations of {obj_name} available from {telescope_name}. \
                     Please try a different telescope or input a different object.")
        
        else:
            instrument_name = st.radio("Please select an available instrument:", options = instruments)

            # get filters
            mask = (observation_table["instrument_name"] == instrument_name) & (observation_table["calib_level"] > 1)
            observation_table = observation_table[mask]
            total_filter_list = list(np.unique(observation_table["filters"]))

            # choose filters
            if len(total_filter_list) == 1:
                st.write(f"Only 1 filter image of {obj_name} is available from {telescope_name}/{instrument_name}, \
                         which is not enough to create an RGB image. Please try a different instrument or input a different object.")
                filter_list = None
            
            elif len(total_filter_list) == 2:
                st.write(f"Only 2 filter images of {obj_name} are available from {telescope_name}/{instrument_name}, \
                         which is not enough to create an RGB image. Please try a different instrument or input a different object.")
                filter_list = None

            else:
                filter_list = st.multiselect(label = f"The following {len(total_filter_list)} filters are available. \
                                             Select which ones to plot:",
                                             options = total_filter_list)
                
                # sort filter list
                for ff in range(len(filter_list)):
                    filt = filter_list[ff]
                    if ";" in filt:
                        new_filter_name = filt.split(";")[1].strip() # get part after semicolon
                        filter_list[ff] = new_filter_name

                filter_list = sorted(filter_list, key = lambda s: int(re.search(r'\d+', s).group()))

                # download filter data
                for filter_name in filter_list:
                    im = query.search_all_images(filter = filter_name, telescope = telescope_name, obj = obj_name)
                    vars()["file_" + filter_name] = query.download_image(im, filter_name, telescope_name, obj_name, 0)
                    st.write("    " + filter_name)

st.markdown("""---""")

#################### SECTION 2 : cleaning images ####################
col1, col2, col3 = st.columns(spec = [0.34, 0.33, 0.33])

with col1:
    st.button("Click here to show individual filter images.", on_click = clicked, args = [2]) # button 2 with callback function

    # explanation
    expander_filter = st.expander("Use the toggles above each individual filter image to adjust the lower and upper bounds.")
    expander_filter.write("These numbers specify the lower and upper bound for clipping and normalizing the data. \
                          In other words, all pixels with values less than the lower bound will be plotted as white \
                          and all pixels with values greater than the higher bound will be plotted as black. \
                          The ideal image will have a white background and an object whose features can be seen by varying shades of gray.")
    expander_filter.image("images/norm_example.png", "Properly normalized image (NIRCAM/F470N of NGC 3132).")
    expander_filter.write("If the background of the image looks too light, consider decreasing the lower bound. \
                          If the whole object is oversaturated in black with no shading, consider increasing the upper bound.")

    ### expander to go between different observing runs to flip through ???

    # plot customizations -- grid and ra/dec coords
    expander_fplot = st.expander("Plot customizations.")
    expander_fplot.write("Grid lines can help visually orient the viewer.")
    grid = expander_fplot.checkbox("Add grid lines.")

    expander_fplot.write("Sometimes it can be helpful for astronomers to use a celestial coordinate system to see \
                         the spatial spread of an object across the sky. The [equatorial coordinate system](%s), \
                         with which an object's location is specified by its right ascension (RA) and declination (DEC), \
                         works like longitude and latitude lines projected on the sky." \
                         % "https://en.wikipedia.org/wiki/Equatorial_coordinate_system")

    hms = expander_fplot.checkbox("Use RA/Dec.")

    st.write("If the individual filter images were taken during different observing runs or reduced differently, \
             they might be slightly misaligned. You can use the button below to reproject all the images against a \
             reference filter image such that they all cover the same spatial area.")

    # align images
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
            try:
                vars()[filt + "_data"] = query.reproject(vars()["file_" + filt], match_header, filt)
            except FITSFixedWarning:
                st.error("The images cannot be aligned due to problems with the FITS headers.")

    else:
        for filt in filter_list:
            vars()[filt + "_header"], vars()[filt + "_data"] = query.load_fits(vars()["file_" + filt])

    # plot single filter images
    with col2:
        col2_index = np.arange(0, len(filter_list), 2)
        for ii in col2_index:
            filt = filter_list[ii]
            filter_data = vars()[filt + "_data"]
       
            if reproj:
                obj_wcs = WCS(match_header)
            elif hms:
                obj_wcs = WCS(vars()[filt + "_header"])
            else:
                obj_wcs = None

            vars()["norm_" + filt + "_data"] = cleaning.show_streamlit(filter_data, filt, wcs_header = obj_wcs, grid = grid)

    with col3:
        col3_index = np.arange(1, len(filter_list), 2)
        for ii in col3_index:
            ff = filter_list[ii]
            filter_data = vars()[ff + "_data"]

            if reproj:
                obj_wcs = WCS(match_header)
            elif hms:
                obj_wcs = WCS(vars()[ff + "_header"])
            else:
                obj_wcs = None

            vars()["norm_" + ff + "_data"] = cleaning.show_streamlit(filter_data, ff, wcs_header = obj_wcs, grid = grid)

st.markdown("""---""")

#################### SECTION 3 : RGB imaging ####################
col1, col2 = st.columns(spec = [0.34, 0.66])

with col1:
    # if there are enough filters selected, show user input for creation of RGB image
    if (len(filter_list) >= 3):
        st.write("Now we are ready to layer all of the individual filter images to create an RGB image!")

        st.write("Assign each color a filter below. This can be based off the filter passband ranges shown below, \
                 or it can be whatever you decide!")

        # explain filter passbands
        instrument, image_name, caption_label = imaging.get_filter_throughput(telescope_name, instrument_name, filter_list)

        expander_filters = st.expander("Filter passband ranges.")
        expander_filters.write(f"Scientists often color filters based on the wavelength range of light \
                               each filter allows to pass through. The image below shows the throughput \
                               for each of the filters on {instrument_name}.")
        expander_filters.image(image_name, caption = caption_label)
        expander_filters.write("Regardless of whether the filters can image visible light, \
                               scientists usually match longer wavelengths to redder colors and \
                               shorter wavelengths to bluer colors to match what our eyes see.")

        # explain rgb imaging
        expander_colors = st.expander(f"Creating RGB images.")
        expander_colors.write("In its simplest form, RGB images are created by layering three distinct images colored red, green, and blue. \
                              To get a more nuanced set of colors, you can also layer an image in multiple colors such that they \
                              are colored in a shade other than red, green, or blue. For example, if you layer an image in both \
                              red and blue, it will show up as magenta.")
        expander_colors.image("images/rgb_venn.jpeg", caption = "Primary colors (red, green, blue) and \
                              secondary colors (yellow, cyan, magenta).")
        expander_colors.write("Using the options below, assign each filter (or multiple filters) to the colors red, green, and blue. \
                              If you have multiple filters assigned to one color, you must specify the comparative weights of each filter. \
                              Additionally, if you are assigning a filter to multiple colors, make sure that it is weighted similarly \
                              in each filter if you want to accurately achieve a secondary color shade.")

        # assign rgb values
        r_filter = st.multiselect(label = "Select the filter(s) for red:", options = filter_list)

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
                    try:
                        r += normalized_weight * vars()["norm_" + rf + "_data"]
                    except ValueError:
                        st.error("Please align the images using the 'Align images.' checkbox above.")

        # green filter
        g_filter = st.multiselect(label = "Select the filter(s) for green:", options = filter_list)

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
                
            for gg in range(g_num):
                gf = g_filter[gg]
                normalized_weight = vars()["gweight_" + gf] / gweight
                if gg == 0:
                    g = normalized_weight * vars()["norm_" + gf + "_data"]
                else:
                    try:
                        g += normalized_weight * vars()["norm_" + gf + "_data"]
                    except ValueError:
                        st.error("Please align the images using the 'Align images.' checkbox above.")

        # blue filter
        b_filter = st.multiselect(label = "Select the filter(s) for blue:", options = filter_list)

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
                    try:
                        b += normalized_weight * vars()["norm_" + bf + "_data"]
                    except ValueError:
                        st.error("Please align the images using the 'Align images.' checkbox above.")

        st.button("Click here to create RGB image.", on_click = clicked, args = [3]) # button 3

        # normalization adjustments -- s and q values
        expander_sq = st.expander("Normalization adjustments.")
        expander_sq.write("Adjust the linear stretch of the image.")
        s_value = expander_sq.slider(label = "Stretch factor", min_value = 0.01, max_value = 2., step = 0.01, value = 0.5)
        expander_sq.write("Adjust the asinh softening parameter.")
        q_value = expander_sq.slider(label = "Q-factor", min_value = 0.1, max_value = 20., step = 0.01, value = 5.)

        # plot customizations
        expander_plot = st.expander("Plot customizations.")
        title_name = expander_plot.text_input("Add a title (e.g. object name) to your image:", placeholder = None)
        plot_filter = expander_plot.checkbox("Add filter names and colors to image.")

        # get color of each filter
        color_list, final_filters = None, None
        if plot_filter:
            final_filters, color_list = [], []
            for ff in range(len(filter_list)):
                filt = filter_list[ff]
                filt_color = "#"
                    
                # get hex color by looping through RGB
                for color_filter in [r_filter, g_filter, b_filter]:
                    if filt in color_filter:
                        filt_color += "FF"
                    else:
                        filt_color += "00"
                    
                if filt_color != "#000000": # filter not used in final image
                    final_filters += [filt]
                    color_list += [filt_color]
            
        st.write("Don't forget to download your finished product at the end! :)")

if st.session_state.clicked[3]: # button 3
    with col2:
        # create rgb image
        try:
            rgb = make_lupton_rgb(r, g, b, stretch = s_value, Q = q_value)
            fig, ax = imaging.plot_rgb(rgb, telescope_full, instrument, final_filters, title = title_name, filter_colors = color_list)
            st.pyplot(fig)
        except NameError:
            st.error("Make sure that you have assigned a filter to each color.")
        except ValueError:
            st.error("Please align the images using the 'Align images.' checkbox above.")

        # download image
        image_name = obj_name + "_" + instrument + ".pdf"
        plt.savefig(image_name, dpi = 600, bbox_inches = "tight")
        st.download_button(label = "Download image.", data = open(image_name, "rb").read(), file_name = image_name, mime = "application/octet-stream")

### more information about object + observing run (PI)