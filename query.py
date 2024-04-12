# import necessary packages
from astroquery.mast import Observations
from reproject import reproject_interp
from astropy.io import fits
import streamlit as st

# query filters
@st.cache_resource
def search_all_images(filter_name, telescope, obj):
    """
    Parameters
        obj_name (str) : object name to query
        query_filter (str) : JWST filter to query
    Output
        file_name (str) : name of downloaded file
    """
    obs = Observations.query_criteria(obs_collection = telescope, filters = filter_name, target_name = obj, intentType = "science", dataproduct_type = "image")
    prod = Observations.get_product_list(obs[0])
    if telescope == "JWST":
        im = prod[(prod['productType'] == 'SCIENCE') & (prod['productSubGroupDescription'] == "I2D") & (prod['calib_level'] > 2)]
    else:
        im = prod[prod['productType'] == 'SCIENCE']

    if len(im) != 1:
        st.write("Error !!!")

    return im

# download data
@st.cache_data
def download_image(_im, extension = 0):
    # choose specific image
    image = _im[extension]
    
    # download file and return filename
    download = Observations.download_products(_im, cache = False)
    file_name = download['Local Path'][0]
    
    return file_name

# open data
@st.cache_data
def load_fits(filename):
    """
    Parameters
    ----------
        filename (str) : filepath to a fits file
        extension (int) : extension of fits file to load, 0 (default)
    Returns
    -------
        file_header : header of given extension fits file
        file_data : data of given extension of fits file
    """
    hdu = fits.open(filename, ignore_missing_simple = True)

    try:
        file_header = hdu[1].header
        file_data = hdu[1].data
    except:
        file_header = hdu[0].header
        file_data = hdu[0].data
    
    return file_header, file_data

# align images using WCS / headers
@st.cache_data
def repoject(file_name, _match_header):
    """Reproject data using a reference header."""
    header, data = load_fits(file_name)
    original_data = (data, header)
    reprojected_data, _ = reproject_interp(input_data = original_data, output_projection = _match_header)
    # reprojected_data, _ = reproject_interp(input_data = fits.open(vars()["file_" + ff])[1], output_projection = match_header)

    return reprojected_data

# create sky map of available observations ??