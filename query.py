# import necessary packages
from astroquery.mast import Observations
from reproject import reproject_interp
from astropy.io import fits
import streamlit as st

# query filters
@st.cache_resource
def search_all_images(filter, telescope, obj):
    """
    Search for list of observations for a specific filter/telescope/object combination through MAST.
    Parameters
    ----------
        filter (str) : specific filter to query
        telescope (str) : telescope to query
        obj (str) : name of obj to query
    Output
        im (astropy table) : list of queried observations
    """
    obs = Observations.query_criteria(obs_collection = telescope, filters = filter, target_name = obj, intentType = "science", dataproduct_type = "image")
    prod = Observations.get_product_list(obs[0])
    if telescope == "JWST":
        im = prod[(prod['productType'] == 'SCIENCE') & (prod['productSubGroupDescription'] == "I2D") & (prod['calib_level'] > 2)]
    else:
        im = prod[prod['productType'] == 'SCIENCE' & (prod['calib_level'] > 2)]

    if len(im) != 1:
        st.write("Error !!!")

    return im

# download data
@st.cache_data # have to specify filter, telescope, obj to create a unique input to correctly cache
def download_image(_im, filter, telescope, obj, extension = 0):
    """
    Download image for a certain observation.
    Parameters
    ----------
        _im (astropy table) : list of queried observations
        filter (str) : filter name (for caching)
        telescope (str) : telescope name (for caching)
        obj (str) : astronomical object name (for caching)
        extension (int) : queried observation index to download
    Output
    ------
        file_name (str) : filepath to fits file
    """
    # choose specific image
    image = _im[extension]
    
    # download file and return filename
    download = Observations.download_products(image, cache = False)
    file_name = download['Local Path'][0]
    
    return file_name

# open data
@st.cache_data
def load_fits(file_name):
    """
    Parameters
    ----------
        file_name (str) : filepath to fits file
    Returns
    -------
        file_header (fits header) : header of given extension fits file
        file_data (fits data) : data of given extension of fits file
    """
    # open image fits file
    hdu = fits.open(file_name, ignore_missing_simple = True)

    # get header and data from file
    try:
        file_header = hdu[1].header
        file_data = hdu[1].data
    except:
        file_header = hdu[0].header
        file_data = hdu[0].data
    
    return file_header, file_data

# align images using WCS defined by headers
@st.cache_data # have to specify match_filter to create a unique input to correctly cache
def reproject(file_name, _match_header, match_filter):
    """Reproject data using a reference header.
    Parameters
    ----------
        file_name (str) : filepath to fits file to realign
        match_filter (str) : filter name to align against (for caching)
        _match_header (fits header) : header of image to align against
    Output
    ------
        reprojected_data (fits data) : realigned data
    """
    # open data to reproject
    header, data = load_fits(file_name)
    original_data = (data, header)

    # reproject data
    reprojected_data, _ = reproject_interp(input_data = original_data, output_projection = _match_header)

    return reprojected_data

### create sky map of available observations ??