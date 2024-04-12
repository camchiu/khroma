# import necessary packages
import matplotlib.pyplot as plt

# method to create and download a full RGB or RGB + CMY astronomical image

# 
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