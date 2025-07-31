import numpy as np
from scipy.ndimage import maximum_filter, minimum_filter
from scipy import signal

def gauss_kern(size, sizey=None):
    """Create a normalized, 2d, Gaussian convolution filter 

       size - gridpoint radius of filter
       sizey - if filter is elliptical, y-axis gridpoint radius
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)

    xx, yy = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(xx**2/float(size)+yy**2/float(sizey)))

    return g / g.sum()

def _gaussian_filter(data, size):
    if size == 0:
        return data

    kernel = gauss_kern(size)
    return signal.convolve2d(data, kernel, 'same')

def apply_max_or_min_filter(data, max_size=None, min_size=None):
        """Applies maximum or minimum filter to each ensemble member."""
        if max_size is not None:
            func = maximum_filter
            size = max_size
        elif min_size is not None:
            func = minimum_filter
            size = min_size

        data_filtered = func(data[:,:], size=size)

        return data_filtered

def nmep(self, data, thresh_dict, gauss_sizes, dx, max_sizes=None, min_sizes=None, plot_type='severe'):
        """
        Compute the Neighborhood Maximum Ensemble Probability (NMEP), but 
        as not as described in Sobash and Schwartz (2017). Instead, 
        the Gaussian smoothing applied prior to computing the ensemble 
        probability to maintain sharpness. Additionally, 
        a custom gaussian filter is used that also maintains sharpness. 
    
        TODO: incorporate the ability to compute the 
              neighborhood minimum ensemble probability (e.g.,
              for min UH). 
    
        Parameters
        --------------
        data : dicts of array-like of shape (NE, NY, NX)
            Ensemble data for a different variables stored 
            with a dictionary. 
        thresh : dicts of float
            Threshold used to binarize the ensemble data 
            for each variable. 
        max_sizes/min_sizes : int or list thereof.
            Maximum or Minimum filter size(s)
        gauss_sizes : int or list thereof.
            Gaussian filter size. 
        dx : int
            Grid spacing in km. 

        Returns
        --------------
        nmep_data : dict 
            NMEP for each variable at different neighborhood sizes. 
        """
        if not isinstance(max_sizes, list):
            max_sizes = [max_sizes]

        if not isinstance(gauss_sizes, list):
            gauss_sizes = [gauss_sizes]

        nmep_data = OrderedDict()
        data_vars = data.keys()
        for v in data_vars:
            threshs = thresh_dict[v.split('__')[0]]

            threshs = check_threshs_for_plot_type(threshs, plot_type)

            if isinstance(threshs, dict):
                threshs = threshs[v.split('__')[1]]

            if not isinstance(threshs, list):
                threshs = [threshs]

            for th in threshs:
                for max_size, gauss_size in zip(max_sizes, gauss_sizes):
                    # Apply the maximum value to each member separately.
                    data_maxed = self.apply_max_or_min_filter(data[v], max_size=max_size)
                    # Convert to binary based on the intensity threshold.
                    data_binary = np.where(data_maxed>=th,1,0)
                    # Compute the ensemble probability 
                    #(fraction of members exceeding the threshold).
                    probs = np.mean(data_binary, axis=0)

                    # Apply Gaussian filter to smooth the contours. 
                    probs= _gaussian_filter(probs, size=gauss_size)

                    th_txt = int(th) if th >= 1 else str(th).replace('.', '_')

                    nmep_data[f'{v}__nmep_{int(max_size*dx):02d}km__thresh_{th_txt}'] = probs

        return nmep_data

