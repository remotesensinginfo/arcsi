"""
Module for making pixel-based classification on Sentinel-2 L1C imagery
"""

import copy
import os
import warnings
import numpy as np

from scipy.ndimage.filters import convolve
from skimage.morphology import disk, dilation
import joblib

from .PixelClassifier import PixelClassifier


warnings.filterwarnings("ignore", category=UserWarning)

S2CLOUDLESS_MODEL_FILENAME = 'pixel_s2_cloud_detector_lightGBM_v0.1.joblib.dat'
install_prefix = __file__[:__file__.find('lib')]
DEFAULT_ARCSI_S2CLOUDLESS_MODEL = os.path.join(install_prefix, "share","arcsi", S2CLOUDLESS_MODEL_FILENAME)

MODEL_EVALSCRIPT = 'return [B01,B02,B04,B05,B08,B8A,B09,B10,B11,B12]'
S2_BANDS_EVALSCRIPT = 'return [B01,B02,B03,B04,B05,B06,B07,B08,B8A,B09,B10,B11,B12]'


class S2PixelCloudDetector:
    """
    Sentinel Hub's pixel-based cloud detector for Sentinel-2 imagery.

    Classifier takes as an input Sentinel-2 image of shape n x m x 13 (all 13 bands)
    or n x m x 10 (bands 1, 2, 4, 5, 8, 8A, 9, 10, 11, 12) and returns a raster
    binary cloud mask of shape n x m, where 0 (1) indicates clear sky (cloudy) pixel.
    The classifier can instead of a raster cloud mask return a cloud probability map
    of shape n x m, where each pixel's value is bound between 0 (clear-sky-like pixel)
    and 1 (cloud-like pixel).

    User can control cloud probability threshold and/or post-processing steps -
    convolution with disk (with user defined filter size) and dilation with disk
    (with user defined filter size).

    :param threshold: Cloud probability threshold. All pixels with cloud probability above
                      threshold value are masked as cloudy pixels. Default is 0.4.
    :type threshold: float
    :param all_bands: Flag specifying that input images will consists of all 13 Sentinel-2 bands.
    :type all_bands: bool
    :param average_over: Size of the disk in pixels for performing convolution (averaging probability
                         over pixels). Value 0 means do not perform this post-processing step.
                         Default is 1.
    :type average_over: int
    :param dilation_size: Size of the disk in pixels for performing dilation. Value 0 means do not perform
                          this post-processing step. Default is 1.
    :type dilation_size: int
    :param model_filename: Location of the serialised model. If None the default model provided with the
                           package is loaded.
    :type model_filename: str or None
    """
    BAND_IDXS = [0, 1, 3, 4, 7, 8, 9, 10, 11, 12]

    # pylint: disable=invalid-name
    def __init__(self, threshold=0.4, all_bands=False, average_over=1, dilation_size=1, model_filename=None):

        self.threshold = threshold
        self.all_bands = all_bands
        self.average_over = average_over
        self.dilation_size = dilation_size

        if model_filename is None:
            model_filename = DEFAULT_ARCSI_S2CLOUDLESS_MODEL

        self._load_classifier(model_filename)

        if average_over > 0:
            self.conv_filter = disk(average_over) / np.sum(disk(average_over))

        if dilation_size > 0:
            self.dilation_filter = disk(dilation_size)

    def _load_classifier(self, filename):
        """
        Loads the classifier.
        """
        self.classifier = PixelClassifier(joblib.load(filename))

    def get_cloud_probability_maps(self, X):
        """
        Runs the cloud detection on the input images (dimension n_images x n x m x 10
        or n_images x n x m x 13) and returns an array of cloud probability maps (dimension
        n_images x n x m). Pixel values close to 0 indicate clear-sky-like pixels, while
        values close to 1 indicate pixels covered with clouds.

        :param X: input Sentinel-2 image obtained with Sentinel-Hub's WMS/WCS request
                  (see https://github.com/sentinel-hub/sentinelhub-py)
        :type X: numpy array (shape n_images x n x m x 10 or n x m x 13)
        :return: cloud probability map
        :rtype: numpy array (shape n_images x n x m)
        """
        band_num = X.shape[-1]
        exp_bands = 13 if self.all_bands else len(self.BAND_IDXS)
        if band_num != exp_bands:
            raise ValueError("Parameter 'all_bands' is set to {}. Therefore expected band data with {} bands, "
                             "got {} bands".format(self.all_bands, exp_bands, band_num))

        if self.all_bands:
            X = X[..., self.BAND_IDXS]

        return self.classifier.image_predict_proba(X)[..., 1]

    def get_cloud_masks(self, X):
        """
        Runs the cloud detection on the input images (dimension n_images x n x m x 10
        or n_images x n x m x 13) and returns the raster cloud mask (dimension n_images x n x m).
        Pixel values equal to 0 indicate pixels classified as clear-sky, while values
        equal to 1 indicate pixels classified as clouds.

        :param X: input Sentinel-2 image obtained with Sentinel-Hub's WMS/WCS request
                  (see https://github.com/sentinel-hub/sentinelhub-py)
        :type X: numpy array (shape n_images x n x m x 10 or n x m x 13)
        :return: raster cloud mask
        :rtype: numpy array (shape n_images x n x m)
        """

        cloud_probs = self.get_cloud_probability_maps(X)

        return self.get_mask_from_prob(cloud_probs)

    def get_mask_from_prob(self, cloud_probs, threshold=None):
        """
        Returns cloud mask by applying morphological operations -- convolution and dilation --
        to input cloud probabilities.

        :param cloud_probs: cloud probability map
        :type cloud_probs: numpy array of cloud probabilities (shape n_images x n x m)
        :param threshold: A float from [0,1] specifying threshold
        :type threshold: float
        :return: raster cloud mask
        :rtype: numpy array (shape n_images x n x m)
        """
        threshold = self.threshold if threshold is None else threshold

        if self.average_over:
            cloud_masks = np.asarray([convolve(cloud_prob, self.conv_filter) > threshold
                                      for cloud_prob in cloud_probs], dtype=np.uint8)
        else:
            cloud_masks = (cloud_probs > threshold).astype(np.uint8)

        if self.dilation_size:
            cloud_masks = np.asarray([dilation(cloud_mask, self.dilation_filter) for cloud_mask in cloud_masks],
                                     dtype=np.uint8)

        return cloud_masks



