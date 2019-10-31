import rsgislib
import rsgislib.imageutils
import rsgislib.imagecalc
import rsgislib.rastergis
import numpy
from rios import applier
from rios import cuiprogress
from .S2PixelCloudDetector import S2PixelCloudDetector

# Using python-fmask (http://pythonfmask.org)
import fmask.config
import fmask.fmask

import os.path


def run_s2cloudless(input_img, out_prob_img, out_cloud_msk, gdalformat, toa_scale_factor=10000.0):
    """
Function which runs the S2Cloudless methods in 256 x 256 pixel blocks 
across an image.

:param input_img: input sentinel-2 image with all 13 bands.
:param output_img: the output image cloud mask
:param gdalformat: the GDAL image file format of the output image file.

"""
    s2_pxl_cloud_detect = S2PixelCloudDetector(all_bands=True)

    infiles = applier.FilenameAssociations()
    infiles.s2image = input_img
    outfiles = applier.FilenameAssociations()
    outfiles.out_prob_img = out_prob_img
    outfiles.out_cloud_msk = out_cloud_msk
    otherargs = applier.OtherInputs()
    otherargs.s2_pxl_cloud_detect = s2_pxl_cloud_detect
    otherargs.toa_scale_factor = toa_scale_factor
    aControls = applier.ApplierControls()
    aControls.progress = cuiprogress.CUIProgressBar()
    aControls.drivername = gdalformat
    aControls.omitPyramids = True
    aControls.calcStats = False

    def _applyS2Cloudless(info, inputs, outputs, otherargs):
        """
        This is an internal rios function
        """
        # Current shape is: [13 x n x m]
        # Image data needs to be in shape [1 x n x m x 13]
        s2img_reshp = numpy.expand_dims(numpy.stack([inputs.s2image[0], inputs.s2image[1], inputs.s2image[2], inputs.s2image[3], inputs.s2image[4], inputs.s2image[5], inputs.s2image[6], inputs.s2image[7], inputs.s2image[8], inputs.s2image[9], inputs.s2image[10], inputs.s2image[11], inputs.s2image[12]], axis=2), axis=0)
        s2img_reshp_toa = s2img_reshp/otherargs.toa_scale_factor
        outputs.out_prob_img = otherargs.s2_pxl_cloud_detect.get_cloud_probability_maps(s2img_reshp_toa)
        outputs.out_cloud_msk = otherargs.s2_pxl_cloud_detect.get_mask_from_prob(outputs.out_prob_img)

    applier.apply(_applyS2Cloudless, infiles, outfiles, otherargs, controls=aControls)


def run_pyfmask_shadow_masking(sen2_toa_img, sen2_sat_img, sen2_view_angles_img, cloud_msk_img, tmp_base_dir, toa_img_scale_factor, out_cloud_cldshad_msk_img):
    anglesInfo = fmask.config.AnglesFileInfo(sen2_view_angles_img, 3, sen2_view_angles_img, 2, sen2_view_angles_img, 1,
                                             sen2_view_angles_img, 0)

    tmp_base_name = os.path.splitext(os.path.basename(sen2_toa_img))[0]

    fmaskCloudsImg = os.path.join(tmp_base_dir, '{}_pyfmask_clouds_result.kea'.format(tmp_base_name))
    fmaskFilenames = fmask.config.FmaskFilenames()
    fmaskFilenames.setTOAReflectanceFile(sen2_toa_img)
    fmaskFilenames.setSaturationMask(sen2_sat_img)
    fmaskFilenames.setOutputCloudMaskFile(fmaskCloudsImg)

    fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_SENTINEL2)
    fmaskConfig.setAnglesInfo(anglesInfo)
    fmaskConfig.setKeepIntermediates(True)
    fmaskConfig.setVerbose(True)
    fmaskConfig.setTempDir(tmp_base_dir)
    fmaskConfig.setTOARefScaling(float(toa_img_scale_factor))
    fmaskConfig.setMinCloudSize(8)

    missingThermal = True

    print("Cloud layer, pass 1")
    (pass1file, Twater, Tlow, Thigh, NIR_17, nonNullCount) = fmask.fmask.doPotentialCloudFirstPass(fmaskFilenames, fmaskConfig, missingThermal)

    print("Potential shadows")
    potentialShadowsFile = fmask.fmask.doPotentialShadows(fmaskFilenames, fmaskConfig, NIR_17)

    print("Clumping clouds")
    (clumps, numClumps) = fmask.fmask.clumpClouds(cloud_msk_img)

    print("Making 3d clouds")
    (cloudShape, cloudBaseTemp, cloudClumpNdx) = fmask.fmask.make3Dclouds(fmaskFilenames, fmaskConfig, clumps, numClumps, missingThermal)

    print("Making cloud shadow shapes")
    shadowShapesDict = fmask.fmask.makeCloudShadowShapes(fmaskFilenames, fmaskConfig, cloudShape, cloudClumpNdx)

    print("Matching shadows")
    interimShadowmask = fmask.fmask.matchShadows(fmaskConfig, cloud_msk_img, potentialShadowsFile, shadowShapesDict, cloudBaseTemp, Tlow, Thigh, pass1file)

    bandDefns = []
    bandDefns.append(rsgislib.imagecalc.BandDefn('cld', cloud_msk_img, 1))
    bandDefns.append(rsgislib.imagecalc.BandDefn('shd', interimShadowmask, 1))
    rsgislib.imagecalc.bandMath(out_cloud_cldshad_msk_img, 'cld==1?1:shd==1?2:0', 'KEA', rsgislib.TYPE_8UINT, bandDefns)


