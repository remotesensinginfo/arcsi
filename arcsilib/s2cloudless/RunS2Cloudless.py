import rsgislib
import rsgislib.imageutils
import rsgislib.imagecalc
import rsgislib.rastergis
import rsgislib.segmentation
import numpy
from rios import applier
from rios import cuiprogress
from .S2PixelCloudDetector import S2PixelCloudDetector

# Using python-fmask (http://pythonfmask.org)
import fmask.config
import fmask.fmask

import os.path

def run_s2cloudless(input_img, out_prob_img, out_cloud_msk, gdalformat, tmp_base_dir, toa_scale_factor=10000.0, min_obj_size=8):
    """
Function which runs the S2Cloudless methods in 256 x 256 pixel blocks 
across an image.

:param input_img: input sentinel-2 image with all 13 bands.
:param out_prob_img: the output probability of cloud image
:param out_cloud_msk: the output image cloud mask
:param gdalformat: the GDAL image file format of the output image files.
:param tmp_base_dir: a file path temporary intermediate files.
:param toa_scale_factor: scale factor to get TOA between 0-1.
:param min_obj_size: the minimum size of a cloud feature in pixels (default 8 pixels), zero to ignore.

"""
    s2_pxl_cloud_detect = S2PixelCloudDetector(all_bands=True)

    out_tmp_cloud_msk = os.path.join(tmp_base_dir, 'tmp_init_s2less_cloud_msk.kea')

    try:
        import tqdm
        progress_bar = rsgislib.TQDMProgressBar()
    except:
        progress_bar = cuiprogress.GDALProgressBar()

    infiles = applier.FilenameAssociations()
    infiles.s2image = input_img
    outfiles = applier.FilenameAssociations()
    outfiles.out_prob_img = out_prob_img
    outfiles.out_cloud_msk = out_tmp_cloud_msk
    otherargs = applier.OtherInputs()
    otherargs.s2_pxl_cloud_detect = s2_pxl_cloud_detect
    otherargs.toa_scale_factor = toa_scale_factor
    aControls = applier.ApplierControls()
    aControls.progress = progress_bar
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

    out_cloud_closed_msk = os.path.join(tmp_base_dir, 'tmp_init_s2less_cloud_msk_clumps_closed.kea')
    out_cloud_closed_msk_tmp = os.path.join(tmp_base_dir, 'tmp_init_s2less_cloud_msk_clumps_closed_tmp.kea')

    morph_operator = os.path.join(tmp_base_dir, 'morph_circ5')
    morph_operator_file = '{}.gmtxt'.format(morph_operator)
    morph_op_size = 5
    rsgislib.imagemorphology.createCircularOp(morph_operator, morph_op_size)

    rsgislib.imagemorphology.imageClosing(out_tmp_cloud_msk, out_cloud_closed_msk, out_cloud_closed_msk_tmp,
                                          morph_operator_file, True, morph_op_size, 'KEA', rsgislib.TYPE_8UINT)
    out_tmp_cloud_msk = out_cloud_closed_msk

    if min_obj_size > 0:
        out_tmp_cloud_msk_clumps = os.path.join(tmp_base_dir, 'tmp_init_s2less_cloud_msk_clumps.kea')
        out_tmp_cloud_msk_clumps_rmsml = os.path.join(tmp_base_dir, 'tmp_init_s2less_cloud_msk_clumps_rmsml.kea')
        rsgislib.segmentation.clump(out_tmp_cloud_msk, out_tmp_cloud_msk_clumps, 'KEA', False, 0, False)
        rsgislib.rastergis.populateStats(clumps=out_tmp_cloud_msk_clumps, addclrtab=True, calcpyramids=False, ignorezero=True)
        rsgislib.segmentation.rmSmallClumps(out_tmp_cloud_msk_clumps, out_tmp_cloud_msk_clumps_rmsml, min_obj_size, 'KEA')
        rsgislib.imagecalc.imageMath(out_tmp_cloud_msk_clumps_rmsml, out_cloud_msk, 'b1>0?1:0', gdalformat, rsgislib.TYPE_8UINT)
    else:
        rsgislib.imageutils.gdal_translate(out_tmp_cloud_msk, out_cloud_msk, gdalformat)


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
    fmaskConfig.setCloudBufferSize(10)
    fmaskConfig.setShadowBufferSize(10)

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


def run_fmask_cloud_msk(sen2_toa_img, sen2_sat_img, sen2_view_angles_img, fmask_cloud_out_msk, tmp_base_dir, toa_img_scale_factor):
    """

    :param sen2_toa_img:
    :param sen2_sat_img:
    :param sen2_view_angles_img:
    :param fmask_cloud_out_msk:
    :param tmp_base_dir:
    :param toa_img_scale_factor:

    """
    anglesInfo = fmask.config.AnglesFileInfo(sen2_view_angles_img, 3, sen2_view_angles_img, 2, sen2_view_angles_img, 1,
                                             sen2_view_angles_img, 0)

    fmaskFilenames = fmask.config.FmaskFilenames()
    fmaskFilenames.setTOAReflectanceFile(sen2_toa_img)
    fmaskFilenames.setSaturationMask(sen2_sat_img)
    fmaskFilenames.setOutputCloudMaskFile(fmask_cloud_out_msk)

    fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_SENTINEL2)
    fmaskConfig.setAnglesInfo(anglesInfo)
    fmaskConfig.setKeepIntermediates(True)
    fmaskConfig.setVerbose(True)
    fmaskConfig.setTempDir(tmp_base_dir)
    fmaskConfig.setTOARefScaling(float(toa_img_scale_factor))
    fmaskConfig.setMinCloudSize(8)
    fmaskConfig.setCloudBufferSize(10)
    fmaskConfig.setShadowBufferSize(10)

    missingThermal = True

    print("FMASK: Cloud layer, pass 1")
    (pass1file, Twater, Tlow, Thigh, NIR_17, nonNullCount) = fmask.fmask.doPotentialCloudFirstPass(fmaskFilenames, fmaskConfig, missingThermal)
    print("\tTwater=", Twater, "Tlow=", Tlow, "Thigh=", Thigh, "NIR_17=", NIR_17, "nonNullCount=", nonNullCount)

    print("FMASK: Cloud layer, pass 2")
    (pass2file, landThreshold) = fmask.fmask.doPotentialCloudSecondPass(fmaskFilenames, fmaskConfig, pass1file, Twater, Tlow, Thigh, missingThermal, nonNullCount)
    print("\tlandThreshold=", landThreshold)

    print("FMASK: Cloud layer, pass 3")
    tmp_out_cloud_msk = fmask.fmask.doCloudLayerFinalPass(fmaskFilenames, fmaskConfig, pass1file, pass2file, landThreshold, Tlow, missingThermal)

    rsgislib.imageutils.gdal_translate(tmp_out_cloud_msk, fmask_cloud_out_msk, gdal_format='KEA')