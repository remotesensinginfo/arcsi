"""
This module lists all externally useful classes and functions
"""

#from .S2PixelCloudDetector import S2PixelCloudDetector, MODEL_EVALSCRIPT, S2_BANDS_EVALSCRIPT
#from .PixelClassifier import PixelClassifier
#from .RunS2Cloudless import run_s2cloudless

import os

import rsgislib
import rsgislib.imagecalc
import rsgislib.imageutils
import rsgislib.classification.classlightgbm
import rsgislib.imagemorphology
import rsgislib.segmentation
import rsgislib.rastergis

# Using python-fmask (http://pythonfmask.org)
import fmask.config
import fmask.fmask

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

S2CLOUDLESS_MODEL_FILENAME = 'pixel_s2_cloud_detector_lightGBM_v0.1.txt'
install_prefix = __file__[:__file__.find('lib')]
DEFAULT_ARCSI_S2CLOUDLESS_MODEL = os.path.join(install_prefix, "share","arcsi", S2CLOUDLESS_MODEL_FILENAME)


def run_fmask_cloud_msk(sen2_toa_img, sen2_sat_img, sen2_view_angles_img, fmask_cloud_out_msk,
                        tmp_base_dir, toa_img_scale_factor, use_frantz_disp=False):
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
    if use_frantz_disp:
        fmaskConfig.setSen2displacementTest(True)  # Frantz et al implementation.
    else:
        fmaskConfig.setSen2displacementTest(False)  # Frantz et al implementation.

    missingThermal = True

    print("FMASK: Cloud layer, pass 1")
    (pass1file, Twater, Tlow, Thigh, NIR_17, nonNullCount) = fmask.fmask.doPotentialCloudFirstPass(fmaskFilenames,
                                                                                                   fmaskConfig,
                                                                                                   missingThermal)
    print("\tTwater={} Tlow={} Thigh= {} NIR_17={} nonNullCount={}".format(Twater, Tlow, Thigh, NIR_17, nonNullCount))

    print("FMASK: Cloud layer, pass 2")
    (pass2file, landThreshold) = fmask.fmask.doPotentialCloudSecondPass(fmaskFilenames, fmaskConfig, pass1file, Twater,
                                                                        Tlow, Thigh, missingThermal, nonNullCount)
    print("\tlandThreshold={}".format(landThreshold))

    print("FMASK: Cloud layer, pass 3")
    tmp_out_cloud_msk = fmask.fmask.doCloudLayerFinalPass(fmaskFilenames, fmaskConfig, pass1file, pass2file,
                                                          landThreshold, Tlow, missingThermal)

    rsgislib.imageutils.gdal_translate(tmp_out_cloud_msk, fmask_cloud_out_msk, gdal_format='KEA')


def run_pyfmask_shadow_masking(sen2_toa_img, sen2_sat_img, sen2_view_angles_img, cloud_msk_img, tmp_base_dir,
                               toa_img_scale_factor, out_cloud_cldshad_msk_img):
    """

    :param sen2_toa_img:
    :param sen2_sat_img:
    :param sen2_view_angles_img:
    :param cloud_msk_img:
    :param tmp_base_dir:
    :param toa_img_scale_factor:
    :param out_cloud_cldshad_msk_img:

    """
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
    (pass1file, Twater, Tlow, Thigh, NIR_17, nonNullCount) = fmask.fmask.doPotentialCloudFirstPass(fmaskFilenames,
                                                                                                   fmaskConfig,
                                                                                                   missingThermal)

    print("Potential shadows")
    potentialShadowsFile = fmask.fmask.doPotentialShadows(fmaskFilenames, fmaskConfig, NIR_17)

    print("Clumping clouds")
    (clumps, numClumps) = fmask.fmask.clumpClouds(cloud_msk_img)

    print("Making 3d clouds")
    (cloudShape, cloudBaseTemp, cloudClumpNdx) = fmask.fmask.make3Dclouds(fmaskFilenames, fmaskConfig, clumps,
                                                                          numClumps, missingThermal)

    print("Making cloud shadow shapes")
    shadowShapesDict = fmask.fmask.makeCloudShadowShapes(fmaskFilenames, fmaskConfig, cloudShape, cloudClumpNdx)

    print("Matching shadows")
    interimShadowmask = fmask.fmask.matchShadows(fmaskConfig, cloud_msk_img, potentialShadowsFile, shadowShapesDict,
                                                 cloudBaseTemp, Tlow, Thigh, pass1file)

    bandDefns = []
    bandDefns.append(rsgislib.imagecalc.BandDefn('cld', cloud_msk_img, 1))
    bandDefns.append(rsgislib.imagecalc.BandDefn('shd', interimShadowmask, 1))
    rsgislib.imagecalc.bandMath(out_cloud_cldshad_msk_img, 'cld==1?1:shd==1?2:0', 'KEA', rsgislib.TYPE_8UINT, bandDefns)


def run_s2cloudless(s2img, out_cloud_msk, s2_vmsk_img, gdalformat, tmp_dir, toa_scale_factor=10000.0, min_obj_size=10,
                    morph_close_size=5, morph_dilate_size=9):
    """

    :param s2img:
    :param out_cloud_msk:
    :param s2_vmsk_img:
    :param gdalformat:
    :param tmp_dir:
    :param toa_scale_factor:
    :param min_obj_size:
    :param morph_close_size:
    :param morph_dilate_size:

    """
    basename = os.path.splitext(os.path.basename(s2img))[0]
    s2_tmp_img = os.path.join(tmp_dir, '{}_s2img_flt.kea'.format(basename))
    rsgislib.imagecalc.imageMath(s2img, s2_tmp_img, 'b1/{}'.format(toa_scale_factor), 'KEA', rsgislib.TYPE_32FLOAT)

    imgs_info = [
        rsgislib.imageutils.ImageBandInfo(fileName=s2_tmp_img, name='sen2', bands=[1, 2, 4, 5, 8, 9, 10, 11, 12, 13])]

    out_score_img = os.path.join(tmp_dir, '{}_s2cls_tmp_score.kea'.format(basename))
    out_cls_img = os.path.join(tmp_dir, '{}_s2cls_tmp_cls.kea'.format(basename))
    rsgislib.classification.classlightgbm.apply_lightgbm_binary_classifier(DEFAULT_ARCSI_S2CLOUDLESS_MODEL, s2_vmsk_img,
                                                                           1, imgs_info, out_score_img, 'KEA',
                                                                           out_cls_img, class_thres=5000)

    morph_close_opt = os.path.join(tmp_dir, 'morph_circ_close')
    morph_close_opt_file = '{}.gmtxt'.format(morph_close_opt)
    rsgislib.imagemorphology.createCircularOp(morph_close_opt, morph_close_size)

    out_cls_morph_img = os.path.join(tmp_dir, '{}_s2cls_tmp_cls_morph.kea'.format(basename))
    out_cls_morph_tmp_img = os.path.join(tmp_dir, '{}_s2cls_tmp_cls_morph_tmp.kea'.format(basename))
    rsgislib.imagemorphology.imageClosing(out_cls_img, out_cls_morph_img, out_cls_morph_tmp_img,
                                          morph_close_opt_file, True, morph_close_size, 'KEA', rsgislib.TYPE_8UINT)

    out_cls_img_clumps = os.path.join(tmp_dir, '{}_s2cls_tmp_cls_clumps.kea'.format(basename))
    out_cls_img_clumps_rmsml = os.path.join(tmp_dir, '{}_s2cls_tmp_cls_clumps_rmsml.kea'.format(basename))
    rsgislib.segmentation.clump(out_cls_morph_img, out_cls_img_clumps, 'KEA', False, 0, False)
    rsgislib.rastergis.populateStats(clumps=out_cls_img_clumps, addclrtab=True, calcpyramids=False, ignorezero=True)
    rsgislib.segmentation.rmSmallClumps(out_cls_img_clumps, out_cls_img_clumps_rmsml, min_obj_size, 'KEA')

    out_cloud_cls_img = os.path.join(tmp_dir, '{}_s2cls_tmp_cls_clumps.kea'.format(basename))
    rsgislib.imagecalc.imageMath(out_cls_img_clumps_rmsml, out_cloud_cls_img, 'b1>0?1:0', gdalformat,
                                 rsgislib.TYPE_8UINT)

    morph_dilate_opt = os.path.join(tmp_dir, 'morph_circ_dilate')
    morph_dilate_opt_file = '{}.gmtxt'.format(morph_dilate_opt)
    rsgislib.imagemorphology.createCircularOp(morph_dilate_opt, morph_dilate_size)
    rsgislib.imagemorphology.imageDilate(out_cloud_cls_img, out_cloud_msk, morph_dilate_opt_file,
                                         True, morph_dilate_size, 'KEA', rsgislib.TYPE_8UINT)


