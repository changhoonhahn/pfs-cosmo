SELECT
 meas.object_id
, meas.parent_id
, meas.tract
, meas.patch
, meas.i_ra
, meas.i_dec
, meas.i_variance_value
, meas.i_footprintarea_value
,meas.skymap_id
,meas2.i_psfflux_mag as meas_i_psfflux_mag

--patch
,patch_qa.patch
,patch_qa.skymap_id
,patch_qa.ra
,patch_qa.dec

--depth
,patch_qa.gmag_psf_depth  as g_depth
,patch_qa.rmag_psf_depth  as r_depth
,patch_qa.imag_psf_depth  as i_depth
,patch_qa.zmag_psf_depth  as z_depth
,patch_qa.ymag_psf_depth  as y_depth

--psf size
,patch_qa.gseeing
,patch_qa.rseeing
,patch_qa.iseeing
,patch_qa.zseeing
,patch_qa.yseeing
    
-- forced measurement (extinction needed, subtract this to get magnitude)
, forced.a_g
, forced.a_r
,forced.a_i
, forced.a_z
, forced.a_y

-- forced CModel magnitudes and fluxes (needed)
, forced.g_cmodel_mag       as forced_g_cmodel_mag
, forced.g_cmodel_magerr    as forced_g_cmodel_magerr

, forced.r_cmodel_mag       as forced_r_cmodel_mag
, forced.r_cmodel_magerr    as forced_r_cmodel_magerr

, forced.i_cmodel_mag       as forced_i_cmodel_mag
, forced.i_cmodel_magerr    as forced_i_cmodel_magerr

, forced.z_cmodel_mag       as forced_z_cmodel_mag
, forced.z_cmodel_magerr    as forced_z_cmodel_magerr
    
, forced.y_cmodel_mag       as forced_y_cmodel_mag
, forced.y_cmodel_magerr    as forced_y_cmodel_magerr

--meas Cmodel mag
, meas.g_cmodel_mag       as meas_g_cmodel_mag
, meas.g_cmodel_magerr    as meas_g_cmodel_magerr

, meas.r_cmodel_mag       as meas_r_cmodel_mag
, meas.r_cmodel_magerr    as meas_r_cmodel_magerr

, meas.i_cmodel_mag       as meas_i_cmodel_mag
, meas.i_cmodel_magerr    as meas_i_cmodel_magerr

, meas.z_cmodel_mag       as meas_z_cmodel_mag
, meas.z_cmodel_magerr    as meas_z_cmodel_magerr
    
, meas.y_cmodel_mag       as meas_y_cmodel_mag
, meas.y_cmodel_magerr    as meas_y_cmodel_magerr


--forced psf
--, forced2.g_psfflux_mag as forced_g_psf_mag
--, forced2.g_psfflux_magerr as forced_g_psf_magerr

--, forced2.r_psfflux_mag as forced_r_psf_mag
--, forced2.r_psfflux_magerr as forced_r_psf_magerr

--, forced2.i_psfflux_mag as forced_i_psf_mag
--, forced2.i_psfflux_magerr as forced_i_psf_magerr

--, forced2.z_psfflux_mag as forced_z_psf_mag
--, forced2.z_psfflux_magerr as forced_z_psf_magerr

--, forced2.y_psfflux_mag as forced_y_psf_mag
--, forced2.y_psfflux_magerr as forced_y_psf_magerr

-- columns which can be used for selection (needed)
--, forced3.g_apertureflux_10_mag as forced_g_apertureflux_10_mag
--, forced3.g_apertureflux_10_magerr as forced_g_apertureflux_10_magerr

--, forced3.r_apertureflux_10_mag as forced_r_apertureflux_10_mag
--, forced3.r_apertureflux_10_magerr as forced_r_apertureflux_10_magerr

--, forced3.i_apertureflux_10_mag as forced_i_apertureflux_10_mag
--, forced3.i_apertureflux_10_magerr as forced_i_apertureflux_10_magerr

--, forced3.z_apertureflux_10_mag as forced_z_apertureflux_10_mag
--, forced3.z_apertureflux_10_magerr as forced_z_apertureflux_10_magerr

--, forced3.y_apertureflux_10_mag as forced_y_apertureflux_10_mag
--, forced3.y_apertureflux_10_magerr as forced_y_apertureflux_10_magerr

FROM
s21a_wide.meas as meas
--LEFT JOIN s18a_dud_u2k.forced as s18a_forced using (object_id)
LEFT JOIN s21a_wide.meas2 as meas2 using (object_id)
LEFT JOIN s21a_wide.meas3 as meas3 using (object_id)
LEFT JOIN s21a_wide.patch_qa as patch_qa on meas.skymap_id=patch_qa.skymap_id
LEFT JOIN s21a_wide.forced as forced using (object_id)
LEFT JOIN s21a_wide.forced2 as forced2 using (object_id)


-- select region by region
-- region1: ra = [126.536996, 227.374634], dec = [-3.266990, 6.347560]
-- region2: ra = [197.866832, 251.622786], dec = [41.072401, 45.624517]
-- region3: ra = [328.803739, 41.067016], dec = [-7.991096, 8.372002]
WHERE
meas.tract  IN  ({$tract})                  AND
NOT meas.i_deblend_skipped                  AND
NOT meas2.i_sdsscentroid_flag               AND
NOT meas.i_pixelflags_edge                  AND
NOT meas.i_pixelflags_interpolatedcenter    AND
NOT meas.i_pixelflags_saturatedcenter       AND
NOT meas.i_pixelflags_crcenter              AND
NOT meas.i_pixelflags_bad                   AND
NOT meas.i_pixelflags_suspectcenter         AND
meas.i_detect_isprimary
--- NOT meas.i_pixelflags_clipped           AND --no need
--- NOT meas2.i_hsmshaperegauss_flag        AND --no need
--- meas2.i_hsmshaperegauss_sigma != 'NaN'  AND --no need

--- select galaxies: extendedness != 0
--- stars: extendedness == 0
--- meas.i_extendedness_value != 0
--meas.i_cmodel_mag - meas2.i_psfflux_mag <-0.15
--meas3.i_kronflux_psf_radius - meas3.i_kronflux_radius < -0.015
ORDER BY meas.object_id