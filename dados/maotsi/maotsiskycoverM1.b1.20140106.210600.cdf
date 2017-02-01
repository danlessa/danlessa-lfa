CDF      
      time          E   command_line      tsi_ingest -s mao -f M1 -R -d      process_version       ingest-tsi-12.2-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20140106210600.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         100    band_hw       
30 pixels      band_top      -30 pixels     
brightness        7      cam_mir_dist      64 mm      camera_arm_width      
20 pixels      camera_head_height        120 pixels     camera_head_width         
55 pixels      ccd_height_factor         	0.013547       ccd_width_factor      	0.013521       center_x      235 pixels     center_y      332 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
42.434132      reference_fitCoef2        
-0.137740      reference_fitCoefMag2         -70.070396     reference_fitCoef3        	0.207757       reference_fitCoefMag3         235.775101     reference_fitCoef4        
-0.217540      reference_fitCoefMag4         -207.707993    reference_magnitudeFactor         	1.000000       reference_magnitudeMaxLimit       400.000000     reference_factor      	1.000000       reference_maxLimit        	0.400000       region_horizon_az         50 degrees     region_horizon_alt        40 degrees     region_horizon_enabled        true       region_sun_enabled        true       region_sun_radius         25 degrees     region_zenith_enabled         true       region_zenith_radius      50 degrees     opaque_thresh         85     sunny_thresh      35     thin_thresh       45     	site_name         mao    latitude      -3.21297 degrees       	longitude         -60.5981 degrees       qc_standards_version      1.0    	qc_method         Standard Mentor QC     
qc_comment       The QC field values are a bit packed representation of true/false values for the tests that may have been performed. A QC value of zero means that none of the tests performed on the value failed.

The QC field values make use of the internal binary format to store the results of the individual QC tests. This allows the representation of multiple QC states in a single value. If the test associated with a particular bit fails the bit is turned on. Turning on the bit equates to adding the integer value of the failed test to the current value of the field. The QC field's value can be interpreted by applying bit logic using bitwise operators, or by examining the QC value's integer representation. A QC field's integer representation is the sum of the individual integer values of the failed tests. The bit and integer equivalents for the first 5 bits are listed below:

bit_1 = 00000001 = 0x01 = 2^0 = 1
bit_2 = 00000010 = 0x02 = 2^1 = 2
bit_3 = 00000100 = 0x04 = 2^2 = 4
bit_4 = 00001000 = 0x08 = 2^3 = 8
bit_5 = 00010000 = 0x10 = 2^4 = 16       qc_bit_1_description      !Value is equal to missing_value.       qc_bit_1_assessment       Bad    qc_bit_2_description      "Value is less than the valid_min.      qc_bit_2_assessment       Bad    qc_bit_3_description      %Value is greater than the valid_max.       qc_bit_3_assessment       Bad    
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-01-08 22:39:14, using ingest-tsi-12.2-0.el6          5   	base_time                string        2014-01-06 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         Ix   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-01-06 00:00:00 0:00          I�   time                	long_name         Time offset from midnight      units         'seconds since 2014-01-06 00:00:00 0:00          I�   qc_time                 	long_name         :Quality check results on field: Time offset from midnight      units         	unitless       description       vThis field contains bit packed values which should be interpreted as listed. No bits set (zero) represents good data.      bit_1_description         9Delta time between current and previous samples is zero.       bit_1_assessment      Indeterminate      bit_2_description         fDelta time between current and previous samples is less than the delta_t_lower_limit field attribute.      bit_2_assessment      Indeterminate      bit_3_description         iDelta time between current and previous samples is greater than the delta_t_upper_limit field attribute.       bit_3_assessment      Indeterminate      delta_t_lower_limit       @=         delta_t_upper_limit       @?         prior_sample_flag                comment       �If the 'prior_sample_flag' is set the first sample time from a new raw file will be compared against the time just previous to it in the stored data. If it is not set the qc_time value for the first sample will be set to 0.         I�   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         I�   qc_percent_opaque                   	long_name         5Quality check results on field: Percent opaque cloud       units         	unitless       description       7See global attributes for individual bit descriptions.          I�   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         I�   qc_percent_thin                 	long_name         3Quality check results on field: Percent thin cloud     units         	unitless       description       7See global attributes for individual bit descriptions.          I�   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          I�   qc_sunny                	long_name         /Quality check results on field: Sunshine meter     units         	unitless       description       7See global attributes for individual bit descriptions.          I�   sun_strength                	long_name         "Relative 'strength' of direct sun      units         	unitless       	valid_min         ��     	valid_max         B�     
resolution               missing_value         �<         I�   qc_sun_strength                 	long_name         BQuality check results on field: Relative 'strength' of direct sun      units         	unitless       description       7See global attributes for individual bit descriptions.          I�   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         I�   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          I�   solar_azimuth                   	long_name         Solar azimuth angle    units         degree     	valid_min                	valid_max         C�     
resolution        ?�     missing_value         �<         I�   qc_solar_azimuth                	long_name         4Quality check results on field: Solar azimuth angle    units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_zenith_count_thin                	long_name         *Pixel count: number thin in zenith circle      units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_zenith_count_thin                 	long_name         JQuality check results on field: Pixel count: number thin in zenith circle      units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_zenith_count_opaque                  	long_name         ,Pixel count: number opaque in zenith circle    units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_zenith_count_opaque                   	long_name         LQuality check results on field: Pixel count: number opaque in zenith circle    units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_zenith_count                 	long_name         +Pixel count: number total in zenith circle     units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_zenith_count                  	long_name         KQuality check results on field: Pixel count: number total in zenith circle     units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_sun_count_thin                   	long_name         'Pixel count: number thin in sun circle     units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_sun_count_thin                	long_name         GQuality check results on field: Pixel count: number thin in sun circle     units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_sun_count_opaque                 	long_name         )Pixel count: number opaque in sun circle       units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_sun_count_opaque                  	long_name         IQuality check results on field: Pixel count: number opaque in sun circle       units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_sun_count                	long_name         (Pixel count: number total in sun circle    units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_sun_count                 	long_name         HQuality check results on field: Pixel count: number total in sun circle    units         	unitless       description       7See global attributes for individual bit descriptions.          I�   region_horizon_count_thin                   	long_name         )Pixel count: number thin in horizon area       units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         I�   qc_region_horizon_count_thin                	long_name         IQuality check results on field: Pixel count: number thin in horizon area       units         	unitless       description       7See global attributes for individual bit descriptions.          J    region_horizon_count_opaque                 	long_name         +Pixel count: number opaque in horizon area     units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         J   qc_region_horizon_count_opaque                  	long_name         KQuality check results on field: Pixel count: number opaque in horizon area     units         	unitless       description       7See global attributes for individual bit descriptions.          J   region_horizon_count                	long_name         *Pixel count: number total in horizon area      units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         J   qc_region_horizon_count                 	long_name         JQuality check results on field: Pixel count: number total in horizon area      units         	unitless       description       7See global attributes for individual bit descriptions.          J   count_sub_proczen                   	long_name         ?Pixel count: number total between horizon and processed circle     units         pixels     	valid_min         ��     	valid_max         H�     
resolution        ?�     missing_value         �<         J   qc_count_sub_proczen                	long_name         _Quality check results on field: Pixel count: number total between horizon and processed circle     units         	unitless       description       7See global attributes for individual bit descriptions.          J   count_opaque                	long_name         !Pixel count: number total opaque       units         pixels     	valid_min                	valid_max          �    
resolution              missing_value         ����        J   qc_count_opaque                 	long_name         AQuality check results on field: Pixel count: number total opaque       units         	unitless       description       7See global attributes for individual bit descriptions.          J    
count_thin                  	long_name         Pixel count: number total thin     units         pixels     	valid_min                	valid_max          �    
resolution              missing_value         ����        J$   qc_count_thin                   	long_name         ?Quality check results on field: Pixel count: number total thin     units         	unitless       description       7See global attributes for individual bit descriptions.          J(   	count_box                   	long_name         0Pixel count: number in box, outside mirror area    units         pixels     	valid_min         ��     	valid_max         H�     
resolution        ?�     missing_value         �<         J,   qc_count_box                	long_name         PQuality check results on field: Pixel count: number in box, outside mirror area    units         	unitless       description       7See global attributes for individual bit descriptions.          J0   	count_sky                   	long_name         .Pixel count: number total in processed circle      units         pixels     	valid_min         ��     	valid_max         H�     
resolution        ?�     missing_value         �<         J4   qc_count_sky                	long_name         NQuality check results on field: Pixel count: number total in processed circle      units         	unitless       description       7See global attributes for individual bit descriptions.          J8   count_unknown                   	long_name         (Pixel count: number total indeterminate    units         pixels     	valid_min                	valid_max         H�     
resolution        ?�     missing_value         �<         J<   qc_count_unknown                	long_name         HQuality check results on field: Pixel count: number total indeterminate    units         	unitless       description       7See global attributes for individual bit descriptions.          J@   
count_mask                  	long_name         1Pixel count: number in camera and sun strip mask       units         pixels     	valid_min         ��     	valid_max         H�     
resolution        ?�     missing_value         �<         JD   qc_count_mask                   	long_name         QQuality check results on field: Pixel count: number in camera and sun strip mask       units         	unitless       description       7See global attributes for individual bit descriptions.          JH   count_sub_horz                  	long_name         +Pixel count: number below horizon in image     units         pixels     	valid_min         ��     	valid_max         H�     
resolution        ?�     missing_value         �<         JL   qc_count_sub_horz                   	long_name         KQuality check results on field: Pixel count: number below horizon in image     units         	unitless       description       7See global attributes for individual bit descriptions.          JP   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            I|   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           I�   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            I�R�� �M�M�rdtBH  @�    @�        B�                �    �<    Az.(    Cw��            GAv     GAv             E�P     E�P             FJ     FJ     F݀      �p            Gvv     G��             F�     G{�     