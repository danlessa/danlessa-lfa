CDF   i   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.3-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20141004210000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      69 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.012922       ccd_width_factor      	0.013000       center_x      239 pixels     center_y      317 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
42.434132      reference_fitCoef2        
-0.137740      reference_fitCoefMag2         -70.070396     reference_fitCoef3        	0.207757       reference_fitCoefMag3         235.775101     reference_fitCoef4        
-0.217540      reference_fitCoefMag4         -207.707993    reference_magnitudeFactor         	1.000000       reference_magnitudeMaxLimit       400.000000     reference_factor      	1.000000       reference_maxLimit        	0.400000       region_horizon_az         50 degrees     region_horizon_alt        40 degrees     region_horizon_enabled        true       region_sun_enabled        true       region_sun_radius         25 degrees     region_zenith_enabled         true       region_zenith_radius      50 degrees     opaque_thresh         85     sunny_thresh      35     thin_thresh       65     	site_name         mao    latitude      -3.21297 degrees       	longitude         -60.5981 degrees       qc_standards_version      1.0    	qc_method         Standard Mentor QC     
qc_comment       The QC field values are a bit packed representation of true/false values for the tests that may have been performed. A QC value of zero means that none of the tests performed on the value failed.

The QC field values make use of the internal binary format to store the results of the individual QC tests. This allows the representation of multiple QC states in a single value. If the test associated with a particular bit fails the bit is turned on. Turning on the bit equates to adding the integer value of the failed test to the current value of the field. The QC field's value can be interpreted by applying bit logic using bitwise operators, or by examining the QC value's integer representation. A QC field's integer representation is the sum of the individual integer values of the failed tests. The bit and integer equivalents for the first 5 bits are listed below:

bit_1 = 00000001 = 0x01 = 2^0 = 1
bit_2 = 00000010 = 0x02 = 2^1 = 2
bit_3 = 00000100 = 0x04 = 2^2 = 4
bit_4 = 00001000 = 0x08 = 2^3 = 8
bit_5 = 00010000 = 0x10 = 2^4 = 16       qc_bit_1_description      !Value is equal to missing_value.       qc_bit_1_assessment       Bad    qc_bit_2_description      "Value is less than the valid_min.      qc_bit_2_assessment       Bad    qc_bit_3_description      %Value is greater than the valid_max.       qc_bit_3_assessment       Bad    
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-10-05 00:49:02, using ingest-tsi-12.3-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20141004.210000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:13 2017 UTC       ANDERS_armtime_timestamp      1488830653     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-10-04 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-10-04 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-10-04 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �T/8�BH  �rdt�M�M@�u     @�u     BaM/  �A�-�    AOC�@�v�    @�v�    B]�  �A���    AME�@�x�    @�x�    B]	�  �A�h@    AKG�@�z�    @�z�    B]�  �Aꄉ    AII�@�|�    @�|�    B[��  �A��o    AGK�@�~`    @�~`    BX/  �A�    AEM�@�@    @�@    BV��  �A�th    ACO�@�     @�     B[ /  �Aհ�    AAQ�@�     @�     B[m>  �Aɥ9    A?Sy@��    @��    B[Ր  �Aƀ�    A=Uo@��    @��    BV_�  �A�z    A;Wg@�    @�    BXy�  �A�9H    A9Y`@�    @�    BWӼ  �A���    A7[Y@�`    @�`    BS�  �A���    A5]T@�@    @�@    BO�0  �A�c�    A3_O@�     @�     BL��  �A��    A1aK@�     @�     BKs�  �A���    A/cI@��    @��    BG�#  �A�<�    A-eG@��    @��    BE  �A�T�    A+gE@�    @�    B@��  �A�v.    A)iE@�    @�    B<(  �A�PC    A'kF@�`    @�`    B7�.  �A��    A%mH@�@    @�@    B7�j  �A��    A#oJ@�     @�     B3  �A��    A!qN@�     @�     B,  �A���    AsR@��    @��    B'��  �A��    AuW@��    @��    B%:  �A�!�    Aw^@�    @�    B#SD  �A���    Aye@�    @�    B�U  �A�<    A{m@�`    @�`    B�9  �A�@�    A}v@�@    @�@    B�  �Ay�C    A�@�     @�     B��  �At�    A��@�     @�     B��  �An�     A��@��    @��    B;�  �Ais�    A��@��    @��    Bs�  �A^�	    A��@�    @�    B!#  �A^z�    A	��@�    @�    A�ם  �AMx�    A��@�`    @�`    B�v  �AW�    A��@�@    @�@    A���  �AA
�    A��@�     @�     A�<	  �A?ˌ    A�@��     @��     A��q  �A>��    @�(3@���    @���    A�[�  �A4��    @�,]@���    @���    A�h�  �A4�    @�0�@�Š    @�Š    A�4Q  �A,;S    @�4�@�ǀ    @�ǀ    A���  �A$/)    @�8�@��`    @��`    A��<  �A&HR    @�=@��@    @��@    A�PM  �A��    @�AK@��     @��     A��  �A	��    @�E�@��     @��     A�ԗ  �A�4    @�I�@���    @���    A���  �@�T    @�M�@���    @���    A��  �@���    @�R/@�Ԡ    @�Ԡ    A��  �@�NI    @�Vm@�ր    @�ր    A�	�  �@�1    @�Z�@��`    @��`    AΆ�  �@Ԟ    @�^�@��@    @��@    Aȍf  �@˻	    @�c3@��     @��     Aȫk  �@�c
    @�gz@��     @��     AÍ�  �@��<    @�k�@���    @���    A�-b  �@�G    @�p@���    @���    A��3  �@���    @�tZ@��    @��    A��  �@���    @�x�@��    @��    A�4[  �@�й    @�|�@��`    @��`    A�*/  �@�l    @��N@��@    @��@    A�5>  �@� �    @���@��     @��     A��  �@�}�    @���@��     @��     A��4  �@�.�    @��V@���    @���    A���  �@w)	    @���@���    @���    A��  �@mbn    @��@��    @��    A��]  �@b��    @��q@��    @��    A��E  �@W]<    @���@��`    @��`    A|4�  �@O�L    @��9@��@    @��@    A���  �@KX@    @���@��     @��     A`�2  �@0S�    @��@��     @��     A`Wf  �@$��    @b�@���    @���    A]J�  �@��    @wk�@���    @���    AJ_  �@
D    @ot�@��    @��    AFǕ  �@

�    @g}�@��    @��    A9i�  �@��    @_�}@�`    @�`    A3�  �?�    @W�l@�@    @�@    A�  �?ԗT    @O�a@�	     @�	     A�-  �?�l�    @G�Z@�     @�     ��  ����      @?�W@��    @��    ��  ����      @7�Z@��    @��    ��  ����      @/�a@��    @��    ��  ����      @'�m@��    @��    ��  ����      @�~@�`    @�`    ��  ����      @ה@�@    @�@    ��  ����      @�@�     @�     ��  ����      @��@�     @�     ��  ����      ?���@��    @��    ��  ����      ?��8@��    @��    ��  ����      ?�
�@��    @��    ��  ����      ?��@�!�    @�!�    ��  ����      ?�/l@�#`    @�#`    ��  ����      ?�A�@�%@    @�%@    ��  ����      ?�Tl@�'     @�'     ��  ����      ?�f�@�)     @�)     ��  ����      ?�y�@�*�    @�*�    ��  ����      ?ap@�,�    @�,�    ��  ����      ?A=�@�.�    @�.�    ��  ����      ?!c=@�0�    @�0�    ��  ����      ?��@�2`    @�2`    ��  ����      >�\�@�4@    @�4@    ��  ����      >��@�6     @�6     ��  ����      >�C@�8     @�8     ��  ����      <�