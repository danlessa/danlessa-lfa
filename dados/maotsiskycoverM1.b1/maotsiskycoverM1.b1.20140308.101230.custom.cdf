CDF     
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.2-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20140308101230.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         100    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      67 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.012922       ccd_width_factor      	0.013000       center_x      236 pixels     center_y      326 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-03-08 13:49:01, using ingest-tsi-12.2-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20140308.101230.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:03:56 2017 UTC       ANDERS_armtime_timestamp      1488830636     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-03-08 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-03-08 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-03-08 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �S]�BH  �rdt�M�M@���    @���    �< ���<     =ŠE@���    @���    �< ���<     >b+�@��@    @��@    �< ���<     >���@��     @��     �< ���<     >�q�@� �    @� �    �< ���<     ?#@��    @��    �< ���<     ?7�a@�@    @�@    �< ���<     ?W��@�     @�     �< ���<     ?w�@��    @��    �< ���<     ?���@��    @��    �< ���<     ?���@�@    @�@    �< ���<     ?��n@�     @�     �< ���<     ?�zN@��    @��    �< ���<     ?�f9@�"�    @�"�    �< ���<     ?�R.@�&@    @�&@    �< ���<     ?�>-@�*     @�*     �< ���<     ?�*7@�-�    @�-�    �< ���<     @�&@�1�    @�1�    �< ���<     @�6@�5@    @�5@    �< ���<     @wJ@�9     @�9     �< ���<     @md@�<�    @�<�    �< ���<     @%c�@�@�    @�@�    �< ���<     @-Y�@�D@    @�D@    �< ���<     @5O�@�H     @�H     �< ���<     @=E�@�K�    @�K�    �<   ��<     @E<2@�O�    @�O�    �<   ��<     @M2j@�S@    @�S@    �<   ��<     @U(�@�W     @�W     �<   ��<     @]�@�Z�    @�Z�    �<   ��<     @e0@�^�    @�^�    �<   ��<     @m|@�b@    @�b@    �<   ��<     @u�@�f     @�f     �<   ��<     @|�#@�i�    @�i�    �<   ��<     @�w?@�m�    @�m�    �<   ��<     @�ro@�q@    @�q@    �<   ��<     @�m�@�u     @�u     �<   ��<     @�h�@�x�    @�x�    �<   ��<     @�d@�|�    @�|�    �<   ��<     @�_F@�@    @�@    �<   ��<     @�Z�@�     @�     �<   ��<     @�U�@��    @��    �<   ��<     @�Q @⋀    @⋀    �<   ��<     @�LC@�@    @�@    �<   ��<     @�G�@�     @�     �<   ��<     @�B�@��    @��    �<   ��<     @�>@⚀    @⚀    �<   ��<     @�9d@�@    @�@    �<   ��<     @�4�@�     @�     �<   ��<     @�0@��    @��    �<   ��<     @�+V@⩀    @⩀    �<   ��<     @�&�@�@    @�@    �<   ��<     @�"@�     @�     �<   ��<     @�[@��    @��    �<   ��<     @��@⸀    @⸀    �<   ��<     @�@�@    @�@    �<   ��<     @�u@��     @��     �<   ��<     @�
�@���    @���    �<   ��<     @�<@�ǀ    @�ǀ    �<   ��<     @��@��@    @��@    �<   ��<     @��@��     @��     �<   ��<     @��v@���    @���    �<   ��<     @���@�ր    @�ր    �<   ��<     @��R@��@    @��@    �<   ��<     @���@��     @��     �<   ��<     @��6@���    @���    �<   ��<     A ��@��    @��    �<   ��<     A�@��@    @��@    �<   ��<     A�N@��     @��     �<   ��<     A�@���    @���    �<   ��<     A��@��    @��    �<   ��<     A
�@��@    @��@    �<   ��<     A�L@��     @��     �<   ��<     A�@���    @���    �<   ��<     A��@��    @��    �<   ��<     Aܕ@�@    @�@    �<   ��<     A�Z@�     @�     �<   ��<     A� @��    @��    �<   ��<     A��@��    @��    �<   ��<     Aӯ@�@    @�@    �<   ��<     A�w@�     @�     �<   ��<     A�A@��    @��    �<   ��<     A �@�!�    @�!�    �<   ��<     A"��@�%@    @�%@    �<   ��<     A$ȥ@�)     @�)     �<   ��<     A&�s@�,�    @�,�    �<   ��<     A(�A@�0�    @�0�    �<   ��<     A*�@�4@    @�4@    �<   ��<     A,��@�8     @�8     �<   ��<     A.��@�;�    @�;�    �<   ��<     A0��@�?�    @�?�    �<   ��<     A2�Y@�C@    @�C@    �<   ��<     A4�-@�G     @�G     �<   ��<     A6�@�J�    @�J�    �<   ��<     A8��@�N�    @�N�    �<   ��<     A:��@�R@    @�R@    �<   ��<     A<��@�V     @�V     �<   ��<     A>�a@�Y�    @�Y�    �<   ��<     A@�;@�]�    @�]�    �<   ��<     AB�@�a@    @�a@    �<   ��<     AD��@�e     @�e     �<   ��<     AF��@�h�    @�h�    �<   ��<     AH��@�l�    @�l�    �<   ��<     AJ��@�p@    @�p@    �<   ��<     AL�h@�t     @�t     �<   ��<     AN�I@�w�    @�w�    �<   ��<     AP�)@�{�    @�{�    �<   ��<     AR�@�@    @�@    �<   ��<     AT��@�     @�     �<   ��<     AV��@��    @��    �<   ��<     AX��@㊀    @㊀    �<   ��<     AZ��@�@    @�@    �<   ��<     A\��@�     @�     �<   ��<     A^�h@��    @��    �<   ��<     A`�P@㙀    @㙀    �<   ��<     Ab�8@�@    @�@    �<   ��<     Ad�"@�     @�     �<   ��<     Af�@��    @��    �<   ��<     Ah�@㨀    @㨀    �<   ��<     Aj}�@�@    @�@    �<   ��<     Al{�@�     @�     �<   ��<     Any�@��    @��    �<   ��<     Apw�@㷀    @㷀    �<   ��<     Aru�@�@    @�@    �<   ��<     Ats�@�     @�     �<   ��<     Avq{@���    @���    �<   ��<     Axol@�ƀ    @�ƀ    �<   ��<     Azm_@��@    @��@    �<   ��<     A|kR@��     @��     �<   ��<     A~iE@���    @���    �<   ��<     A�3�@�Հ    @�Հ    �<   ��<     A�2�@��@    @��@    �<   ��<     A�1�@��     @��     �<   ��<     A�0�@���    @���    �<   ��<     A�/�@��    @��    �<   ��<     A�.�@��@    @��@    �<   ��<     A�-�@��     @��     �<   ��<     A�,�@���    @���    �<   ��<     A�+}@��    @��    �<   ��<     A�*{@��@    @��@    �<   ��<     A�)y@��     @��     �<   ��<     A�(w@���    @���    �<   ��<     A�'v@��    @��    �<   ��<     A�&u@�@    @�@    �<   ��<     A�%u@�
     @�
     �<   ��<     A�$u@��    @��    �<   ��<     A�#u@��    @��    �<   ��<     A�"v@�@    @�@    �<   ��<     A�!w@�     @�     �<   ��<     A� x@��    @��    �<   ��<     A�z@� �    @� �    �<   ��<     A�|@�$@    @�$@    �<   ��<     A�~@�(     @�(     �<   ��<     A��@�+�    @�+�    �<   ��<     A��@�/�    @�/�    �<   ��<     A��@�3@    @�3@    �<   ��<     A��@�7     @�7     �<   ��<     A��@�:�    @�:�    �<   ��<     A��@�>�    @�>�    �<   ��<     A��@�B@    @�B@    �<   ��<     A��@�F     @�F     �<   ��<     A��@�I�    @�I�    �<   ��<     A��@�M�    @�M�    �<   ��<     A��@�Q@    @�Q@    �<   ��<     A��@�U     @�U     �<   ��<     A��@�X�    @�X�    �<   ��<     A��@�\�    @�\�    �<   ��<     A��@�`@    @�`@    �<   ��<     A��@�d     @�d     �<   ��<     A��@�g�    @�g�    �<   ��<     A��@�k�    @�k�    �<   ��<     A�
�@�o@    @�o@    �<   ��<     A�	�@�s     @�s     �<   ��<     A��@�v�    @�v�    �<   ��<     A�@�z�    @�z�    �<   ��<     A�@�~@    @�~@    �<   ��<     A�@�     @�     �<   ��<     A�%@��    @��    �<   ��<     A�/@䉀    @䉀    �<   ��<     A�;@�@    @�@    �<   ��<     A�F@�     @�     �<   ��<     A�R@��    @��    �<   ��<     A� ^@䘀    @䘀    �<   ��<     A��j@�@    @�@    �<   ��<     A��v@�     @�     �<   ��<     A���@��    @��    �<   ��<     A���@䧀    @䧀    �<   ��<     A���@�@    @�@    �<   ��<     A���@�     @�     �<   ��<     A���@��    @��    �<   ��<     A���@䶀    @䶀    �<   ��<     A���@�@    @�@    �<   ��<     A���@�     @�     �<   ��<     A���@���    @���    �<   ��<     A��@�ŀ    @�ŀ    �<   ��<     A��@��@    @��@    �<   ��<     A��#@��     @��     �<   ��<     A��3@���    @���    �<   ��<     A��D@�Ԁ    @�Ԁ    �<   ��<     A��U@��@    @��@    �<   ��<     A��f@��     @��     �<   ��<     A��w@���    @���    �<   ��<     A��@��    @��    �<   ��<     A��@��@    @��@    �<   ��<     A��@��     @��     �<   ��<     A��@���    @���    �<   ��<     A���@��    @��    �<   ��<     A���@��@    @��@    �<   ��<     A���@��     @��     �<   ��<     A��
@���    @���    �<   ��<     A��@��    @��    �<   ��<     A��2@�@    @�@    �<   ��<     A��F@�	     @�	     �<   ��<     A��Z@��    @��    �<   ��<     A��o@��    @��    �<   ��<     A��@�     @�     �<   ��<     A�߮@��    @��    �<   ��<     A���@��    @��    �<   ��<     A���@�#@    @�#@    �<   ��<     A���@�'     @�'     �<   ��<     A��@�*�    @�*�    �<   ��<     A��@�.�    @�.�    �<   ��<     A��2@�2@    @�2@    �<   ��<     A��I@�6     @�6     �<   ��<     A��`@�9�    @�9�    �<   ��<     A��x@�=�    @�=�    �<   ��<     A�֏@�A@    @�A@    �<   ��<     A�է@�E     @�E     �<   ��<     A�Կ@�H�    @�H�    �<   ��<     A���@�L�    @�L�    �<   ��<     A���@�P@    @�P@    �<   ��<     A��	@�T     @�T     �<   ��<     A��!@�W�    @�W�    �<   ��<     A��;@�[�    @�[�    �<   ��<     A��T@�_@    @�_@    �<   ��<     A��m@�c     @�c     �<  ��<     A�͇@�f�    @�f�    �<   ��<     A�̡@�j�    @�j�    �<  ��<     A�˻@�n@    @�n@    �<   ��<     A���@�r     @�r     �<   ��<     A���@�u�    @�u�    �<   ��<     A��@�y�    @�y�    �<   ��<     A��&@�}@    @�}@    �<   ��<     A��A@�     @�     �<   ��<     A��\@��    @��    �<   ��<     A��x@刀    @刀    �<   ��<     A�ē@�@    @�@    �<   ��<     A�ï@�     @�     �<   ��<     A���@��    @��    �<   ��<     A���@�     @�     �<   ��<     B\@���    @���    �<   ��<     B�a@��    @��    �<   ��<     BU�@��@    @��@    �<   ��<     BՁ@��     @��     �<   ��<     B	U@���    @���    �<   ��<     B	Ԣ@�)�    @�)�    �<   ��<     B�l@�-�    @�-�    �<   ��<     BN�@�1@    @�1@    �<  ��<     BΏ@�5     @�5     �<   ��<     BN @�8�    @�8�    �<   ��<     Bͱ@�<�    @�<�    �<   ��<     BMC@�@@    @�@@    �<   ��<     B��@�D     @�D     �<   ��<     BLf@�G�    @�G�    �<   ��<     B��@�K�    @�K�    �<   ��<     BK�@�O@    @�O@    �<   ��<     B�@�S     @�S     �<   ��<     BJ�@�V�    @�V�    �<   ��<     B�?@�Z�    @�Z�    �<   ��<     BI�@�^@    @�^@    �<   ��<     B�b@�b     @�b     �<   ��<     BH�@�e�    @�e�    �<   ��<     Bȇ@�i�    @�i�    �<   ��<     BH@�m@    @�m@    �<   ��<     Bǫ@�q     @�q     �<   ��<     BG=@�t�    @�t�    �<   ��<     B��@�x�    @�x�    �<   ��<     BFa@�|@    @�|@    �<   ��<     B��@�     @�     �<   ��<     BE�@��    @��    �<   ��<     B�@懀    @懀    �<   ��<     BD�@�@    @�@    �<   ��<     B�=@�     @�     �<   ��<     BC�@��    @��    �<   ��<     B�b@斀    @斀    �<   ��<     BB�@�@    @�@    �<   ��<     B@�     @�     �<   ��<     BB@��    @��    �<   ��<     B��@楀    @楀    �<   ��<     B A?@�@    @�@    �<   ��<     B ��@�     @�     �<   ��<     B!@e@��    @��    �<   ��<     B!��@洀    @洀    �<   ��<     B"?�@�@    @�@    �<   ��<     B"�@�     @�     �<   ��<     B#>�@��    @��    �<   ��<     B#�C@�À    @�À    �<   ��<     B$=�@��@    @��@    �<   ��<     B$�i@��     @��     �<   ��<     B%<�@���    @���    �<   ��<     B%��@�Ҁ    @�Ҁ    �<   ��<     B&<"@��     @��     �<   ��<     B';H@���    @���    �<   ��<     B'��@��    @��    �<   ��<     B(:n@��@    @��@    �<   ��<     B(�@��     @��     �<   ��<     B)9�@���    @���    �<   ��<     B)�'@���    @���    �<   ��<     B*8�@��@    @��@    �<   ��<     B*�M@��     @��     �<   ��<     B+7�@���    @���    �<   ��<     B+�s@���    @���    �<   ��<     B,7@�@    @�@    �<   ��<     B,��@�     @�     �<   ��<     B-6-@�
�    @�
�    �<   ��<     B-��@��    @��    �<   ��<     B.5S@�@    @�@    �<   ��<     B.��@�     @�     �<   ��<     B/4y@��    @��    �<   ��<     B/�@��    @��    �<   ��<     B03�@�!@    @�!@    �<   ��<     B0�2@�%     @�%     �<   ��<     B12�@�(�    @�(�    �<   ��<     B1�Y@�,�    @�,�    �<   ��<     B21�@�0@    @�0@    �<   ��<     B2�@�4     @�4     �<   ��<     B31@�7�    @�7�    �<   ��<     B3��@�;�    @�;�    �<   ��<     B408@�?@    @�?@    �<   ��<     B4��@�C     @�C     �<   ��<     B5/^@�F�    @�F�    �<   ��<     B5��@�J�    @�J�    �<   ��<     B6.�@�N@    @�N@    �<   ��<     B6�@�R     @�R     �<   ��<     B7-�@�U�    @�U�    �<   ��<     B7�=@�Y�    @�Y�    �<   ��<     B8,�@�]@    @�]@    �<   ��<     B8�c@�a     @�a     �<   ��<     B9+�@�d�    @�d�    �<   ��<     B9��@�h�    @�h�    �<   ��<     B:+@�l@    @�l@    �<   ��<     B:��@�p     @�p     �<   ��<     B;*A@�s�    @�s�    �<   ��<     B;��@�w�    @�w�    �<   ��<     B<)f@�{@    @�{@    �<   ��<     B<��@�     @�     �<   ��<     B=(�@��    @��    �<   ��<     B=�@熀    @熀    �<   ��<     B>'�@�@    @�@    �<   ��<     B>�C@�     @�     �<   ��<     B?&�@��    @��    �<   ��<     B?�h@畀    @畀    �<   ��<     B@%�@�@    @�@    �<   ��<     B@��@�     @�     �<   ��<     BA%@��    @��    �<   ��<     BA��@礀    @礀    �<   ��<     BB$D@�@    @�@    �<   ��<     BB��@�     @�     �<   ��<     BC#h@��    @��    �<   ��<     BC��@糀    @糀    �<   ��<     BD"�@�@    @�@    �<   ��<     BD�@�     @�     �<   ��<     BE!�@��    @��    �<   ��<     BE�B@�    @�    �<   ��<     BF �@��@    @��@    �<   ��<     BF�e@��     @��     �<   ��<     BG�@���    @���    �<   ��<     BG��@�р    @�р    �<   ��<     BH@��@    @��@    �<   ��<     BH��@��     @��     �<   ��<     BI<@���    @���    �<   ��<     BI��@���    @���    �<   ��<     BJ_@��@    @��@    �<   ��<     BJ��@��     @��     �<   ��<     BK�@���    @���    �<   ��<     BK�@��    @��    �<   ��<     BL�@��@    @��@    �<   ��<     BL�4@��     @��     �<   ��<     BM�@���    @���    �<   ��<     BM�U@���    @���    �<   ��<     BN�@�@    @�@    �<   ��<     BN�v@�     @�     �<   ��<     BO@�	�    @�	�    �<   ��<     BO��@��    @��    �<   ��<     BP'@�@    @�@    �<   ��<     BP��@�     @�     �<   ��<     BQF@��    @��    �<   ��<     BQ��@��    @��    �<   ��<     BRf@� @    @� @    �<   ��<     BR��@�$     @�$     �<   ��<     BS�@�'�    @�'�    �<   ��<     BS�@�+�    @�+�    �<   ��<     BT�@�/@    @�/@    �<   ��<     BT�3@�3     @�3     �<   ��<     BU�@�6�    @�6�    �<   ��<     BU�Q@�:�    @�:�    �<   ��<     BV�@�>@    @�>@    �<   ��<     BV�n@�B     @�B     �<   ��<     BW�@�E�    @�E�    �<   ��<     BW��@�I�    @�I�    �<   ��<     BX@�M@    @�M@    �<   ��<     BX��@�Q     @�Q     �<   ��<     BY6@�T�    @�T�    �<   ��<     BY��@�X�    @�X�    �<   ��<     BZQ@�\@    @�\@    �<   ��<     BZ��@�`     @�`     �<   ��<     B[l@�c�    @�c�    �<   ��<     B[��@�g�    @�g�    �<   ��<     B\�@�k@    @�k@    �<   ��<     B\�@�o     @�o     �<   ��<     B]�@�r�    @�r�    �<   ��<     B]�-@�v�    @�v�    �<   ��<     B^�@�z@    @�z@    �<   ��<     B^�F@�~     @�~     �<   ��<     B_
�@��    @��    �<   ��<     B_�^@腀    @腀    �<   ��<     B`	�@�@    @�@    �<   ��<     B`�v@�     @�     �<   ��<     Ba	@��    @��    �<   ��<     Ba��@蔀    @蔀    �<   ��<     Bb@�@    @�@    �<   ��<     Bb��@�     @�     �<   ��<     Bc-@��    @��    �<   ��<     Bc��@裀    @裀    �<   ��<     BdB@�@    @�@    �<   ��<     Bd��@�     @�     �<   ��<     BeV@��    @��    �<   ��<     Be��@貀    @貀    �<   ��<     Bfi@�@    @�@    �<   ��<     Bf��@�     @�     �<   ��<     Bg{@��    @��    �<   ��<     Bg�@���    @���    �<   ��<     Bh�@��@    @��@    �<   ��<     Bh�@��     @��     �<   ��<     Bi�@���    @���    �<   ��<     Bi�%@�Ѐ    @�Ѐ    �<   ��<     Bj �@��@    @��@    �<   ��<     Bj�4@��     @��     �<   ��<     Bj��@���    @���    �<   ��<     BkB@�߀    @�߀    �<   ��<     Bk��@��@    @��@    �<   ��<     Bl~O@��     @��     �<   ��<     Bl��@���    @���    �<   ��<     Bm}[@��    @��    �<   ��<     Bm��@��@    @��@    �<   ��<     Bn|f@��     @��     �<   ��<     Bn��@���    @���    �<   ��<     Bo{p@���    @���    �<   ��<     Bo��@�@    @�@    �<   ��<     Bpzx@�     @�     �<   ��<     Bp��@��    @��    �<   ��<     Bqy�@��    @��    �<   ��<     Bq�@�@    @�@    �<   ��<     Brx�@�     @�     �<   ��<     Br�	@��    @��    �<   ��<     Bsw�@��    @��    �<   ��<     Bs�@�@    @�@    �<   ��<     Btv�@�#     @�#     �<   ��<     Bt�@�&�    @�&�    �<   ��<     Buu�@�*�    @�*�    �<   ��<     Bu�@�.@    @�.@    �<   ��<     Bvt�@�2     @�2     �<   ��<     Bv�@�5�    @�5�    �<   ��<     Bws�@�9�    @�9�    �<   ��<     Bw�@�=@    @�=@    �<   ��<     Bxr�@�A     @�A     �<   ��<     Bx�@�D�    @�D�    �<   ��<     Byq�@�H�    @�H�    �<   ��<     By�@�L@    @�L@    �<   ��<     Bzp�@�P     @�P     �<   ��<     Bz�@�S�    @�S�    �<   ��<     B{o~@�W�    @�W�    �<   ��<     B{��@�[@    @�[@    �<   ��<     B|nv@�_     @�_     �<   ��<     B|��@�b�    @�b�    �<   ��<     B}mk@�f�    @�f�    �<   ��<     B}��@�j@    @�j@    �<   ��<     B~l_@�n     @�n     �<   ��<     B~��@�q�    @�q�    �<   ��<     BkQ@�u�    @�u�    �<   ��<     B��@�y@    @�y@    �<   ��<     B�5 @�}     @�}     �<   ��<     B�t�@��    @��    �<   ��<     B���@鄀    @鄀    �<   ��<     B��S@�@    @�@    �<   ��<     B�4@�     @�     �<   ��<     B�s�@��    @��    �<   ��<     B���@铀    @铀    �<   ��<     B��<@�@    @�@    �<   ��<     B�2�@�     @�     �<   ��<     B�r�@��    @��    �<   ��<     B��i@颀    @颀    �<   ��<     B��"@�@    @�@    �<   ��<     B�1�@�     @�     �<   ��<     B�q�@��    @��    �<   ��<     B��J@鱀    @鱀    �<   ��<     B��@�@    @�@    �<   ��<     B�0�@�     @�     �<   ��<     B�pp@��    @��    �<   ��<     B��'@���    @���    �<   ��<     B���@��@    @��@    �<   ��<     B�/�@��     @��     �<   ��<     B�oI@���    @���    �<   ��<     B���@�π    @�π    �<   ��<     B��@��@    @��@    �<   ��<     B�.h@��     @��     �<   ��<     B�n@���    @���    �<   ��<     B���@�ހ    @�ހ    �<   ��<     B��@��@    @��@    �<   ��<     B�-6@��     @��     �<   ��<     B�l�@���    @���    �<   ��<     B���@��    @��    �<   ��<     B��M@��@    @��@    �<   ��<     B�+�@��     @��     �<   ��<     B�k�@���    @���    �<   ��<     B��`@���    @���    �<   ��<     B��@� @    @� @    �<   ��<     B�*�@�     @�     �<   ��<     B�jo@��    @��    �<   ��<     B��@��    @��    �<   ��<     B���@�@    @�@    �<   ��<     B�)y@�     @�     �<   ��<     B�i'@��    @��    �<   ��<     B���@��    @��    �<   ��<     B��@�@    @�@    �<   ��<     B�(+@�"     @�"     �<   ��<     B�g�@�%�    @�%�    �<   ��<     B���@�)�    @�)�    �<   ��<     B��+@�-@    @�-@    �<   ��<     B�&�@�1     @�1     �<   ��<     B�f}@�4�    @�4�    �<   ��<     B��%@�8�    @�8�    �<   ��<     B���@�<@    @�<@    �<   ��<     B�%t@�@     @�@     �<   ��<     B�e@�C�    @�C�    �<   ��<     B���@�G�    @�G�    �<   ��<     B��e@�K@    @�K@    �<   ��<     B�$
@�O     @�O     �<   ��<     B�c�@�R�    @�R�    �<   ��<     B��P@�V�    @�V�    �<   ��<     B���@�^     @�^     �<   ��<     B�b5@�a�    @�a�    �<   ��<     B���@�e�    @�e�    �<   ��<     B��u@�i@    @�i@    �<   ��<     B�!@�m     @�m     �<   ��<     B�`�@�p�    @�p�    �<   ��<     B��N@�t�    @�t�    �<   ��<     B���@�x@    @�x@    �<   ��<     B��@�|     @�|     �<   ��<     B�_@��    @��    �<   ��<     B���@ꃀ    @ꃀ    �<   ��<     B��Q@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�]�@��    @��    �<   ��<     B��@ꒀ    @ꒀ    �<   ��<     B�ܪ@�@    @�@    �<   ��<     B�>@�     @�     �<   ��<     B�[�@��    @��    �<   ��<     B��b@ꡀ    @ꡀ    �<   ��<     B���@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�Z@��    @��    �<   ��<     B���@가    @가    �<   ��<     B��(@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�X;@��    @��    �<   ��<     B���@꿀    @꿀    �<   ��<     B��J@��@    @��@    �<   ��<     B��@��     @��     �<   ��<     B�VS@���    @���    �<   ��<     B���@�΀    @�΀    �<   ��<     B��V@��@    @��@    �<   ��<     B��@��     @��     �<   ��<     B�TS@���    @���    �<   ��<     B���@�݀    @�݀    �<   ��<     B��J@��@    @��@    �<   ��<     B��@��     @��     �<   ��<     B�R9@���    @���    �<   ��<     B���@��    @��    �<   ��<     B��"@��@    @��@    �<   ��<     B��@��     @��     �<   ��<     B�P@���    @���    �<   ��<     B��q@���    @���    �<   ��<     B���@��@    @��@    �<   ��<     B�F@�     @�     �<   ��<     B�M�@��    @��    �<   ��<     B��@�
�    @�
�    �<   ��<     B��v@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�K4@��    @��    �<   ��<     B���@��    @��    �<   ��<     B���@�@    @�@    �<   ��<     B�	?@�!     @�!     �<   ��<     B�H�@�$�    @�$�    �<   ��<     B���@�(�    @�(�    �<   ��<     B��1@�,@    @�,@    �<   ��<     B�|@�0     @�0     �<   ��<     B�E�@�3�    @�3�    �<   ��<     B��@�7�    @�7�    �<   ��<     B��I@�;@    @�;@    �<   ��<     B��@�?     @�?     �<   ��<     B�B�@�B�    @�B�    �<   ��<     B���@�F�    @�F�    �<   ��<     B��)@�J@    @�J@    �<   ��<     B� X@�N     @�N     �<   ��<     B�?�@�Q�    @�Q�    �<   ��<     B�~�@�U�    @�U�    �<   ��<     B���@�Y@    @�Y@    �<   ��<     B���@�]     @�]     �<   ��<     B�< @�`�    @�`�    �<   ��<     B�{@�d�    @�d�    �<   ��<     B��#@�h@    @�h@    �<   ��<     B��,@�l     @�l     �<   ��<     B�80@�o�    @�o�    �<   ��<     B�w.@�s�    @�s�    �<   ��<     B��&@�w@    @�w@    �<   ��<     B��@�{     @�{     �<   ��<     B�4@�~�    @�~�    �<   ��<     B�r�@낀    @낀    �<   ��<     B���@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�/j@��    @��    �<   ��<     B�n1@둀    @둀    �<   ��<     B���@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�*Q@��    @��    �<   ��<     B�h�@렀    @렀    �<   ��<     B���@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B�$�@��    @��    �<   ��<     B�c@므    @므    �<   ��<     B���@�@    @�@    �<   ��<     B���@�     @�     �<   ��<     B�1@��    @��    �<   ��<     B�\t@뾀    @뾀    �<   ��<     B���@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B��@���    @���    �<   ��<     B�T�@�̀    @�̀    �<   ��<     B���@��@    @��@    �<   ��<     B�Ь@��     @��     �<   ��<     B�r@���    @���    �<   ��<     B�L"@�܀    @�܀    �<   ��<     B���@��@    @��@    �<   ��<     B��:@��     @��     �<   ��<     B��@���    @���    �<   ��<     B�A�@��    @��    �<   ��<     B�@��@    @��@    �<   ��<     B��@��     @��     �<   ��<     B���@���    @���    �<   ��<     B�5�@���    @���    �<   ��<     B�rP@��@    @��@    �<   ��<     B���@�     @�     �<   ��<     B���@��    @��    �<   ��<     B�' @�	�    @�	�    �<   ��<     B�b�@�@    @�@    �<   ��<     B��i@�     @�     �<   ��<     B�ٿ@��    @��    �<   ��<     B��@��    @��    �<   ��<     B�O�@�@    @�@    �<   ��<     B��@�      @�      �<   ��<     B��@�#�    @�#�    �<   ��<     B���@�'�    @�'�    �<   ��<     B�7@�+@    @�+@    �<   ��<     B�o�@�/     @�/     �<   ��<     B��@�2�    @�2�    �<   ��<     B�߷@�6�    @�6�    �<   ��<     B��@�:@    @�:@    �<   ��<     B�M
@�>     @�>     �<   ��<     B���@�A�    @�A�    �<   ��<     B��$@�E�    @�E�    �<   ��<     B��@�I@    @�I@    �<   ��<     B�2@�M     @�M     �<   ��<     B�Na@�P�    @�P�    �<   ��<     B�~!@�T�    @�T�    �<   ��<     B��@@�X@    @�X@    �<   ��<     B�؊@�\     @�\     �<   ��<     B��@�_�    @�_�    �<   ��<     B�*�@�c�    @�c�    �<   ��<     B�O�@�g@    @�g@    �<   ��<     B�r"@�k     @�k     �<   ��<     B��@�n�    @�n�    �<   ��<     B��d@�r�    @�r�    �<   ��<     B�æ@�v@    @�v@    �<   ��<     B�և@�z     @�z     �<   ��<     B��@�}�    @�}�    �<   ��<     B���@쁀    @쁀    �<   ��<     B��@�@    @�@    �<   ��<     B���@�     @�     �<   ��<     B��@��    @��    �<   ��<     B��d@쐀    @쐀    �<   ��<     B��7@�@    @�@    �<   ��<     B��y@�     @�     �<   ��<     B��}@��    @��    �<   ��<     B���@쟀    @쟀    �<   ��<     B�d<@�@    @�@    �<   ��<     B�@�@�     @�     �<   ��<     B�Q@��    @��    �<   ��<     B��o@쮀    @쮀    �<   ��<     B��S@�@    @�@    �<   ��<     B��>@�     @�     �<   ��<     B�jk@��    @��    �<   ��<     B�:@콀    @콀    �<   ��<     B�P@��@    @��@    �<   ��<     B��]@��     @��     �<   ��<     B��V@���    @���    �<   ��<     B�lZ@�̀    @�̀    �<   ��<     B�6�@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B�Ș@���    @���    �<   ��<     B���@�ۀ    @�ۀ    �<   ��<     B�X2@��@    @��@    �<   ��<     B�7@��     @��     �<   ��<     B���@���    @���    �<   ��<     B���@��    @��    �<   ��<     B�q�@��@    @��@    �<   ��<     B�7@��     @��     �<   ��<     B��:@���    @���    �<   ��<     B��@���    @���    �<   ��<     B���@��@    @��@    �<   ��<     B�I�@�     @�     �<   ��<     B� @��    @��    �<   ��<     B���@��    @��    �<   ��<     B���@�@    @�@    �<   ��<     B�Y@�     @�     �<   ��<     B�k@��    @��    �<   ��<     B�ߜ@��    @��    �<   ��<     B���@�@    @�@    �<   ��<     B�e�@�     @�     �<   ��<     B�(\@�"�    @�"�    �<   ��<     B��@�&�    @�&�    �<   ��<     B���@�*@    @�*@    �<   ��<     B�p@�.     @�.     �<   ��<     B�2l@�1�    @�1�    �<   ��<     B���@�5�    @�5�    �<   ��<     B���@�9@    @�9@    �<   ��<     B�x�@�=     @�=     �<   ��<     B�;@�@�    @�@�    �<   ��<     B���@�D�    @�D�    �<   ��<     B���@�H@    @�H@    �<   ��<     B���@�L     @�L     �<   ��<     B�B}@�O�    @�O�    �<   ��<     B�4@�S�    @�S�    �<   ��<     B���@�W@    @�W@    �<   ��<     B��x@�[     @�[     �<   ��<     B�I@�^�    @�^�    �<   ��<     B�
�@�b�    @�b�    �<   ��<     B��@�f@    @�f@    �<   ��<     B��o@�j     @�j     �<   ��<     B�N�@�m�    @�m�    �<   ��<     B�*@�q�    @�q�    �<   ��<     B��y@�u@    @�u@    �<   ��<     B���@�y     @�y     �<   ��<     B�S�@�|�    @�|�    �<   ��<     B�3@퀀    @퀀    �<   ��<     B��a@�@    @�@    �<   ��<     B���@�     @�     �<   ��<     B�X�@��    @��    �<   ��<     B��@폀    @폀    �<   ��<     B���@�@    @�@    �<   ��<     B���@�     @�     �<   ��<     B�\�@��    @��    �<   ��<     B��@힀    @힀    �<   ��<     B���@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B�`�@���    @���    �<   ��<     B�!�@���    @���    �<   ��<     B��@��@    @��@    �<   ��<     B��k@��     @��     �<   ��<     B�dG@���    @���    �<   ��<     B�%@���    @���    �<   ��<     B���@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B�g�@���    @���    �<   ��<     B�(T@�ˀ    @�ˀ    �<   ��<     B��@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B�j�@���    @���    �<   ��<     B�+O@�ڀ    @�ڀ    �<   ��<     B��@��@    @��@    �<   ��<     B���