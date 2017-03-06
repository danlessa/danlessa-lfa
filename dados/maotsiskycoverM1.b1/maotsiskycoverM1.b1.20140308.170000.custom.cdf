CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.2-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20140308170000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         100    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      66 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.012922       ccd_width_factor      	0.013000       center_x      240 pixels     center_y      323 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-03-08 19:49:01, using ingest-tsi-12.2-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20140308.170000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:03:56 2017 UTC       ANDERS_armtime_timestamp      1488830636     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-03-08 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-03-08 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-03-08 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �S]�BH  �rdt�M�M@��     @��     �<   ��<     B�mi@���    @���    �<   ��<     B�.@��    @��    �<   ��<     B���@��@    @��@    �<   ��<     B��j@��     @��     �<   ��<     B�p@���    @���    �<   ��<     B�0�@���    @���    �<   ��<     B��R@��@    @��@    �<   ��<     B���@�      @�      �<   ��<     B�r�@��    @��    �<   ��<     B�3%@��    @��    �<   ��<     B��@�@    @�@    �<   ��<     B��Q@�     @�     �<   ��<     B�t�@��    @��    �<   ��<     B�5u@��    @��    �<   ��<     B��@�@    @�@    �<   ��<     B���@�     @�     �<   ��<     B�w@�!�    @�!�    �<   ��<     B�7�@�%�    @�%�    �<   ��<     B��-@�)@    @�)@    �<   ��<     B���@�-     @�-     �<   ��<     B�y7@�0�    @�0�    �<   ��<     B�9�@�4�    @�4�    �<   ��<     B��:@�8@    @�8@    �<   ��<     B���@�<     @�<     �<   ��<     B�{7@�?�    @�?�    �<   ��<     B�;�@�C�    @�C�    �<   ��<     B��/@�G@    @�G@    �<   ��<     B���@�K     @�K     �<   ��<     B�} @�N�    @�N�    �<   ��<     B�=�@�R�    @�R�    �<   ��<     B��@�V@    @�V@    �<   ��<     B���@�Z     @�Z     �<   ��<     B�~�@�]�    @�]�    �<   ��<     B�?f@�a�    @�a�    �<   ��<     B���@�e@    @�e@    �<   ��<     B��F@�i     @�i     �<   ��<     B���@�l�    @�l�    �<   ��<     B�A"@�p�    @�p�    �<   ��<     B��@�t@    @�t@    �<   ��<     B���@�x     @�x     �<   ��<     B��d@�{�    @�{�    �<   ��<     B�B�@��    @��    �<   ��<     B�5@�@    @�@    �<   ��<     B�Ü@�     @�     �<   ��<     B��@��    @��    �<   ��<     B�Dh@    @    �<   ��<     B��@�@    @�@    �<   ��<     B��0@�     @�     �<   ��<     B���@��    @��    �<   ��<     B�E�@    @    �<   ��<     B�V@�@    @�@    �<   ��<     B�ƶ@�     @�     �<   ��<     B��@��    @��    �<   ��<     B�Gt@    @    �<   ��<     B��@�@    @�@    �<   ��<     B��/@�     @�     �<   ��<     B���@��    @��    �<   ��<     B�H�@    @    �<   ��<     B�	C@�@    @�@    �<   ��<     B�ɝ@��     @��     �<   ��<     B���@���    @���    �<   ��<     B�JP@�ʀ    @�ʀ    �<   ��<     B�
�@��@    @��@    �<   ��<     B�� @��     @��     �<   ��<     B��W@���    @���    �<   ��<     B�K�@�ـ    @�ـ    �<   ��<     B�@��@    @��@    �<   ��<     B��Y@��     @��     �<   ��<     B���@���    @���    �<   ��<     B�M@��    @��    �<   ��<     B�V@��@    @��@    �<   ��<     B�ͩ@��     @��     �<   ��<     B���@���    @���    �<   ��<     B�NN@���    @���    �<   ��<     B��@��@    @��@    �<   ��<     B���@��     @��     �<   ��<     B��A@��    @��    �<   ��<     B�O�@��    @��    �<   ��<     B��@�
@    @�
@    �<   ��<     B��0@�     @�     �<   ��<     B��@��    @��    �<   ��<     B�P�@��    @��    �<   ��<     B�@�@    @�@    �<   ��<     B��i@�     @�     �<   ��<     B���@� �    @� �    �<   ��<     B�R@�$�    @�$�    �<   ��<     B�N@�(@    @�(@    �<   ��<     B�Қ@�,     @�,     �<   ��<     B���@�/�    @�/�    �<   ��<     B�S0@�3�    @�3�    �<   ��<     B�{@�7@    @�7@    �<   ��<     B���@�;     @�;     �<   ��<     B��@�>�    @�>�    �<   ��<     B�TX@�B�    @�B�    �<   ��<     B��@�F@    @�F@    �<   ��<     B���@�J     @�J     �<   ��<     B��2@�M�    @�M�    �<   ��<     B�U{@�Q�    @�Q�    �<   ��<     B��@�U@    @�U@    �<   ��<     B��
@�Y     @�Y     �<   ��<     B��Q@�\�    @�\�    �<   ��<     B�V�@�`�    @�`�    �<   ��<     B��@�d@    @�d@    �<   ��<     B��$@�h     @�h     �<   ��<     B��j@�k�    @�k�    �<   ��<     B�W�@�o�    @�o�    �<   ��<     B��@�s@    @�s@    �<   ��<     B��9@�w     @�w     �<   ��<     B��~@�z�    @�z�    �<   ��<     B�X�@�~�    @�~�    �<   ��<     B�@�@    @�@    �<   ��<     B��J@�     @�     �<   ��<     B���@��    @��    �<   ��<     B�Y�@    @    �<   ��<     B�@�@    @�@    �<   ��<     B��@�     @�     �<   ��<     B53@��    @��    �<   ��<     B~��@    @    �<   ��<     B~6<@�@    @�@    �<   ��<     B}��@�     @�     �<   ��<     B}7C@��    @��    �<   ��<     B|��@變    @變    �<   ��<     B|8H@�@    @�@    �<   ��<     B{��@�     @�     �<   ��<     B{9K@��    @��    �<   ��<     Bz��@ﺀ    @ﺀ    �<   ��<     Bz:L@�@    @�@    �<   ��<     By��@��     @��     �<   ��<     By;L@���    @���    �<   ��<     Bx��@�ɀ    @�ɀ    �<   ��<     Bx<J@��@    @��@    �<   ��<     Bw��@��     @��     �<   ��<     Bw=F@���    @���    �<   ��<     Bv��@�؀    @�؀    �<   ��<     Bv>A@��@    @��@    �<   ��<     Bu��@��     @��     �<   ��<     Bu?;@���    @���    �<   ��<     Bt��@��    @��    �<   ��<     Bt@3@��@    @��@    �<   ��<     Bs��@��     @��     �<   ��<     BsA)@���    @���    �<   ��<     Br��@���    @���    �<   ��<     BrB@��@    @��@    �<   ��<     Bq@��     @��     �<   ��<     BqC@� �    @� �    �<   ��<     BpÌ@��    @��    �<   ��<     BpD@��    @��    �<   ��<     Bo�~@��    @��    �<   ��<     BoD�@�`    @�`    �<   ��<     Bn�o@�
@    @�
@    �<   ��<     BnE�@�     @�     �<   ��<     Bm�^@�     @�     �<   ��<     BmF�@��    @��    �<   ��<     Bl�M@��    @��    �<   ��<     BlG�@��    @��    �<   ��<     Bk�:@��    @��    �<   ��<     BkH�@�`    @�`    �<   ��<     Bj�&@�@    @�@    �<   ��<     BjI�@�     @�     �<   ��<     Bi�@�     @�     �<   ��<     BiJ�@��    @��    �<   ��<     Bh��@� �    @� �    �<   ��<     BhKp@�"�    @�"�    �<   ��<     Bg��@�$�    @�$�    �<   ��<     BgLY@�&`    @�&`    �<   ��<     Bf��@�(@    @�(@    �<   ��<     BfM@@�*     @�*     �<   ��<     Beʹ@�,     @�,     �<   ��<     BeN'@�-�    @�-�    �<   ��<     BdΚ@�/�    @�/�    �<   ��<     BdO@�1�    @�1�    �<   ��<     Bcπ@�3�    @�3�    �<   ��<     BcO�@�5`    @�5`    �<   ��<     Bb�d@�7@    @�7@    �<   ��<     BbP�@�9     @�9     �<   ��<     Ba�H@�;     @�;     �<   ��<     BaQ�@�<�    @�<�    �<   ��<     B`�+@�>�    @�>�    �<   ��<     B`R�@�@�    @�@�    �<   ��<     B_�@�B�    @�B�    �<   ��<     B_S~@�D`    @�D`    �<   ��<     B^��@�F@    @�F@    �<   ��<     B^T_@�H     @�H     �<   ��<     B]��@�J     @�J     �<   ��<     B]U?@�K�    @�K�    �<   ��<     B\կ@�M�    @�M�    �<   ��<     B\V@�O�    @�O�    �<   ��<     B[֎@�Q�    @�Q�    �<   ��<     B[V�@�S`    @�S`    �<   ��<     BZ�m@�U@    @�U@    �<   ��<     BZW�@�W     @�W     �<   ��<     BY�J@�Y     @�Y     �<   ��<     BYX�@�Z�    @�Z�    �<   ��<     BX�(@�\�    @�\�    �<   ��<     BXY�@�^�    @�^�    �<   ��<     BW�@�`�    @�`�    �<   ��<     BWZr@�b`    @�b`    �<   ��<     BV��@�d@    @�d@    �<   ��<     BV[N@�f     @�f     �<   ��<     BUۼ@�h     @�h     �<   ��<     BU\)@�i�    @�i�    �<   ��<     BTܗ@�k�    @�k�    �<   ��<     BT]@�m�    @�m�    �<   ��<     BS�q@�o�    @�o�    �<   ��<     BS]�@�q`    @�q`    �<   ��<     BR�K@�s@    @�s@    �<   ��<     BR^�@�u     @�u     �<   ��<     BQ�%@�w     @�w     �<   ��<     BQ_�@�x�    @�x�    �<   ��<     BP��@�z�    @�z�    �<   ��<     BP`j@�|�    @�|�    �<   ��<     BO��@�~�    @�~�    �<   ��<     BOaB@��`    @��`    �<   ��<     BN�@��@    @��@    �<   ��<     BNb@��     @��     �<   ��<     BM�@��     @��     �<   ��<     BMb�@���    @���    �<   ��<     BL�\@���    @���    �<   ��<     BLc�@���    @���    �<   ��<     BK�3@���    @���    �<   ��<     BKd�@��`    @��`    �<   ��<     BJ�	@�@    @�@    �<   ��<     BJeu@�     @�     �<   ��<     BI��@�     @�     �<   ��<     BIfJ@��    @��    �<   ��<     BH�@��    @��    �<   ��<     BHg @�    @�    �<   ��<     BG�@�    @�    �<   ��<     BGg�@�`    @�`    �<   ��<     BF�_@�@    @�@    �<   ��<     BFh�@�     @�     �<   ��<     BE�4@�     @�     �<   ��<     BEi�@��    @��    �<   ��<     BD�@��    @��    �<   ��<     BDjr@�    @�    �<   ��<     BC��@�    @�    �<   ��<     BCkF@�`    @�`    �<   ��<     BB�@�@    @�@    �<   ��<     BBl@�     @�     �<   ��<     BA�@�u     @�u     �<   ��<     A��?@�v�    @�v�    �<   ��<     A��,@�x�    @�x�    �<   ��<     A��@�z�    @�z�    �<   ��<     A��@�|�    @�|�    �<   ��<     A���@�~`    @�~`    �<   ��<     A���@�@    @�@    �<   ��<     A���@�     @�     �<   ��<     A���@�     @�     �<   ��<     A���@��    @��    �<   ��<     A���@��    @��    �<   ��<     A���@�    @�    �<   ��<     A��@�    @�    �<   ��<     A��p@�`    @�`    �<   ��<     A��`@�@    @�@    �<   ��<     A��Q@�     @�     �<   ��<     A��C@�     @�     �<   ��<     A��5@��    @��    �<   ��<     A��'@��    @��    �<   ��<     A��@�    @�    �<   ��<     A��@�    @�    �<   ��<     Am�@�`    @�`    �<   ��<     A}o�@�@    @�@    �<   ��<     A{q�@�     @�     �<   ��<     Ays�@�     @�     �<   ��<     Awu�@��    @��    �<   ��<     Auw�@��    @��    �<   ��<     Asyq@�    @�    �<   ��<     Aq{\@�    @�    �<   ��<     Ao}G@�`    @�`    �<   ��<     Am4@�@    @�@    �<   ��<     Ak�!@�     @�     �<   ��<     Ai�@�     @�     �<   ��<     Ag��@��    @��    �<   ��<     Ae��@��    @��    �<   ��<     Ac��@�    @�    �<   ��<     Aa��@�    @�    �<   ��<     A_��@�`    @�`    �<   ��<     A]��@�@    @�@    �<   ��<     A[��@�     @�     �<   ��<     AY��@��     @��     �<   ��<     AW��@���    @���    �<   ��<     AU�@���    @���    �<   ��<     AS�t@�Š    @�Š    �<   ��<     AQ�k@�ǀ    @�ǀ    �<   ��<     AO�b@��`    @��`    �<   ��<     AM�Z@��@    @��@    �<   ��<     AK�S@��     @��     �<   ��<     AI�L@��     @��     �<   ��<     AG�G@���    @���    �<   ��<     AE�B@���    @���    �<   ��<     AC�=@�Ԡ    @�Ԡ    �<   ��<     AA�:@�ր    @�ր    �<   ��<     A?�7@��`    @��`    �<   ��<     A=�5@��@    @��@    �<   ��<     A;�4@��     @��     �<   ��<     A9�4@��     @��     �<   ��<     A7�5@���    @���    �<   ��<     A5�6@���    @���    �<   ��<     A3�8@��    @��    �<   ��<     A1�;@��    @��    �<   ��<     A/�?@��`    @��`    �<   ��<     A-�C@��@    @��@    �<   ��<     A+�I@��     @��     �<   ��<     A)�O@��     @��     �<   ��<     A'�V@���    @���    �<   ��<     A%�^@���    @���    �<   ��<     A#�g@��    @��    �<   ��<     A!�p@��    @��    �<   ��<     A�{@��`    @��`    �<   ��<     AΆ@��@    @��@    �<   ��<     AВ@��     @��     �<   ��<     Aҟ@��     @��     �<   ��<     Aԭ@���    @���    �<   ��<     Aּ@���    @���    �<   ��<     A��@��    @��    �<   ��<     A��@��    @��    �<   ��<     A��@�`    @�`    �<   ��<     A� @�@    @�@    �<   ��<     A�@�	     @�	     �<   ��<     A	�'@�     @�     �<   ��<     A�<@��    @��    �<   ��<     A�Q@��    @��    �<   ��<     A�h@��    @��    �<   ��<     A�@��    @��    �<   ��<     @��1@�`    @�`    �<   ��<     @��d@�@    @�@    �<   ��<     @��@�     @�     �<   ��<     @���@�     @�     �<   ��<     @��@��    @��    �<   ��<     @��B@��    @��    �<   ��<     @��@��    @��    �<   ��<     @���@�!�    @�!�    �<   ��<     @���@�#`    @�#`    �<   ��<     @�?@�%@    @�%@    �<   ��<     @��@�'     @�'     �<   ��<     @�	�@�)     @�)     �<   ��<     @�@�*�    @�*�    �<   ��<     @�Z@�,�    @�,�    �<   ��<     @��@�.�    @�.�    �<   ��<     @��@�0�    @�0�    �<   ��<     @�C@�2`    @�2`    �<   ��<     @�#�@�4@    @�4@    �<   ��<     @�'�@�6     @�6     �<   ��<     @�,>@�8     @�8     �<   ��<     @�0�@�9�    @�9�    �<   ��<     @�4�@�;�    @�;�    �<   ��<     @�9K@�=�    @�=�    �<   ��<     @�=�@�?�    @�?�    �<   ��<     @�B	@�A`    @�A`    �<   ��<     @�Fk@�C@    @�C@    �<   ��<     @�J�@�E     @�E     �<   ��<     @�O5@�G     @�G     �<   ��<     @�S�@�H�    @�H�    �<   ��<     @�X@�J�    @�J�    �<   ��<     @�\s@�L�    @�L�    �<   ��<     @�`�@�N�    @�N�    �<   ��<     @�eR@�P`    @�P`    �<   ��<     @xӊ@�R@    @�R@    �<   ��<     @p�s@�T     @�T     �<   ��<     @h�a@�V     @�V     �<   ��<     @`�T@�W�    @�W�    �<   ��<     @X�K@�Y�    @�Y�    �<   ��<     @Q F@�[�    @�[�    �<   ��<     @I	E@�]�    @�]�    �<   ��<     @AI@�_`    @�_`    �< ���<     @9R@�a@    @�a@    �< ���<     @1$_@�c     @�c     �< ���<     @)-p@�e     @�e     �< ���<     @!6�@�f�    @�f�    �< ���<     @?�@�h�    @�h�    �< ���<     @H�@�j�    @�j�    �< ���<     @	Q�@�l�    @�l�    �< ���<     @[@�n`    @�n`    �< ���<     ?��p@�p@    @�p@    �< ���<     ?���@�r     @�r     �< ���<     ?��?@�t     @�t     �< ���<     ?���@�u�    @�u�    �< ���<     ?�3@�w�    @�w�    �< ���<     ?�$�@�y�    @�y�    �< ���<     ?�7M@�{�    @�{�    �< ���<     ?�I�@�}`    @�}`    �< ���<     ?f�@�@    @�@    �< ���<     ?F�y@�     @�     �< ���<     ?'�@�     @�     �< ���<     ?)n@��    @��    �< ���<     >Ξ@��    @��    �< ���<     >��c@�    @�    �< ���<     >i�@�    @�    �< ���<     <�y