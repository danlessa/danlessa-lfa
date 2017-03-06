CDF   �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.4-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20151013094700.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      61 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      231 pixels     center_y      314 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Zcreated by user dsmgr on machine ruby at 2015-10-13 11:49:00, using ingest-tsi-12.4-0.el6      ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20151013.094700.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:45 2017 UTC       ANDERS_armtime_timestamp      1488830685     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2015-10-13 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2015-10-13 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2015-10-13 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �VI�BH  �rdt�M�M@�2�    @�2�    ��  ����      :��@�6@    @�6@    ��  ����      > ��@�:     @�:     ��  ����      >,�@�=�    @�=�    ��  ����      >��(@�A�    @�A�    ��  ����      >�B/@�E@    @�E@    ��  ����      ?�;@�I     @�I     ��  ����      ?>w�@�L�    @�L�    ��  ����      ?^"�@�P�    @�P�    ��  ����      ?}�h@�T@    @�T@    ��  ����      ?��@�X     @�X     ��  ����      ?���@�[�    @�[�    ��  ����      ?�h�@�_�    @�_�    ��  ����      ?�>�@�c@    @�c@    ��  ����      ?��@�g     @�g     ��  ����      ?��@�j�    @�j�    ��  ����      ?��5@�n�    @�n�    ��  ����      ?��s@�r@    @�r@    ��  ����      @��@�v     @�v     ��  ����      @�@�y�    @�y�    ��  ����      @�B@�}�    @�}�    ��  ����      @x~@�@    @�@    ��  ����      @&c�@�     @�     ��  ����      @.O@��    @��    ��  ����      @6:_@ጀ    @ጀ    ��  ����      @>%�@�@    @�@    A,��  �@�c(    @F@�     @�     A4O  �@��    @M��@��    @��    A.Ƞ  �@�L�    @U��@ᛀ    @ᛀ    A'��  �A��    @]�d@�@    @�@    A-"�  �AUt�    @e��@�     @�     Abx  �AR?�    @m�c@��    @��    A�w  �Aj�2    @u��@᪀    @᪀    A�!�  �Ax13    @}�~@�@    @�@    A�S  �A~�[    @���@�     @�     A��.  �A�j�    @��Z@��    @��    A��a  �A�)    @��,@Ṁ    @Ṁ    A�.�  �A���    @��@�@    @�@    A���  �AJu�    @���@��     @��     A�Q  �AJ��    @���@���    @���    A��  �A:7�    @�y�@�Ȁ    @�Ȁ    Aݿ�  �A+��    @�oy@��@    @��@    A�G�  �A�6    @�e^@��     @��     A�[�  �A6    @�[G@���    @���    A��  �@�Vt    @�Q3@�׀    @�׀    A���  �@�y�    @�G#@���    @���    B�3  �@�S    @�)@��    @��    BS  �@�?     @��@��@    @��@    B	  �@��d    @��@��     @��     Bb�  �@���    @�
�@���    @���    B��  �A��    @�@���    @���    B�,  �Ae�    @��	@��@    @��@    BU$  �A4F�    @��@��     @��     B�B  �AXP?    @�� @� �    @� �    B�   �Av��    @��0@��    @��    B�n  �A��%    @��C@�@    @�@    B�Y  �A�b�    @��X@�     @�     B+�  �A�Ț    @�q@��    @��    B�  �A�*    @鱌@��    @��    B~T  �A�H�    @���@�@    @�@    B
l�  �Aۺ    @��@�     @�     BBN  �A��    @���@��    @��    B��  �BE,    @��@�"�    @�"�    A�T�  �Bg9    @��<@�&@    @�&@    A���  �B�    A �4@�*     @�*     A�(S  �B6�    A�K@�-�    @�-�    A�]4  �B#�>    A�c@�1�    @�1�    A���  �B,K�    A�}@�5@    @�5@    A��K  �B/��    A��@�9     @�9     A���  �B9��    A
��@�<�    @�<�    A���  �B?��    A��@�@�    @�@�    A��  �B>�d    A��@�D@    @�D@    A��  �B=5b    A�@�H     @�H     A���  �B<��    A�1@�K�    @�K�    A���  �B<;z    A�S@�O�    @�O�    A�B8  �B7��    A�v@�S@    @�S@    A��  �B4��    A��@�W     @�W     A��s  �B6m�    A{�@�Z�    @�Z�    Bd�  �B3��    Av�@�^�    @�^�    B�m  �B5k(    Ar@�b@    @�b@    B	
t  �B6�]    A m:@�f     @�f     B*�  �B4��    A"hd@�i�    @�i�    B��  �B3�    A$c�@�m�    @�m�    B�R  �B3L�    A&^�@�q@    @�q@    By  �B4�X    A(Y�@�u     @�u     B��  �B7[�    A*U@�x�    @�x�    BCX  �B7�    A,PI@�|�    @�|�    B%��  �B5��    A.K{@�@    @�@    B*�#  �B7Z{    A0F�@�     @�     B2��  �B6w�    A2A�@��    @��    B8+  �B0�U    A4=@⋀    @⋀    B>�3  �B-P�    A68J@�@    @�@    BB�  �B){�    A83�@�     @�     BD�g  �B(W�    A:.�@��    @��    BIt�  �B$��    A<)�@⚀    @⚀    BJU�  �B$��    A>%*@�@    @�@    BL<3  �B"��    A@ d@�     @�     BK�  �B%#    AB�@��    @��    BF ^  �B+��    AD�@⩀    @⩀    B@�1  �B4-    AF@�@    @�@    B@A�  �B3~�    AHW@�     @�     B8�  �B<�\    AJ�@��    @��    B0��  �BB��    AL�@⸀    @⸀    B3@�  �BA�`    AM�@�@    @�@    B28�  �BD�)    AO�Z@��     @��     B0�$  �BF�&    AQ��@���    @���    B��  �B[�    AS��@�ǀ    @�ǀ    B�  �B^��    AU�%@��@    @��@    B�  �BhU=    AW�j@��     @��     B�  �Bn��    AY�@���    @���    B�?  �Ba��    A[��@�ր    @�ր    B-�T  �BM�    A]�@@��@    @��@    B5��  �BE�p    A_ԉ@��     @��     BA�g  �B:��    Aa��@���    @���    BG��  �B4��    Ac�@��    @��    BEv3  �B6<�    Ae�h@��@    @��@    BG�  �B4��    Ag��@��     @��     B6��  �BF�    Ai�@���    @���    B+��  �BSD�    Ak�N@��    @��    B6�  �BG��    Am��@��@    @��@    B:�  �BE��    Ao��@��     @��     B6��  �BI��    Aq�;@���    @���    B/e�  �BR�    As��@��    @��    B)��  �BY     Au��@�@    @�@    B~$  �Be�p    Aw�/@�     @�     BA�  �Bo��    Ay��@��    @��    B��  �Bs�f    A{��@��    @��    Bؠ  �Bu�<    A}�)@�@    @�@    B�%  �Bw7    A�~@�     @�     B$�  �ByB    A��i@��    @��    B�_  �B��    A��@�!�    @�!�    A�  �B�Z�    A���@�%@    @�%@    A�9�  �B��3    A��k@�)     @�)     A��  �B��/    A��@�,�    @�,�    A�e  �B�a    A���@�0�    @�0�    A���  �B�    A��p@�4@    @�4@    A��*  �B��    A��@�8     @�8     A�`�  �B���    A���@�;�    @�;�    A�#�  �B���    A��w@�?�    @�?�    A�9%  �B�-�    A��%@�C@    @�C@    A���  �B��C    A���@�G     @�G     Au;�  �B��Q    A���@�J�    @�J�    A^�  �B�?    A��/@�N�    @�N�    AIl  �B��    A���@�R@    @�R@    A3��  �B�Il    A���