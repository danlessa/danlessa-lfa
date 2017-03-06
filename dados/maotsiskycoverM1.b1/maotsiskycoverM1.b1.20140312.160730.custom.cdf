CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.2-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20140312160730.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         100    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      69 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.012922       ccd_width_factor      	0.013000       center_x      240 pixels     center_y      320 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-03-12 19:49:02, using ingest-tsi-12.2-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20140312.160730.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:03:56 2017 UTC       ANDERS_armtime_timestamp      1488830636     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-03-12 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-03-12 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-03-12 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �S��BH  �rdt�M�M@�X@    @�X@    A���  �B��x    B���@�\     @�\     B�\  �B��    B��@�_�    @�_�    B�+  �B�k�    B�2�@�c�    @�c�    A�  �B��*    B�rc@�g@    @�g@    A�Z`  �B�AJ    B��1@�k     @�k     A�ͼ  �B�̎    B���@�n�    @�n�    A��f  �B�(�    B�1�@�r�    @�r�    B�  �B��?    B�q@�v@    @�v@    A�  �B���    B���@�z     @�z     A��H  �B��H    B��@�}�    @�}�    A��  �B�o�    B��J@쁀    @쁀    A�_  �B��,    B���@�@    @�@    A�6  �B��:    B�L�@�     @�     B[�  �B�m    B��@��    @��    A��3  �B�z�    B��@@쐀    @쐀    A��F  �B�i    B��w@�@    @�@    A�:c  �B�a    B�M�@�     @�     A�`5  �B���    B��@��    @��    A���  �B��    B���@쟀    @쟀    A뒃  �B�HA    B��@�@    @�@    B�d  �B�H�    B�N0@�     @�     B�  �B�._    B�N@��    @��    A��  �B��    B��j@쮀    @쮀    B$  �B��r    B���@�@    @�@    B��  �B��1    B�N�@�     @�     B�  �B��    B��@��    @��    B9�  �B|[    B���@콀    @콀    B�  �B�;E    B���@��@    @��@    A�f  �B�    B�O@��     @��     A�V�  �B�}    B� @���    @���    A�\  �B���    B��8@�̀    @�̀    B��  �B�^�    B��Q@��@    @��@    B �  �B�Q\    B�Oi@��     @��     B}�  �B�1    B��@���    @���    B 
�  �B�N�    B�ϙ@�ۀ    @�ۀ    A�:�  �B�Ɩ    B���@��@    @��@    A��  �B��P    B�O�@��     @��     A���  �B��E    B��@���    @���    A�74  �B�ʜ    B���@��    @��    A튥  �B���    B��@��@    @��@    BP  �B��    B�P&@��     @��     Bo�  �B���    B�=@���    @���    Bq�  �B��    B��T@���    @���    BM�  �B��    B��k@��@    @��@    B  �B��J    B�P�@�     @�     B&G  �B��I    B��@��    @��    B�-  �B�?�    B�б@��    @��    B(  �B��    B���@�@    @�@    A�͹  �B���    B�P�@�     @�     A�1�  �B�$�    B��@��    @��    A�hE  �B� �    B��@��    @��    B��  �B��u    B��#@�@    @�@    B"0  �B���    B�Q:@�     @�     Bq�  �B|b�    B�Q@�"�    @�"�    B	��  �B���    B��h@�&�    @�&�    A��t  �B���    B��~@�*@    @�*@    A�)h  �B��C    B�Q�@�.     @�.     B1�  �B�an    B��@�1�    @�1�    BE�  �B��l    B���@�5�    @�5�    A��  �B�U�    B���@�9@    @�9@    A�y,  �B���    B�Q�@�=     @�=     B�D  �B��h    B�@�@�    @�@�    B	�  �B~�    B��@�D�    @�D�    BP�  �Bz�    B��4@�H@    @�H@    B�  �B�~�    B�RK@�L     @�L     A�u  �B��    B�a@�O�    @�O�    A屖  �B�e    B��x@�S�    @�S�    A���  �B���    B���@�W@    @�W@    A��  �B��t    B�R�@�[     @�[     A�G�  �B��    B��@�^�    @�^�    A��w  �B�ד    B���@�b�    @�b�    B�;  �B���    B���@�f@    @�f@    B��  �B��    B�R�@�j     @�j     B��  �Bv�R    B�@�m�    @�m�    B�  �B�e    B��-@�q�    @�q�    B�C  �B�h�    B��C@�u@    @�u@    B̊  �B�I\    B�SZ@�y     @�y     B\  �B��b    B�p@�|�    @�|�    B��  �B�U
    B�Ӈ@퀀    @퀀    B	Z  �B���    B���@�@    @�@    B �  �B���    B�S�@�     @�     B�`  �B��H    B��@��    @��    B�z  �B��o    B���@폀    @폀    B�u  �B���    B���@�@    @�@    B
},  �B�t    B�T@�     @�     Bg�  �B�:    B�$@��    @��    B�  �B���    B��;@힀    @힀    B ؂  �B�C�    B��R@��@    @��@    A��'  �B��!    B�Th@��     @��     B��  �B�(    B�@���    @���    A�ִ  �B��P    B�ԕ@���    @���    A���  �B��    B���@��@    @��@    A�.�  �B���    B�T�@��     @��     A�ڢ  �B�x    B��@���    @���    A���  �B��\    B���@���    @���    A�~�  �B���    B��@��@    @��@    A��J  �B��T    B�U@��     @��     A�=U  �B�:*    B�3@���    @���    A� �  �B��    B��I@�ˀ    @�ˀ    A֙  �B��v    B��`@��@    @��@    A��  �B�Ǽ    B�Uv@��     @��     A�h�  �B��    B��@���    @���    A���  �B�/�    B�գ@�ڀ    @�ڀ    A�    �B�7    B���@��@    @��@    A׊q  �B���    B�U�@��     @��     A�!  �B�f�    B��@���    @���    A���  �B�:�    B���@��    @��    A��"  �B�    B��@��@    @��@    A�x�  �B�?^    B�V*@��     @��     A��c  �B�)�    B�A@���    @���    A��;  �B���    B��X@���    @���    B �'  �B�t    B��n@��@    @��@    A�F4  �B��    B�V�@�      @�      A���  �B��5    B��@��    @��    A�"  �B��    B�ֲ@��    @��    A�`�  �B�h�    B���@�@    @�@    A�?  �B�އ    B�V�@�     @�     A�aE  �B�/�    B��@��    @��    A��8  �B��v    B��@��    @��    A�`�  �B�~�    B��"@�@    @�@    A�P�  �B�O�    B�W9@�     @�     Aؙ�  �B�J�    B�P@�!�    @�!�    Aӯ  �B�lU    B��f@�%�    @�%�    Aɻ�  �B�    B��}@�)@    @�)@    A�+O  �B��	    B�W�@�-     @�-     A��  �B�`F    B��@�0�    @�0�    A��.  �B�n    B���@�4�    @�4�    A�  �B�'V    B���@�8@    @�8@    A�J�  �B��p    B�W�@�<     @�<     A��  �B��    B�@�?�    @�?�    A��x  �B���    B��@�C�    @�C�    A��  �B�k�    B��1@�G@    @�G@    A�   �B���    B�XH@�K     @�K     A�U�  �B�^!    B�_@�N�    @�N�    A���  �B���    B��u@�R�    @�R�    A�<�  �B��`    B���@�V@    @�V@    A���  �B��n    B�X�@�Z     @�Z     A�mj  �B�&�    B��@�]�    @�]�    A�(�  �B��1    B���@�a�    @�a�    A��  �B���    B���@�e@    @�e@    Aܴ@  �B�F�    B�X�@�i     @�i     A�y�  �B�I�    B�@�l�    @�l�    A��  �B�V�    B��*@�p�    @�p�    A��  �B�3�    B��A@�t@    @�t@    A�/�  �B�"D    B�YX@�x     @�x     A��Z  �B���    B�n@�{�    @�{�    A��i  �B�O�    B�م@��    @��    A�Q�  �B�7y    B���@�@    @�@    A��z  �B��U    B�Y�@�     @�     A�hB  �B��    B��@��    @��    A��  �B���    B���@    @    A�%  �B���    B���@�@    @�@    A��z  �B�ji    B�Z@�     @�     A��  �B�bW    B�$@��    @��    Bg�  �B��V    B��;@    @    B 0�  �B�~�    B��Q@�@    @�@    Bgk  �B�Zh    B�Zh@�     @�     BÏ  �B|�~    B�@��    @��    BI  �By�    B�ږ@    @    B�  �Bs?�    B���@�@    @�@    B�  �B��e    B�Z�@�     @�     B&s  �B�     B��@��    @��    B D=  �B�i<    B���@    @    B	��  �B���    B��@�@    @�@    B�W  �B}�    B�[@��     @��     B�K  �By~\    B�5@���    @���    B��  �BuKg    B��L@�ʀ    @�ʀ    B=+  �B~�    B��c@��@    @��@    B�(  �B�H�    B�[y@��     @��     B�  �B�    B��@���    @���    B�/  �B�"�    B�ۧ@�ـ    @�ـ    B�  �B��     B���@��@    @��@    B��  �B��S    B�[�@��     @��     B	7_  �B��;    B��@���    @���    B<z  �B�V�    B��@��    @��    B��  �B��    B��@��@    @��@    B��  �B{�7    B�\0@��     @��     B�r  �By�    B�G@���    @���    B��  �Bs�b    B��^@���    @���    BN�  �B~�    B��u@��@    @��@    Bd�  �B�n,    B�\�@��     @��     BYO  �B�    B��@��    @��    B��  �B���    B�ܺ@��    @��    A��o  �B��l    B���@�
@    @�
@    A�\�  �B���    B�\�@�     @�     A։V  �B��y    B��@��    @��    A�u  �B�/G    B��@��    @��    A�O  �B�0O    B��,@�@    @�@    A��  �B��    B�]C@�     @�     A���  �B���    B�Z@� �    @� �    B�%  �B��G    B��q@�$�    @�$�    B�  �B���    B���@�(@    @�(@    Bc  �B~�{    B�]�@�,     @�,     B+��  �Bb]    B��@�/�    @�/�    B>#  �Bq��    B���@�3�    @�3�    B��  �Bp,�    B���@�7@    @�7@    B�Y  �Bn3�    B�]�@�;     @�;     B�  �Br�    B�@�>�    @�>�    Bc  �Bo�j    B��)@�B�    @�B�    B�  �Br�    B��A@�F@    @�F@    B��  �B{D�    B�^X@�J     @�J     BH�  �B��w    B�o@�M�    @�M�    B��  �B��    B�ކ@�Q�    @�Q�    BF�  �B}R�    B���@�U@    @�U@    B�   �B{ل    B�^�@�Y     @�Y     B�1  �Bv�T    B��@�\�    @�\�    B:   �B{��    B���@�`�    @�`�    B��  �By�    B���@�d@    @�d@    B/h  �Bs�j    B�_@�h     @�h     BM
  �Br��    B�(@�k�    @�k�    B�  �BwA�    B��?@�o�    @�o�    B��  �B|+�    B��V@�s@    @�s@    A�܂  �B���    B�_m@�w     @�w     A�3D  �B���    B��@�z�    @�z�    B/R  �B�    B�ߜ@�~�    @�~�    B �  �B~�f    B���@�@    @�@    BZ�  �B��    B�_�@�     @�     Ba*  �Bt}    B��@��    @��    B*  �B�y8    B��@    @    B�r  �B|r    B@ @�@    @�@    B	�Y  �B�@F    B~�O@�     @�     B  �B�`    B~@}@��    @��    B�_  �B,�    B}��@    @    B�u  �B}�:    B}@�@�     @�     B�  �B���    B|A8@��    @��    B{0  �B�[T    B{�g@變    @變    B%  �B��1    B{A�@�@    @�@    B��  �B��    Bz��@�     @�     B@[  �B��A    BzA�@��    @��    Bp�  �B~�F    By�"@ﺀ    @ﺀ    B��  �B{X�    ByBP@�@    @�@    B��  �B~��    Bx�@��     @��     B�8  �B�Um    BxB�@���    @���    B��  �B�3    Bw��@�ɀ    @�ɀ    B	�  �B�X+    BwC@��@    @��@    B��  �By��    Bv�;@��     @��     B��  �B��p    BvCj@���    @���    BV<  �B���    BuÙ@�؀    @�؀    B��  �B�V-    BuC�@��@    @��@    B	>  �B���    Bt��@��     @��     A���  �B���    BtD&@���    @���    B !H  �B�!�    Bs�U@��    @��    B	T;  �B�~�    BsD�@��@    @��@    Bb�  �B���    Brĳ@��     @��     B�  �B��7    BrD�@���    @���    B��  �Bz��    Bq�@���    @���    B�  �B|0�    BqEA@��@    @��@    B�p  �BxsM    Bp�p@��     @��     B��  �B�T�    BpE�@� �    @� �    B��  �B~��    Bo��@��    @��    B��  �B���    BoE�@��    @��    B�w  �B�V�    Bn�.@��    @��    B�w  �B~V�    BnF]@�`    @�`    B�o  �B�=�    Bmƌ@�
@    @�
@    B�  �B���    BmF�@�     @�     B��  �B�Q�    Bl��@�     @�     B��  �B�U%    BlG@��    @��    B4I  �B���    Bk�J@��    @��    B$�  �B�-&    BkGz@��    @��    BS�  �B|8�    Bjǩ@��    @��    B �  �Be�    BjG�@�`    @�`    B��  �Bu�    Bi�	@�@    @�@    B�  �BtyQ    BiH8@�     @�     B�  �B}��    Bh�h@�     @�     B�  �Bv�p    BhH�@��    @��    B`2  �Bz4�    Bg��@� �    @� �    B��  �B��r    BgH�@�"�    @�"�    B7+  �BzX^    Bf�'@�$�    @�$�    B$�  �Buu    BfIW@�&`    @�&`    Bc�  �By�    Beɇ@�(@    @�(@    B$�  �B{e�    BeI�@�*     @�*     B�C  �B~�    Bd��@�,     @�,     B-  �Bz��    BdJ@�-�    @�-�    B\�  �Bz5�    Bc�F@�/�    @�/�    B-�  �BtU�    BcJw@�1�    @�1�    B��  �B|��    Bbʧ@�3�    @�3�    BSw  �Bz5�    BbJ�@�5`    @�5`    B`  �By�D    Ba�@�7@    @�7@    B�9  �B{��    BaK7@�9     @�9     B4  �B��3    B`�g@�;     @�;     B[U  �Bp9�    B`K�@�<�    @�<�    B!�  �Bx`�    B_��@�>�    @�>�    B��  �By�{    B_K�@�@�    @�@�    B�r  �Bu��    B^�(@�B�    @�B�    B��  �Bz    B^LX@�D`    @�D`    B��  �B��8    B]̉@�F@    @�F@    B�p  �Bw��    B]L�@�H     @�H     B!�!  �Bl��    B\��@�J     @�J     B {�  �BnG    B\M@�K�    @�K�    B%QY  �Bi;    B[�K@�M�    @�M�    B&E�  �BhJc    B[M{@�O�    @�O�    B%��  �Bh��    BZͬ@�Q�    @�Q�    B �  �Bn�9    BZM�@�S`    @�S`    B �  �Bn��    BY�@�U@    @�U@    B#��  �Bk�    BYN=@�W     @�W     B,��  �Ba�}    BX�n@�Y     @�Y     B1j�  �B],t    BXN�@�Z�    @�Z�    B5��  �BX�Z    BW��@�\�    @�\�    B03Q  �B^XG    BWO @�^�    @�^�    B,~  �BbD>    BV�1@�`�    @�`�    B*�  �Bd��    BVOb@�b`    @�b`    B,pt  �Bb/m    BUϓ@�d@    @�d@    B+��  �Bb�n    BUO�@�f     @�f     B+n  �Bc�V    BT��@�h     @�h     B-��  �B`��    BTP&@�i�    @�i�    B587  �BYK<    BS�W@�k�    @�k�    B2�  �B\��    BSP�@�m�    @�m�    B71  �BWk4    BRй@�o�    @�o�    B0_  �B^N�    BRP�@�q`    @�q`    B1L�  �B].	    BQ�@�s@    @�s@    B1�E  �B\�<    BQQM@�u     @�u     B0��  �B]��    BP�~@�w     @�w     B3z�  �B[l    BPQ�@�x�    @�x�    B9T  �BU}�    BO��@�z�    @�z�    BA�%  �BL�4    BOR@�|�    @�|�    B8%x  �BVe    BN�C@�~�    @�~�    B1��  �B\�r    BNRu@��`    @��`    B.&�  �B`p�    BMҦ@��@    @��@    B0��  �B]�&    BMR�@��     @��     B2?�  �B\R�    BL�	@��     @��     B2��  �B[��    BLS;@���    @���    B.�o  �B_�,    BK�m@���    @���    B/��  �B^�&    BKS�@���    @���    B-7D  �Baj	    BJ��@���    @���    B,��  �Ba��    BJT@��`    @��`    B1+  �B]Y�    BI�4@�@    @�@    B0PH  �B^A;    BITf@�     @�     B.7�  �B`OF    BHԗ@�     @�     B+��  �Bb��    BHT�@��    @��    B-�  �Ba��    BG��@��    @��    B6F;  �BXU�    BGU-@�    @�    B4��  �BY�g    BF�_@�    @�    B4��  �BY�	    BFU�@�`    @�`    B-:�  �Ba9�    BE��@�@    @�@    B-NN  �BaI�    BEU�@�     @�     B)[�  �Be4\    BD�(@�     @�     B.Z�  �B`'�    BDVZ@��    @��    B,ZX  �Bb/m    BC֌@��    @��    B.��  �B`�    BCV�@�    @�    B/��  �B^�    BB��@�    @�    B1�+  �B\�q    BBW#@�`    @�`    B2L�  �B\SU    BA�V@�@    @�@    B.��  �B_��    BAW�@�     @�     B.�  �B_��    B@׻@�     @�     B1�R  �B\َ    B@W�@��    @��    B5�o  �BX�    B?� @��    @��    B1��  �B]    B?XS@�    @�    B8f�  �BV/d    B>؅@�    @�    B6k�  �BX*z    B>X�@�`    @�`    B8o�  �BV(H    B=��@�@    @�@    B8��  �BU��    B=Y@��     @��     B7
�  �BWp    B<�Q@��     @��     B7�E  �BV�.    B<Y�@���    @���    B8}�  �BV�    B;ٷ@���    @���    B=��  �BP��    B;Y�@�Ǡ    @�Ǡ    B;�  �BR��    B:�@�ɀ    @�ɀ    B=�  �BP��    B:ZP@��`    @��`    BD�  �BIՓ    B9ڃ@��@    @��@    BF5G  �BHe�    B9Z�@��     @��     BD��  �BI_X    B8��@��     @��     BH�a  �BE�S    B8[@���    @���    BIV  �BE_P    B7�P@���    @���    BMO�  �BA;�    B7[�@�֠    @�֠    BM@�  �BAHc    B6۷@�؀    @�؀    BQ�]  �B<�X    B6[�@��`    @��`    BJE  �BD(<    B5�@��@    @��@    BI�M  �BD�g    B5\Q@��     @��     BP(  �B>$�    B4܅@��     @��     BU�  �B8��    B4\�@���    @���    BW$  �B7i�    B3��@���    @���    BZ-�  �B4`n    B3] @��    @��    B\m�  �B2�    B2�T@��    @��    B^T.  �B0�    B2]�@��`    @��`    B_%w  �B/0    B1ݼ@��@    @��@    BbR�  �B,�    B1]�@��     @��     B^I�  �B0�    B0�$@��     @��     Bfܿ  �B'��    B0^X@���    @���    Bc2�  �B++�    B/ތ@���    @���    Bd�  �B)�/    B/^�@���    @���    Bf��  �B'֢    B.��@���    @���    Bf�w  �B'�    B._(@��`    @��`    Bfr  �B'Ү    B-�]@��@    @��@    Bj+�  �B$�    B-_�@��     @��     Bk�:  �B"�2    B,��@��     @��     Bh��  �B%�{    B,_�@���    @���    Bg�  �B&�n    B+�.@��    @��    Bl�  �B"@�    B+`c@��    @��    Bj9N  �B$�    B*��@��    @��    Bp�4  �B��    B*`�@�`    @�`    Bo#Y  �B22    B)�@�	@    @�	@    Bn��  �B��    B)a5@�     @�     BmO�  �B!&    B(�j@�     @�     BhMM  �B&�    B(a�@��    @��    Bn^�  �B�i    B'��@��    @��    Bjr�  �B#�p    B'b	@��    @��    Bkpv  �B"�=    B&�>@��    @��    Bma4  �B �G    B&bs@�`    @�`    Bl��  �B!{d    B%�@�@    @�@    Bo�o  �B��    B%b�@�     @�     Bp�o  �B�_    B$�@�     @�     Br^3  �B�    B$cH@��    @��    Bs��  �B��    B#�}@��    @��    BsB  �BB    B#c�@�!�    @�!�    Bm�@  �B Q�    B"��@�#�    @�#�    Bt��  �B�    B"d@�%`    @�%`    Bt�  �B�*    B!�S@�'@    @�'@    Bv�  �B�    B!d�@�)     @�)     B}y�  �B��    B �@�+     @�+     By�  �B�.    B d�@�,�    @�,�    Bv�\  �B=�    B�*@�.�    @�.�    Bq��  �Bvj    Be`@�0�    @�0�    Bw�  �B `    B�@�2�    @�2�    Bw�Q  �Bc�    Be�@�4`    @�4`    Bxu�  �B��    B�@�6@    @�6@    B}qz  �B��    Bf8@�8     @�8     B{��  �B��    B�n@�:     @�:     Bz8�  �B,\    Bf�@�;�    @�;�    ByCY  �B�    B��@�=�    @�=�    By"�  �B��    Bg@�?�    @�?�    Bid  �B$�v    B�G@�A�    @�A�    Bx�'  �BL�    Bg~@�C`    @�C`    B|��  �BI%    B�@�E@    @�E@    Bt��  �B��    Bg�@�G     @�G     Bm"�  �B!    B�!@�I     @�I     BkS�  �B"�    BhX@�J�    @�J�    Bns:  �B�    B�@�L�    @�L�    Bu��  �B�f    Bh�@�N�    @�N�    Bzy�  �B��    B��@�P�    @�P�    Be��  �B(�    Bi3@�R`    @�R`    Bi�?  �B$>>    B�j@�T@    @�T@    BU��  �B8%�    Bi�@�V     @�V     BZ��  �B3e�    B��@�X     @�X     B`dN  �B-�d    Bj@�Y�    @�Y�    BYu  �B4��    B�G@�[�    @�[�    BHx|  �BEØ    Bj~@�]�    @�]�    BI��  �BDd    B�@�_�    @�_�    BL�0  �BA�4    Bj�@�a`    @�a`    BS��  �B:��    B�$@�c@    @�c@    BT�8  �B9��    Bk\@�e     @�e     BF�  �BG��    B�@�g     @�g     BH�  �BFOp    Bk�@�h�    @�h�    BGZ'  �BF�    B�@�j�    @�j�    BI�
  �BD�    Bl:@�l�    @�l�    BK��  �BBf�    B�r@�n�    @�n�    BUN�  �B9�    Bl�@�p`    @�p`    B_�  �B/)�    B��@�r@    @�r@    B]�  �B1u    Bm@�t     @�t     BV�N  �B7|~    B�R@�v     @�v     BV�(  �B6�;    Bm�@�w�    @�w�    BR]�  �B;m�    B��@�y�    @�y�    BG&  �BF�7    Bm�@�{�    @�{�    BE,�  �BH��    B
�3@�}�    @�}�    BD$�  �BJ0    B
nl@�`    @�`    BF�T  �BGS�    B	�@�@    @�@    BA��  �BL�>    B	n�@�     @�     BC��  �BJ��    B�@�     @�     BJ��  �BC�t    BoN@��    @��    BKI  �BB��    B�@��    @��    B>J�  �BO�V    Bo�@�    @�    B-��  �B`^�    B��@�    @�    B4��  �BYS�    Bp2@�`    @�`    B<�7  �BQ�:    B�k@�@    @�@    B9�  �BT[
    Bp�@�     @�     B0��  �B]�    B��@�     @�     BΔ  �Bp��    Bq@��    @��    B��  �Bt��    B�P@��    @��    BU�  �BvU    Bq�@�    @�    B��  �Br��    B��@�    @�    B�{  �Bz��    Bq�@�`    @�`    B=�  �Bw!n    B�6@�@    @�@    B�?  �Bt/    Bro@�     @�     B4(  �B}��    B �@�     @�     BӃ  �Buub    B r�@��    @��    B!�8  �BlXv    A��:@��    @��    B�o  �Bu��    A��@�    @�    B"O5  �Bk�{    A��"@�    @�    B)u�  �Bd�j    A��@�`    @�`    B-   �Ba�    A��@�@    @�@    B"6�  �Bk��    A��@�     @�     B�"  �BpA�    A���@�     @�     BB�  �By
J    A��i@��    @��    B��  �Bww�    A���@��    @��    B'�  �Bx�    A��S@�    @�    B�  �Bu�    A���@�    @�    B�  �Br�    A��>@�`    @�`    B��  �BtN�    A��@�@    @�@    B�  �BuM�    A��)@�     @�     B4  �Bp7    A��@��     @��     B#  �Bj�K    A��@���    @���    B!΂  �BlK�    A��@���    @���    B$�R  �Bil    A��@�Ơ    @�Ơ    B'7  �Bf�    A��x@�Ȁ    @�Ȁ    B+�@  �Ba��    A���@��`    @��`    B(P1  �Be�R    A��f@��@    @��@    B/�  �B^�R    A���@��     @��     B-�N  �B`'    A��T@��     @��     B3N  �BZ��    A���@���    @���    B939  �BT�c    A��C@���    @���    B9y�  �BT}B    A��@�ՠ    @�ՠ    B>��  �BN�    A��3@�׀    @�׀    B<<T  �BQ��    A��@��`    @��`    B?�  �BN�1    A��#@��@    @��@    BAg�  �BLo?    A��@��     @��     BE*�  �BH�P    A��@��     @��     BFo�  �BGT;    A��@���    @���    BH>  �BE˻    A��@���    @���    BK��  �BBJ�    A��~@��    @��    BN�  �B>�c    A���@��    @��    BNٟ  �B>�2    A��p@��`    @��`    BS>  �B:��    A���@��@    @��@    BV$}  �B7l%    A��c@��     @��     BX#�  �B5�    A���@��     @��     BW�  �B5�V    A��W@���    @���    B^G  �B/1�    A���@���    @���    B`��  �B,��    A��L@��    @��    Bce  �B)��    A���@���    @���    Bgnt  �B%�    A��A@��`    @��`    Bgd.  �B&6�    A���@��@    @��@    Bk�O  �B!��    A��7@��     @��     Bo�w  �B��    A���@��     @��     Bp�  �Bm    A��-@���    @���    Bt�M  �B�N    A���@� �    @� �    Bu�z  �B��    A��%@��    @��    B|�  �Bp�    A���@��    @��    B�+  �B�    A��@�`    @�`    B��I  �B�    A���@�@    @�@    B��  �B�    A��@�
     @�
     B���  �B7;    A���@�     @�     B�'  �B��    A� @��    @��    B�j�  �A���    A� �@��    @��    B���  �A�ت    A�	@��    @��    B��>  �A�$�    A��@��    @��    B��@  �A�Q�    A�@�`    @�`    B�x  �A�K/    A��@�@    @�@    B��l  �A߉�    A��@�     @�     B�C�  �A��     A�}@�     @�     B�2�  �A�/"    A��@��    @��    B��d  �A�$�    A�z@��    @��    B�[�  �A��7    A��@� �    @� �    B�	�  �A�I"    A�w@�"�    @�"�    B��&  �Aʶ�    A��@�$`    @�$`    B���  �A��    A�v@�&@    @�&@    B�ɶ  �A�L�    A��@�(     @�(     B���  �A�O�    A�t@�*     @�*     B��/  �A�|�    A��@�+�    @�+�    B���  �A�TA    A�t@�-�    @�-�    B�  �A�l�    A��@�/�    @�/�    B�ٷ  �A��3    A�	u@�1�    @�1�    B��S  �A�Z�    A�	�@�3`    @�3`    B��"  �B ��    A�
v@�5@    @�5@    B��  �B&D    A�
�@�7     @�7     B�MU  �B.�    A�x@�9     @�9     B��  �B�    A��@�:�    @�:�    B|�-  �B��    A�{@�<�    @�<�    Bw�|  �B��    A��@�>�    @�>�    Bu��  �B�    A�@�@�    @�@�    Bo#E  �B��    A�@�B`    @�B`    Bk�N  �B ��    A��@�D@    @�D@    Bk^v  �B!Ȝ    A�@�F     @�F     BhER  �B$e�    A��@�H     @�H     BhcQ  �B$n<    A�@�I�    @�I�    Bg?O  �B%�"    A��@�K�    @�K�    Bb��  �B*��    A�@�M�    @�M�    B`  �B-5    A��@�O�    @�O�    B^�  �B/!�    A�@�Q`    @�Q`    B]�  �B0hE    A��@�S@    @�S@    B\N`  �B0�Q    A�"@�U     @�U     B[@�  �B2%    A��@�W     @�W     BY@�  �B3��    A�+@�X�    @�X�    BU�}  �B7��    A��@�Z�    @�Z�    BU��  �B7�m    A�4@�\�    @�\�    BS��  �B9m�    A��@�^�    @�^�    BP��  �B<tn    A�?@�``    @�``    BM�{  �B?��    A��@�b@    @�b@    BJ��  �BB�q    A�K@�d     @�d     BM�l  �B?�5    A��@�f     @�f     BI�+  �BC�/    A�W@�g�    @�g�    BF7�  �BG's    A��@�i�    @�i�    BG1�  �BFf,    A�e@�k�    @�k�    BB��  �BJ�Q    A��@�m�    @�m�    BAt  �BK�p    A�s@�o`    @�o`    B@�]  �BLJ�    A��@�q@    @�q@    B>�5  �BN��    A��@�s     @�s     B7}�  �BU��    A�
@�u     @�u     B7^  �BU��    A��@�v�    @�v�    B4H  �BX��    A�@�x�    @�x�    B3�k  �BY��    A��@�z�    @�z�    B3�%  �BX��    A�,@�|�    @�|�    B36�  �BY��    A��@�~`    @�~`    B.E�  �B^��    A�>@�@    @�@    B)8�  �Bc��    A��@�     @�     B*V�  �Bb`|    A� Q@�     @�     B*�N  �Bb�    A� �@��    @��    B$��  �Bh    A�!e@��    @��    B)_�  �Bb�
    A�!�@�    @�    B*�Z  �Baue    A�"z@�    @�    B)�E  �Ba��    A�#@�`    @�`    B%˾  �Bf�    A�#�@�@    @�@    B$F�  �Bh?Q    A�$@�     @�     B&;�  �BeI�    A�$�@�     @�     B'��  �Bc��    A�%3@��    @��    B&�Z  �Bd��    A�%�@��    @��    B'G2  �Bd]Q    A~L�@�    @�    B&�  �Bd�#    A|M�@�    @�    B(��  �Bb@\    AzN�@�`    @�`    B(�  �Bb��    AxO�@�@    @�@    B)A�  �BbDO    AvP�@�     @�     B+��  �B_R�    AtR@�     @�     B,��  �B^i�    ArS6@��    @��    B/C  �B[�    ApTR@��    @��    B0��  �BZ.g    AnUn@�    @�    B2�V  �BX8    AlV�@�    @�    B2��  �BXM�    AjW�@�`    @�`    B5=  �BUV    AhX�@�@    @�@    B7��  �BRz�    AfY�@�     @�     B6b�  �BS|�    Ad[@�     @�     B9  �BQ,�    Ab\%@��    @��    B6_  �BS�s    A`]E@��    @��    B8}  �BQ��    A^^f@�    @�    B9;�  �BP�    A\_�@�    @�    B:�@  �BN��    AZ`�@�`    @�`    B:8a  �BO��    AXa�@�@    @�@    B<�'  �BL�W    AVb�@�     @�     B<f�  �BL��    ATd@��     @��     B>$  �BK�    ARe4@���    @���    B?h  �BJ(�    APfY@���    @���    BB-  �BF��    ANg~@�Š    @�Š    BA2  �BG�R    ALh�@�ǀ    @�ǀ    BDP�  �BD�Y    AJi�@��`    @��`    BE��  �BCV�    AHj�@��@    @��@    BF&t  �BBA�    AFl@��     @��     BF�r  �BA�i    ADm>@��     @��     BI=�  �B?:�    ABnf@���    @���    BH�i  �B>�    A@o�@���    @���    BK_  �B<�e    A>p�@�Ԡ    @�Ԡ    BN�  �B9q    A<q�@�ր    @�ր    BNd*  �B9�(    A:s@��`    @��`    BOŖ  �B7��    A8t7@��@    @��@    BQKN  �B6��    A6ub@��     @��     BQ��  �B5��    A4v�@��     @��     BS��  �B2�    A2w�@���    @���    BS��  �B2�    A0x�@���    @���    BRf�  �B3     A.z@��    @��    BS#  �B1v    A,{D@��    @��    BRΦ  �B1�    A*|s@��`    @��`    BPB<  �B2-    A(}�@��@    @��@    BO�� �B3S    A&~�@��     @��     BN:� �B2'8    A$�@��     @��     BM~� �B4l    A"�4@���    @���    BI� �B53�    A �f@���    @���    BH=* �B5�    A��@��    @��    BF�o �B6c�    A��@��    @��    BO`4 �B-�)    A��@��`    @��`    BE� �B5�    A�3@��@    @��@    BFQ� �B5g    A�g@��     @��     BF7� �B3�Q    A��@��     @��     BCzS �B6�9    A��@���    @���    BF�� �B19    A�	@���    @���    B@l� �B5b�    A�@@��    @��    BA}� �B5ok    A�x@��    @��    BC�� �B1�    A
��@�`    @�`    BE�� �B.:=    A��@�@    @�@    BKo� �B(��    A�#@�	     @�	     BMi� �B$�~    A�]@�     @�     BMrV �B%�Z    A��@��    @��    BK*� �B&e�    A ��@��    @��    BI0 �B%�    @�.@��    @��    BB�f �B+�    @�0�@��    @��    BEϫ �B& �    @�3@�`    @�`    BD	 �B%1�    @�5�@�@    @�@    BA�� �B#�'    @�8@�     @�     BGX �B�    @�:�@�     @�     BJY� �B�3    @�=	@��    @��    BS�$  �BT.    @�?�@��    @��    Bb��  �Bi�    @�B@��    @��    B^9�  �B�    @�D�@�!�    @�!�    B_�  �B��    @�G@�#`    @�#`    B]uV  �B��    @�I�@�%@    @�%@    B\��  �B�    @�L#@�'     @�'     B[�w  �B    @�N�@�)     @�)     BZi�  �B 7�    @�Q6@�*�    @�*�    BY�  �A�    @�S�@�,�    @�,�    BV��  �A�,    @�VO@�.�    @�.�    BR��  �A��    @�X�@�0�    @�0�    BOU�  �A�{U    @�[l@�2`    @�2`    BM�8  �B ��    @�]�@�4@    @�4@    BH;�  �A���    @�`�@�6     @�6     BE��  �A���    @�c#@�8     @�8     BC`�  �A�<�    @�e�@�9�    @�9�    B@C�  �A��}    @�hO@�;�    @�;�    B=�Q  �A��L    @�j�@�=�    @�=�    B="�  �A�ۧ    @�m�@�?�    @�?�    B:��  �A�q    @�p@�A`    @�A`    B;C<  �A�?,    @�r�@�C@    @�C@    B8AM  �A܁~    @�uT@�E     @�E     B7-�  �A�&�    @�w�@�G     @�G     B6G�  �A��z    @�z�@�H�    @�H�    B3�T  �A��M    @�}5@�J�    @�J�    B38L  �A��K    @z��@�L�    @�L�    B/��  �A�?    @s�@�N�    @�N�    B+�_  �A��E    @k
E@�P`    @�P`    B*��  �A���    @c�@�R@    @�R@    B%ɥ  �A�_4    @[�@�T     @�T     B!   �A�o    @S;@�V     @�V     B�v  �A��    @K�@�W�    @�W�    BvD  �A�D�    @C$�@�Y�    @�Y�    ��  ����      @;*K@�[�    @�[�    ��  ����      @3/�@�]�    @�]�    ��  ����      @+5@�_`    @�_`    ��  ����      @#:v@�a@    @�a@    ��  ����      @?�@�c     @�c     ��  ����      @EL@�e     @�e     ��  ����      @J�@�f�    @�f�    ��  ����      @P/@�h�    @�h�    ��  ����      ?��J@�j�    @�j�    ��  ����      ?�<@�l�    @�l�    ��  ����      ?��4@�n`    @�n`    ��  ����      ?��3@�p@    @�p@    ��  ����      ?��7@�r     @�r     ��  ����      ?��B@�t     @�t     ��  ����      ?��T@�u�    @�u�    ��  ����      ?��k@�w�    @�w�    ��  ����      ?n@�y�    @�y�    ��  ����      ?N[@�{�    @�{�    ��  ����      ?.3�@�}`    @�}`    ��  ����      ?J@�@    @�@    ��  ����      >��@�     @�     ��  ����      >���@�     @�     ��  ����      >:6$@��    @��    ��  ����      =jA�