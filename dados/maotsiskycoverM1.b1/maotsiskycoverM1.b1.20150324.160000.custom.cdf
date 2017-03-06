CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.4-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20150324160000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      64 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      224 pixels     center_y      314 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Zcreated by user dsmgr on machine ruby at 2015-03-24 18:49:02, using ingest-tsi-12.4-0.el6      ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20150324.160000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:28 2017 UTC       ANDERS_armtime_timestamp      1488830668     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2015-03-24 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2015-03-24 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2015-03-24 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �U� BH  �rdt�M�M@�      @�      B�.  �B}�*    B���@�#�    @�#�    B(g  �B��F    B��!@�'�    @�'�    B}U  �B��w    B��@�+@    @�+@    B��  �Bz�q    B���@�/     @�/     B
�S  �B~��    B��@�2�    @�2�    A��F  �B�R    B�"�@�6�    @�6�    B ��  �B�7�    B�6	@�:@    @�:@    B}�  �B}A�    B�G�@�>     @�>     B�d  �B� �    B�X,@�A�    @�A�    A� �  �B���    B�f�@�E�    @�E�    A�_z  �B���    B�t@�I@    @�I@    A�f�  �B��m    B��@�M     @�M     A�R�  �B�=    B���@�P�    @�P�    A��  �B���    B���@�T�    @�T�    A��L  �B��    B��`@�X@    @�X@    B˹  �B�	    B��F@�\     @�\     A�1  �B���    B��y@�_�    @�_�    B�  �Bl��    B���@�c�    @�c�    B�A  �B���    B���@�g@    @�g@    B��  �Bq��    B���@�k     @�k     B�t  �B�4�    B��7@�n�    @�n�    A�9�  �B��T    B���@�r�    @�r�    A�'j  �B��    B���@�v@    @�v@    B
�@  �By�    B��;@�z     @�z     B(>  �Bl��    B�|�@�}�    @�}�    BGz  �Bh�    B�p�@쁀    @쁀    A���  �B��    B�ch@�@    @�@    B�>  �B~׮    B�TE@�     @�     B�E  �B_�    B�C�@��    @��    B�  �Bi.x    B�1`@쐀    @쐀    B$  �Bxia    B��@�@    @�@    A��Z  �B�7@    B��@�     @�     A�-$  �B��    B���@��    @��    B��  �Bd�    B��@쟀    @쟀    B&�  �B`e�    B���@�@    @�@    B��  �B\�g    B��#@�     @�     B9�  �BU    B��H@��    @��    BX  �BSgx    B�m2@쮀    @쮀    B
�  �BU\�    B�N�@�@    @�@    A�c&  �Bn1    B�/�@�     @�     Bn�  �BG]�    B��@��    @��    B�  �Bd��    B��_@콀    @콀    A��  �Bu�'    B�ʽ@��@    @��@    B�  �B]�    B��@��     @��     A�6Y  �Bz�    B���@���    @���    A�7z  �B���    B�] @�̀    @�̀    A���  �B��Y    B�6�@��@    @��@    A�l�  �Bz?f    B�S@��     @��     B8�  �BV�    B��=@���    @���    A�eT  �Bw�    B��[@�ۀ    @�ۀ    A롛  �BqP�    B���@��@    @��@    A��X  �B|�(    B�jV@��     @��     Aȥ�  �B��    B�?A@���    @���    A�0]  �B��F    B�@��    @��    A���  �B��'    B��@��@    @��@    Aß�  �B�6�    B��@��     @��     A�vF  �B��~    B��l@���    @���    A��\  �B��    B�^7@���    @���    A�M�  �B��(    B�/s@��@    @��@    A�f�  �B�]    B� '@�     @�     A�+�  �B�Ld    B��Y@��    @��    A��  �B���    B��@��    @��    A�&  �B��    B�oG@�@    @�@    A���  �B��Z    B�>@�     @�     A��;  �By�0    B�c@��    @��    A�P>  �Bg}     B��N@��    @��    B�\  �BNz�    B���@�@    @�@    A��  �Bi2    B�t�@�     @�     B�1  �BX�l    B�A�@�"�    @�"�    A�V�  �B�4    B�@�&�    @�&�    A���  �B|�    B��#@�*@    @�*@    A��g  �B���    B���@�.     @�.     A�[�  �B�6R    B�q>@�1�    @�1�    A���  �B���    B�<S@�5�    @�5�    A�gc  �B��5    B�@�9@    @�9@    A�  �B�n    B�ћ@�=     @�=     A��  �B��+    B���@�@�    @�@�    A�c�  �Bb��    B�e�@�D�    @�D�    A��  �B_K    B�/z@�H@    @�H@    A��l  �Bm�    B���@�L     @�L     A���  �B���    B��"@�O�    @�O�    A���  �Bp�    B��@�S�    @�S�    A�M  �BwsA    B�S�@�W@    @�W@    A��w  �B���    B�h@�[     @�[     A�˩  �B�N�    B��@�^�    @�^�    A��  �B��2    B���@�b�    @�b�    A�z�  �B��    B�t�@�f@    @�f@    A]�"  �B�p    B�<�@�j     @�j     AZ�@  �B��9    B�@�m�    @�m�    A_7  �B�2�    B�ˀ@�q�    @�q�    Ah�  �B��    B���@�u@    @�u@    Aj��  �B��    B�Y�@�y     @�y     A��  �B��    B� �@�|�    @�|�    At�a  �B���    B��i@퀀    @퀀    A��m  �B�Z�    B�� @�@    @�@    A�U  �B�i|    B�ts@�     @�     A�L  �B�=�    B�:�@��    @��    A�O  �B�Zc    B� �@폀    @폀    A�F  �Br�    B���@�@    @�@    A���  �B�*�    B���@�     @�     AuC2  �B��    B�R�@��    @��    Ay��  �B��    B�S@힀    @힀    A}*�  �B���    B���@��@    @��@    A�r�  �B��    B��O@��     @��     A�֭  �B~��    B�h�@���    @���    A�7.  �Bp]H    B�-�@���    @���    A��/  �B�.c    B���@��@    @��@    A�uO  �Bg��    B���@��     @��     A��  �Bi��    B�|�@���    @���    A�e�  �Bdq    B�A�@���    @���    A��  �Bd?�    B�s@��@    @��@    A��  �B_�(    B��@��     @��     A��  �B_�    B���@���    @���    A͈   �Be�M    B�T@�ˀ    @�ˀ    A��u  �Bs     B�|@��@    @��@    A�5  �BZd�    B���@��     @��     A�!a  �BPV�    B��@���    @���    Aۊ   �B5E    B�e(@�ڀ    @�ڀ    AߓX  �B*yQ    B�);@��@    @��@    Aܶ�  �B#��    B��:@��     @��     A���  �B�    B��(@���    @���    A�dx  �BLU    B�u@��    @��    A��  �B:�    B�8�@��@    @��@    A�س  �B�}    B���@��     @��     A�y  �B3[    B��3@���    @���    AѤ�  �B�K    B���@���    @���    A�Y   �BU�    B�GW@��@    @��@    A�iy  �B~�    B�
�@�      @�      A���  �B%y�    B��?@��    @��    A�z1  �B2��    B���@��    @��    A�6�  �BBY�    B�T�@�@    @�@    A�7!  �BM�*    B�/@�     @�     A���  �BU�m    B��c@��    @��    A��  �B^��    B���@��    @��    A�׏  �BjY    B�a�@�@    @�@    A��  �B���    B�$�@�     @�     A���  �B���    B��@�!�    @�!�    A���  �B��]    B���@�%�    @�%�    Ax�s  �B��T    B�m�@�)@    @�)@    Ae�+  �B���    B�0u@�-     @�-     Ag:  �B�\�    B��I@�0�    @�0�    A\q  �B��    B��@�4�    @�4�    A�U  �B�T�    B�x�@�8@    @�8@    AP��  �B�YJ    B�;�@�<     @�<     A<
�  �B��E    B��-@�?�    @�?�    A�V  �B��    B���@�C�    @�C�    A&x  �B�w�    B��b@�G@    @�G@    A�U  �B�_�    B�E�@�K     @�K     AF~  �B�CC    B�r@�N�    @�N�    A��  �B���    B���@�R�    @�R�    A�D  �B��     B��\@�V@    @�V@    A �  �B�9�    B�O�@�Z     @�Z     A^V  �B���    B�#@�]�    @�]�    A��  �B��J    B��z@�a�    @�a�    Aк  �B�\l    B���@�e@    @�e@    A��  �B��    B�Y@�i     @�i     A'�  �B�T    B�N@�l�    @�l�    A�  �B�    B�݅@�p�    @�p�    A�r  �B�I�    B���@�t@    @�t@    A%	�  �B�^     B�a�@�x     @�x     A`  �B��    B�#�@�{�    @�{�    A*�  �B�C"    B��@��    @��    A/�  �B���    B��(@�@    @�@    A8U�  �B��    B�j3@�     @�     A7�[  �B��z    B�,8@��    @��    A>$t  �B��v    B��5@    @    A7��  �B��
    B��-@�@    @�@    AM�  �B��U    B�r@�     @�     A���  �B��    B�4	@��    @��    Aw6  �B�@�    B���@    @    Aht�  �B�6�    B���@�@    @�@    A]�  �B���    B�y�@�     @�     A{Fy  �B��%    B�;w@��    @��    A[u;  �B�P�    B��E@    @    AF}�  �B�p�    B��@�@    @�@    AIQ�  �B��    B���@�     @�     AM��  �B��b    B�B�@��    @��    AX�  �B���    B�B@    @    Ad�O  �B�    B���@�@    @�@    A�u  �B���    B���@��     @��     A�xR  �B���    B�II@���    @���    A��C  �B���    B�
�@�ʀ    @�ʀ    A�Na  �B�jX    B�̊@��@    @��@    A�@"  �B�H�    B��#@��     @��     A��7  �B�V�    B�O�@���    @���    Aʾ  �B��    B�H@�ـ    @�ـ    A���  �Bq��    B���@��@    @��@    B�1  �BeO    B��[@��     @��     B��  �B_0s    B�U�@���    @���    B�  �BW �    B�[@��    @��    B��  �BS{�    B���@��@    @��@    Br�  �BK:�    B��K@��     @��     B3�  �BD�    B�[�@���    @���    B    �BD�    B�+@���    @���    BI�  �BK��    B�ޕ@��@    @��@    B��  �BE	�    B���@��     @��     B�,  �BG=    B�a]@��    @��    A�S6  �BS�=    B�"�@��    @��    A��$  �BM/&    B��@�
@    @�
@    A�,�  �BVL*    B��l@�     @�     A�C�  �BHK�    B�f�@��    @��    A�bC  �BND�    B�(@��    @��    A�b  �BR      B��[@�@    @�@    A��+  �BG��    B���@�     @�     B�  �B:�|    B�k�@� �    @� �    BX�  �B5�M    B�-,@�$�    @�$�    Bϙ  �B=�    B��j@�(@    @�(@    A���  �B?�    B���@�,     @�,     A���  �BB�s    B�p�@�/�    @�/�    A�)  �BM
U    B�2@�3�    @�3�    A�zC  �B\ �    B��F@�7@    @�7@    A�̺  �Ba��    B��u@�;     @�;     A�\  �Bl�    B�u�@�>�    @�>�    A��H  �Bz�>    B�6�@�B�    @�B�    A��  �B{��    B���@�F@    @�F@    A��  �B���    B��@�J     @�J     A�Z  �B���    B�z4@�M�    @�M�    A��S  �B�4�    B�;R@�Q�    @�Q�    A�_�  �B��z    B��m@�U@    @�U@    A�p�  �B��    B���@�Y     @�Y     A��  �B�i�    B�~�@�\�    @�\�    A�.&  �B��*    B�?�@�`�    @�`�    A��[  �B��B    B� �@�d@    @�d@    A�R�  �B���    B��@�h     @�h     A�v�  �B�.�    B�@�k�    @�k�    A��5  �B� �    B~��@�o�    @�o�    A��  �B�/�    B~	�@�s@    @�s@    A��  �B���    B}��@�w     @�w     A�k  �B�H�    B}�@�z�    @�z�    A��}  �B�K�    B|��@�~�    @�~�    A���  �B�ږ    B|�@�@    @�@    A�M�  �B�d3    B{��@�     @�     A�q�  �B��d    B{�@��    @��    A�2�  �B��    Bz��@    @    Aɇm  �B���    Bz�@�@    @�@    A��  �Bj�P    By�o@�     @�     A�o  �Bi��    ByO@��    @��    A��  �Bq��    Bx�*@    @    A���  �B^lt    Bx!@�@    @�@    BG  �BP��    Bw��@�     @�     B)  �BO&    Bw$�@��    @��    BX�  �BQݜ    Bv�n@變    @變    Bv�  �BH��    Bv(4@�@    @�@    B �  �B;p-    Bu��@�     @�     B&O.  �B;ɖ    Bu+�@��    @��    B(��  �B9�x    Bt�r@ﺀ    @ﺀ    B S�  �BJ�$    Bt/*@�@    @�@    B$<j  �BFE    Bs��@��     @��     B!fs  �BE�    Bs2�@���    @���    B#H  �BA@I    Br�;@�ɀ    @�ɀ    B�L  �BmkT    Br5�@��@    @��@    A�;M  �B|bs    Bq��@��     @��     A��  �B~��    Bq9,@���    @���    A�Bc  �B��     Bp��@�؀    @�؀    A�pd  �B��=    Bp<f@��@    @��@    A�#�  �B�ܟ    Bo��@��     @��     A��*  �B�Sb    Bo?�@���    @���    A�}�  �B�=    Bn�$@��    @��    A�F�  �B���    BnB�@��@    @��@    A�t�  �B~��    Bm�=@��     @��     B Yc  �Bu��    BmE�@���    @���    A���  �Bx    Bl�J@���    @���    A�  �Bz��    BlH�@��@    @��@    B��  �Bn�    Bk�J@��     @��     A��X  �Bn1f    BkK�@� �    @� �    B�j  �Bc�c    Bj�?@��    @��    Bby  �Bb$    BjN�@��    @��    B'~  �B`B    Bi�(@��    @��    A��  �BebL    BiQ�@�`    @�`    A�n  �BgB�    Bh�@�
@    @�
@    A��  �Bi��    BhTo@�     @�     A�C�  �BkG    Bg��@�     @�     A��  �Bd��    BgW<@��    @��    A�i=  �B^M�    Bf؞@��    @��    B�  �B?��    BfY�@��    @��    B��  �BE��    Be�[@��    @��    By  �B<�    Be\�@�`    @�`    B	��  �B?6�    Bd�@�@    @�@    B��  �BH�    Bd_b@�     @�     A��  �BG�O    Bc�@�     @�     A�#�  �BCL�    Bcb@��    @��    A�  �B>�    Bb�S@� �    @� �    A�s  �B;�P    Bbd�@�"�    @�"�    A�4�  �B4�h    Ba��@�$�    @�$�    B �-  �B*v�    Bag.@�&`    @�&`    B�	  �B�    B`�r@�(@    @�(@    BG�  �B!7V    B`i�@�*     @�*     A�M�  �B++    B_��@�,     @�,     B   �BYv    B_l1@�-�    @�-�    B��  �B&G7    B^�l@�/�    @�/�    B
��  �B'�    B^n�@�1�    @�1�    B��  �B)і    B]��@�3�    @�3�    A��  �B6�Z    B]q@�5`    @�5`    A��t  �B4��    B\�B@�7@    @�7@    Aϵ�  �BQ>~    B\sr@�9     @�9     Aĸ:  �BO�f    B[��@�;     @�;     A���  �BL'    B[u�@�<�    @�<�    A��  �BL��    BZ��@�>�    @�>�    A�4]  �BI�1    BZx@�@�    @�@�    A�n�  �BCD    BY�D@�B�    @�B�    A�7\  �BA�    BYzg@�D`    @�D`    A�z  �B84    BX��@�F@    @�F@    A��O  �B4.�    BX|�@�H     @�H     A�;�  �B4�    BW��@�J     @�J     A��  �B�g    BW~�@�K�    @�K�    A�Ӵ  �Bw�    BV��@�M�    @�M�    A�E�  �B�c    BV�@�O�    @�O�    A�N_  �B%    BV,@�Q�    @�Q�    A���  �B�f    BU�@@�S`    @�S`    A٢`  �A�=    BUS@�U@    @�U@    A��  �A�K�    BT�d@�W     @�W     A���  �A���    BTs@�Y     @�Y     A�l  �A�0�    BS��@�Z�    @�Z�    A��  �A��p    BS�@�\�    @�\�    A���  �A窅    BR��@�^�    @�^�    A���  �A�I�    BR
�@�`�    @�`�    A�=  �Aۧ�    BQ��@�b`    @�b`    A�,�  �A��f    BQ�@�d@    @�d@    A�')  �A�fG    BP��@�f     @�f     A�+N  �A�]�    BP�@�h     @�h     A��5  �A���    BO��@�i�    @�i�    A�wS  �A�E�    BO�@�k�    @�k�    A��t  �A���    BN��@�m�    @�m�    A��h  �A��    BN�@�o�    @�o�    A�X�  �A��A    BM��@�q`    @�q`    A�b�  �A�M�    BM�@�s@    @�s@    A�C�  �A�>�    BL��@�u     @�u     A��  �A��_    BL@�w     @�w     A��4  �A�P�    BK�r@�x�    @�x�    Az�C  �A��    BKd@�z�    @�z�    A�{�  �A��    BJ�T@�|�    @�|�    A���  �A�IP    BJC@�~�    @�~�    A�%   �A��h    BI�1@��`    @��`    A��,  �A���    BI@��@    @��@    A�A�  �A�|�    BH�@��     @��     A�n  �A�	6    BH�@��     @��     A��  �A�|�    BG��@���    @���    A��  �A��    BG�@���    @���    A}��  �A���    BF��@���    @���    Av��  �A}�#    BF!�@���    @���    A�"�  �Al�    BE�i@��`    @��`    A�#�  �Aa�*    BE#J@�@    @�@    A�K�  �AZ�t    BD�)@�     @�     A���  �AYix    BD%@�     @�     A��m  �AQ[    BC��@��    @��    A� �  �AM�k    BC&�@��    @��    A�:  �AJ�r    BB��@�    @�    A�Ĩ  �AH��    BB(s@�    @�    A��v  �AA�d    BA�K@�`    @�`    A��Z  �AE��    BA*"@�@    @�@    A��  �AE�y    B@��@�     @�     A��  �AE�E    B@+�@�     @�     A��  �AJO    B?��@��    @��    A���  �AN�4    B?-p@��    @��    A��  �AN�2    B>�@@�    @�    A���  �AI��    B>/@�    @�    A��8  �AL�K    B=��@�`    @�`    A�P�  �Aifj    B=0�@�@    @�@    A�Y�  �A���    B<�w@�     @�     A�$y  �A�i�    B<2A@�     @�     AÚ�  �A�X\    B;�@��    @��    A�g  �A��*    B;3�@��    @��    Ač�  �A�_�    B:��@�    @�    A@  �A�!H    B:5a@�    @�    A�4-  �A��7    B9�&@�`    @�`    A���  �A�	�    B96�@�@    @�@    A��K  �B��    B8��@��     @��     A�<!  �Bma    B88p@��     @��     A�ٓ  �Be�    B7�0@���    @���    A���  �B'n:    B79�@���    @���    A���  �B5��    B6��@�Ǡ    @�Ǡ    A���  �BG�[    B6;m@�ɀ    @�ɀ    A�e�  �B?c.    B5�*@��`    @��`    A���  �B>*�    B5<�@��@    @��@    A��  �BA�    B4��@��     @��     A��  �BB��    B4>Z@��     @��     A�9=  �BG׳    B3�@���    @���    A�%�  �BOB�    B3?�@���    @���    A�H:  �BH�g    B2��@�֠    @�֠    A���  �BGֈ    B2A7@�؀    @�؀    A�D�  �BV�    B1��@��`    @��`    A���  �BaX�    B1B�@��@    @��@    A��}  �BaG�    B0�S@��     @��     Ax�z  �Bht�    B0D@��     @��     AZq�  �Bp�6    B/ķ@���    @���    A^�x  �Bs��    B/Eg@���    @���    Adbt  �Bt�    B.�@��    @��    Aj�o  �Bzr�    B.F�@��    @��    A��4  �B}$.    B-�s@��`    @��`    A���  �B~�    B-H@��@    @��@    A�y�  �B}yc    B,��@��     @��     A��  �Bv�!    B,Iv@��     @��     A�0�  �Bs;3    B+� @���    @���    A�0�  �Bv[]    B+J�@���    @���    A�Uv  �B|��    B*�r@���    @���    Av�8  �B��    B*L@���    @���    A}AZ  �B~|�    B)��@��`    @��`    A�%~  �Bw(    B)Mf@��@    @��@    A���  �Bnh`    B(�@��     @��     A��  �BfO�    B(N�@��     @��     A��/  �B\�    B'�S@���    @���    A�~�  �BT}�    B'O�@��    @��    A���  �BIhf    B&З@��    @��    A�|�  �BAG`    B&Q8@��    @��    A�ڷ  �B>`]    B%��@�`    @�`    A��  �B7�    B%Rx@�	@    @�	@    A��  �B1�    B$�@�     @�     A�yI  �B/Ƀ    B$S�@�     @�     A�i�  �B,8:    B#�R@��    @��    B��  �B# j    B#T�@��    @��    B�\  �B7X    B"Պ@��    @��    B
Ǎ  �B	�    B"V%@��    @��    B9�  �B-�    B!ֿ@�`    @�`    BYf  �B�m    B!WY@�@    @�@    B��  �A��J    B ��@�     @�     BZ�  �A��L    B X�@�     @�     B\�  �A��    B�!@��    @��    B��  �A�iO    BY�@��    @��    B�V  �A�X�    B�M@�!�    @�!�    BH  �A�E�    BZ�@�#�    @�#�    B��  �A��    B�w@�%`    @�%`    B	�  �A��+    B\@�'@    @�'@    B�  �A�<�    Bܞ@�)     @�)     B �  �A�(�    B]0@�+     @�+     A�?�  �A�gy    B��@�,�    @�,�    A�Q�  �A�M    B^S@�.�    @�.�    A�~  �A�y    B��@�0�    @�0�    A�@q �A���    B_s@�2�    @�2�    A��� �B��    B�@�4`    @�4`    A�� �B8�    B`�@�6@    @�6@    A��X �B
�    B�@�8     @�8     A�"�  �B��    Ba�@�:     @�:     A�}�  �B��    B�8@�;�    @�;�    A���  �B�X    Bb�@�=�    @�=�    A�T  �B�~    B�P@�?�    @�?�    A�.�  �B��    Bc�@�A�    @�A�    A�=?  �B!!�    B�d@�C`    @�C`    A��  �B)9�    Bd�@�E@    @�E@    A�*z  �B>D!    B�w@�G     @�G     AՂ�  �BB+~    Be�@�I     @�I     A�  �BA$�    B�@�J�    @�J�    A�I�  �B@�    Bg@�L�    @�L�    A��Z  �B=%W    B�@�N�    @�N�    AՂ�  �B9%�    Bh@�P�    @�P�    A�͕  �B48    B�@�R`    @�R`    Aኛ  �B+��    Bi$@�T@    @�T@    A��  �B%_�    B�@�V     @�V     A���  �B    Bj+@�X     @�X     B�  �B��    B�@�Y�    @�Y�    B��  �B-    Bk1@�[�    @�[�    B�$  �B�    B�@�]�    @�]�    B3j  �B.    Bl4@�_�    @�_�    B�l  �A���    B�@�a`    @�a`    B
ں �A�    Bm5@�c@    @�c@    B�S �A��    B��@�e     @�e     Ba� �A�B    Bn3@�g     @�g     B)  �A޺�    B�@�h�    @�h�    Bh� �Aێ�    Bo0@�j�    @�j�    B�* �A؄�    B
�@�l�    @�l�    B �A��    B
p*@�n�    @�n�    A�a� �AѰ$    B	�@�p`    @�p`    B� �A��    B	q#@�r@    @�r@    Bk. �A��    B�@�t     @�t     B�� �Aʝ�    Br@�v     @�v     BE6 �A��    B�@�w�    @�w�    B� �A�A�    Bs@�y�    @�y�    B� �A���    B�@�{�    @�{�    A��� �A��    Bt @�}�    @�}�    A�#R �Aĕ    B�x@�`    @�`    A��� �A��^    Bt�@�@    @�@    A��5 �Aĸa    B�h@�     @�     A��4 �A���    Bu�@�     @�     A�g( �Aǐ�    B�U@��    @��    A�8 �A��    Bv�@��    @��    A뼭 �A��.    B�A@�    @�    A�{� �A�(�    Bw�@�    @�    A�h� �A�3j    B�*@�`    @�`    A�
� �A���    Bx�@�@    @�@    A�`? �A��<    B �@�     @�     A�� �A��    B y�@�     @�     B�p �AϜ�    A���@��    @��    BZG �A�At    A���@��    @��    B c� �A�\�    A���@�    @�    B+
 �AͲ�    A���@�    @�    B9 �A�=�    A��}@�`    @�`    B�� �AжL    A��^@�@    @�@    B	�� �A�Q�    A��>@�     @�     B
 �A�~    A��@�     @�     BӁ �A�0    A���@��    @��    Bi� �AնU    A���@��    @��    Bv� �A�4�    A���@�    @�    B�S �A�E    A���@�    @�    B~� �A�+�    A��m@�`    @�`    B�� �A���    A��G@�@    @�@    B& �A�g    A�  @�     @�     B	{� �A�=    A� �@�     @�     B
�	 �A��    A��@��    @��    B	=y �A�c'    A��@��    @��    B�B �A���    A�}@�    @�    B�� �A�p�    A�R@�    @�    B�� �A�
�    A�'@�`    @�`    B	:� �A�;c    A��@�@    @�@    B�+ �A���    A��@�     @�     B  �B 5    A��@��     @��     Bl�  �B�    A�q@���    @���    B1�  �A��Q    A�	A@���    @���    B
<  �A��=    A�
@�Ơ    @�Ơ    BE  �A��    A�
�@�Ȁ    @�Ȁ    B	�  �A�e�    A��@��`    @��`    B
$}  �A���    A�|@��@    @��@    B͜  �A���    A�H@��     @��     Bo  �BD    A�@��     @��     A�  �B&    A��@���    @���    A��N  �B��    A��@���    @���    A�1�  �B?-    A�t@�ՠ    @�ՠ    A���  �A��h    A�=@�׀    @�׀    A��l �A��)    A�@��`    @��`    A�7W �A    A��@��@    @��@    B �> �A�x�    A��@��     @��     B� �A�M�    A�Z@��     @��     A��� �A���    A�@���    @���    A��  �A�}�    A��@���    @���    A�9  �A�#�    A��@��    @��    A�D  �A�;4    A�k@��    @��    A�� �A��B    A�.@��`    @��`    A�Α  �A��g    A��@��@    @��@    A��5 �A�I�    A��@��     @��     A�ad �A�+1    A�r@��     @��     A�V� �AԖ�    A�2@���    @���    A�'
 �A��    A��@���    @���    A�/� �A��    A��@��    @��    A�� �Aϔl    A�n@���    @���    A�� �A̴    A�+@��`    @��`    A�& �A���    A��@��@    @��@    A�E	 �AŜc    A��@��     @��     A�X^ �Aā�    A� _@��     @��     A��� �A��t    A�!@���    @���    A鴘 �A���    A�!�@� �    @� �    A�i� �A��    A�"�@��    @��    A�, �A��[    A�#F@��    @��    A�!` �A��.    A�#�@�`    @�`    A�`q �A�V�    A�$�@�@    @�@    A⿆ �A��E    A�%m@�
     @�
     A�o �A���    A�&#@�     @�     A��� �A��/    A�&�@��    @��    A�%� �A�~(    A�'�@��    @��    A�QK �A���    A�(C@��    @��    A� �A�|�    A�(�@��    @��    A��� �A��@    A�)�@�`    @�`    A�V �A���    A�*]@�@    @�@    Aϣ� �A�2G    A�+@�     @�     A�o9 �A�ߩ    A�+�@�     @�     A��7 �Al�8    A�,r@��    @��    A�� �Ad�$    A�-"@��    @��    A��o �A^��    A�-�@� �    @� �    A��� �AZ.�    A�.�@�"�    @�"�    A�M� �A^)_    A�/0@�$`    @�$`    A�_ �AX�    A�/�@�&@    @�&@    A��L �AQ̠    A�0�@�(     @�(     A��� �AH(�    A�18@�*     @�*     A�QM �AC.�    A�1�@�+�    @�+�    A��� �A>�5    A�2�@�-�    @�-�    A�Ή �A<��    A�3<@�/�    @�/�    A��& �A:�    A�3�@�1�    @�1�    A��� �A;�v    A�4�@�3`    @�3`    A�]x �A8h�    A�5;@�5@    @�5@    A�B �A0L�    A�5�@�7     @�7     A�� �A0F�    A�6�@�9     @�9     A�� �A,�&    A�75@�:�    @�:�    A� �A/    A�7�@�<�    @�<�    A�� �A1��    A�8�@�>�    @�>�    A�]~ �A7n�    A�9*@�@�    @�@�    A��z �A9��    A�9�@�B`    @�B`    A�d� �A3��    A�:v@�D@    @�D@    A��N �A+��    A�;@�F     @�F     A��e �A#RU    A�;�@�H     @�H     A�,� �A��    A�<c@�I�    @�I�    A�C� �A8    A�=@�K�    @�K�    A��� �A�    A�=�@�M�    @�M�    A�[� �A�    A�>L@�O�    @�O�    A��� �A�D    A�>�@�Q`    @�Q`    A��� �A �}    A�?�@�S@    @�S@    A��< �@��o    A�@0@�U     @�U     A��� �@�C�    A�@�@�W     @�W     A�BI �@�M%    A�Aq@�X�    @�X�    A��� �@���    A�B@�Z�    @�Z�    A��I �@�w�    A�B�@�\�    @�\�    A��� �@��Q    A�CN@�^�    @�^�    A��� �@�U�    A�C�@�``    @�``    A�-� �@��s    A�D�@�b@    @�b@    A�= �@�>    A�E'@�d     @�d     A�!� �@�&    A�E�@�f     @�f     A��� �@��x    A�F`@�g�    @�g�    A�p� �@�7�    A�F�@�i�    @�i�    A�5� �@��    A�G�@�k�    @�k�    A��� �@�s�    A�H2@�m�    @�m�    A��w �@�v,    A�H�@�o`    @�o`    A��� �@�;8    A�If@�q@    @�q@    A�Ћ �@�N�    A�J @�s     @�s     A��� �@�w�    A�J�@�u     @�u     A��F �@Φ�    A�K2@�v�    @�v�    A��� �@ʟ    A�K�@�x�    @�x�    A�C� �@�y�    A�Lb@�z�    @�z�    A�� �@�z`    A�L�@�|�    @�|�    A��� �@���    A�M�@�~`    @�~`    A��� �@�<�    A�N&@�@    @�@    A��� �@��o    A�N�@�     @�     A��B �@��    A�OR@�     @�     A�g` �@��    A�O�@��    @��    A��V �@�P    A~��@��    @��    A�!{ �@�4�    A|� @�    @�    A�l� �@�S3    Az�H@�    @�    A�o� �@��    Ax�o@�`    @�`    Ay�% �@���    Av��@�@    @�@    Ap�F �@�c�    At��@�     @�     AtI; �@��    Ar��@�     @�     A^�n �@�h�    Ap�@��    @��    AQ�A �@���    An�&@��    @��    AH�� �@�L    Al�H@�    @�    AB�� �@�?�    Aj�i@�    @�    A:�= �@�ݨ    Ah��@�`    @�`    A8�� �@��    Af��@�@    @�@    A4 �@��    Ad��@�     @�     A,p� �@��    Ab��@�     @�     A-P� �@��    A`�@��    @��    A1Q� �@�cf    A^�!@��    @��    A.�� �@�bL    A\�=@�    @�    A+eb �@�48    AZ�X@�    @�    A-�� �@�K    AX�r@�`    @�`    A.^� �@�z    AV��@�@    @�@    A0� �@��    AT��@�     @�     A5�� �@�L    AR��@�     @�     A5�� �@��    AP��@��    @��    A8� �A��    AN��@��    @��    A<�2 �@�VV    AL�@�    @�    A?� �@���    AJ�@�    @�    A;� �@���    AH�*@�`    @�`    A9�� �A0    AF�>@�@    @�@    A2� �A�^    AD�Q@�     @�     A+�B �A�[    AB�c@��     @��     A'f� �A]�    A@�u@���    @���    A%%x �A�    A>Ć@���    @���    A'�J �A|�    A<Ŗ@�Š    @�Š    A&R1 �A��    A:ƥ@�ǀ    @�ǀ    A&�o �A��    A8Ǵ@��`    @��`    A(DO �AQ    A6��@��@    @��@    A/ �A�t    A4��@��     @��     A,"8 �A `    A2��@��     @��     A+�. �A&*Q    A0��@���    @���    A+� �A'f�    A.��@���    @���    A%S �A(�x    A,��@�Ԡ    @�Ԡ    A!Se �A+��    A*�@�ր    @�ր    A#� �A7�Y    A(�@��`    @��`    Ax �A/�    A&�@��@    @��@    A�x �A2,�    A$�"@��     @��     A�Z �A4F    A"�)@��     @��     A� �AG!�    A �0@���    @���    A�� �AB>    A�6@���    @���    A�/ �AFJ�    A�;@��    @��    A�� �AE3    A�@@��    @��    A�M �A@H    A�D@��`    @��`    AX� �A@|R    A�H@��@    @��@    Afm �AD��    A�J@��     @��     A
 �ALх    A�M@��     @��     Ar, �AH{�    A�N@���    @���    A	#X �AG?    A�O@���    @���    A�j �AG�,    A�O@��    @��    A�] �AK8�    A
�O@��    @��    A� �AIF�    A�N@��`    @��`    A�6 �AI��    A�M@��@    @��@    A2. �AF}    A�K@��     @��     AF� �AC
�    A�H@��     @��     A�� �A>�    A �E@���    @���    A!� �A@�    @�ʁ@���    @���    A"�� �A=�1    @��x@��    @��    A%�� �A;k�    @��n@��    @��    A)�q �A5�[    @��c@�`    @�`    A*[: �A*�t    @��V@�@    @�@    A-�X �A#L�    @��H@�	     @�	     A.�, �A1�    @��:@�     @�     A0� �ACf    @��*@��    @��    A)�0 �Aam    @��@��    @��    A*F: �A��    @��@��    @��    A*�� �@���    @���@��    @��    A(H �@��    @���@�`    @�`    A&�� �@�W    @���@�@    @�@    A˞ �@إK    @��@�     @�     A�� �@�I2    @��@�     @�     A�� �@�(    @��@��    @��    A �@�    @��g@��    @��    AX� �@�g�    @��L@��    @��    A� �@��    @��/@�!�    @�!�    @�(� �@��    @��@�#`    @�#`    @�l� �@�:    @���@�%@    @�%@    @׷ �@�.    @���@�'     @�'     @�M �@�    @���@�)     @�)     @��� �@���    @���@�*�    @�*�    @�o� �@�M�    @��o@�,�    @�,�    @�� �@�L-    @��L@�.�    @�.�    @ҬC �@�l^    @��'@�0�    @�0�    @�� �@���    @��@�2`    @�2`    @�1D �@�۶    @���@�4@    @�4@    @�_� �@��    @��@�6     @�6     @�lt �A8�    @��@�8     @�8     @ׅb��A��    @�`@�9�    @�9�    @����A'�    @|j@�;�    @�;�    @�
K��A2 !    @t@�=�    @�=�    A	:+��A++�    @l�@�?�    @�?�    A�D��A+1�    @d]@�A`    @�A`    A����A \�    @\�@�C@    @�C@    A$����A��    @T �@�E     @�E     Ao���A9i    @L$>@�G     @�G     A����@���    @D'�@�H�    @�H�    ��  ����      @<+u@�J�    @�J�    ��  ����      @4/@�L�    @�L�    ��  ����      @,2�@�N�    @�N�    ��  ����      @$6:@�P`    @�P`    ��  ����      @9�@�R@    @�R@    ��  ����      @=^@�T     @�T     ��  ����      @@�@�V     @�V     ��  ����      @D{@�W�    @�W�    ��  ����      ?��@�Y�    @�Y�    ��  ����      ?� @�[�    @�[�    ��  ����      ?؞/@�]�    @�]�    ��  ����      ?ȥ;@�_`    @�_`    ��  ����      ?��D@�a@    @�a@    ��  ����      ?��I@�c     @�c     ��  ����      ?��J@�e     @�e     ��  ����      ?��H@�f�    @�f�    ��  ����      ?q��@�h�    @�h�    ��  ����      ?Q�r@�j�    @�j�    ��  ����      ?1�Y@�l�    @�l�    ��  ����      ?�8@�n`    @�n`    ��  ����      >�!@�p@    @�p@    ��  ����      >���@�r     @�r     ��  ����      >G��@�t     @�t     ��  ����      =���