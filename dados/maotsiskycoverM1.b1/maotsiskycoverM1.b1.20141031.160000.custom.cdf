CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.3-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20141031160000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      63 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      238 pixels     center_y      318 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-10-31 23:49:02, using ingest-tsi-12.3-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20141031.160000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:16 2017 UTC       ANDERS_armtime_timestamp      1488830656     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-10-31 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-10-31 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-10-31 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �TR� BH  �rdt�M�M@�      @�      B;�%  �A�dc    B��<@�#�    @�#�    B9�  �A���    B��@�'�    @�'�    B7�g  �A�    B��/@�+@    @�+@    B/ܶ  �Aϳv    B���@�/     @�/     B�  �A׭�    B���@�2�    @�2�    B��  �Aܯy    B��@�6�    @�6�    B�a  �A�     B�o�@�:@    @�:@    B)�A  �A���    B�Y @�>     @�>     B*ez  �A��@    B�A�@�A�    @�A�    B�  �A|u    B�)�@�E�    @�E�    B9h  �Aqӣ    B�d@�I@    @�I@    B�  �Am��    B��^@�M     @�M     Be  �A~^_    B���@�P�    @�P�    A��   �AX�    B�ı@�T�    @�T�    B|  �AZ��    B��@�X@    @�X@    B4H  �AX��    B���@�\     @�\     B�1  �AO?�    B�s8@�_�    @�_�    B��  �AF��    B�W@�c�    @�c�    B	�M  �AD�    B�:W@�g@    @�g@    B��  �A;�v    B�'@�k     @�k     B�� �A9uQ    B��z@�n�    @�n�    B��  �A@��    B��P@�r�    @�r�    A��� �AD�    B�­@�v@    @�v@    A�c� �AB
,    B���@�z     @�z     A�S  �A,�    B�� @�}�    @�}�    B�� �A2X    B�c�@쁀    @쁀    B
[D �A�a    B�C@�@    @�@    B� �AX�    B�"�@�     @�     Bm9  �A��    B�8@��    @��    BwG  �Aa~    B��o@쐀    @쐀    Bg�  �A'�    B��:@�@    @�@    BKY �A)�    B���@�     @�     B
Z �A+u�    B�w�@��    @��    Bh�  �A%�2    B�T@쟀    @쟀    B� �A��    B�0H@�@    @�@    A��% �A��    B�@�     @�     A�ȼ  �A%    B��o@��    @��    A�+Q  �Ac^    B��p@쮀    @쮀    A��  �@���    B��@�@    @�@    A�� �@�_    B�wU@�     @�     A�>�  �@��    B�Q<@��    @��    A��  �@���    B�*�@콀    @콀    A�f) �A��    B��@��@    @��@    A��� �@��(    B���@��     @��     A�ڍ �@���    B��V@���    @���    A�`r  �@��    B���@�̀    @�̀    A٣ �@�Q1    B�e]@��@    @��@    A��h �@ْ�    B�<�@��     @��     A͋�  �@��    B�@���    @���    A�s  �@��}    B���@�ۀ    @�ۀ    A���  �@���    B���@��@    @��@    A���  �@樹    B���@��     @��     A�� �@�Ǧ    B�m�@���    @���    A�U  �@٘�    B�C�@��    @��    A��^  �@�C�    B��@��@    @��@    A���  �@�D�    B��@��     @��     A��  �@��    B���@���    @���    A��%  �@�Э    B���@���    @���    A�?4  �@�?�    B�k�@��@    @��@    A�҉  �@��	    B�?�@�     @�     A��C  �@���    B��@��    @��    A�-�  �@�KZ    B��;@��    @��    A��  �@�$�    B���@�@    @�@    A��W  �@��    B���@�     @�     A��_  �@�7�    B�`W@��    @��    A�~C  �@�j�    B�2�@��    @��    A���  �@�a�    B�6@�@    @�@    A���  �@�_    B��M@�     @�     A��v  �@�4�    B��)@�"�    @�"�    A��  �@�Y�    B�z�@�&�    @�&�    A���  �@��/    B�L8@�*@    @�*@    A��  �@�m�    B�l@�.     @�.     A���  �@ϵ*    B��i@�1�    @�1�    A��5  �@��    B��2@�5�    @�5�    A���  �@؟�    B���@�9@    @�9@    A÷�  �@�n    B�`%@�=     @�=     A�*�  �@�6�    B�0S@�@�    @�@�    A�^�  �@���    B� N@�D�    @�D�    A�y�  �@�7�    B��@�H@    @�H@    A��  �@�(�    B���@�L     @�L     A�S�  �@���    B�o@�O�    @�O�    A��  �@�o.    B�>X@�S�    @�S�    A�C�  �@�r    B�e@�W@    @�W@    A���  �@�5�    B��E@�[     @�[     A� <  �@�Sx    B���@�^�    @�^�    A�l  �@�    B�y�@�b�    @�b�    AѸ�  �@��@    B�G�@�f@    @�f@    Ã�  �@�    B�@�j     @�j     A���  �@�m�    B��@�m�    @�m�    A�7)  �@�-�    B���@�q�    @�q�    A�.  �@��    B��@�u@    @�u@    A���  �@�    B�M@@�y     @�y     A���  �@�    B��@�|�    @�|�    A��   �@�I�    B���@퀀    @퀀    A��R  �@���    B��@�@    @�@    A�*�  �A �	    B���@�     @�     A�U  �A8�    B�N�@��    @��    A�$�  �A��    B��@폀    @폀    A�!S  �A�S    B��@�@    @�@    A�ک  �AZ    B��t@�     @�     AɎ�  �A�X    B���@��    @��    A�j�  �A2M    B�L�@힀    @힀    AΆo  �AWa    B��@��@    @��@    A���  �A'Ӗ    B���@��     @��     A��  �A..&    B���@���    @���    A��%  �A7��    B�|&@���    @���    Aؑ�  �A@1Y    B�G�@��@    @��@    A�*�  �AN3�    B�@��     @��     A���  �ALk�    B��R@���    @���    A�i�  �AF��    B��z@���    @���    A�[�  �A0�    B�t�@��@    @��@    A��O  �A(by    B�?r@��     @��     A�K  �A-�Z    B�
B@���    @���    A��r  �A$�    B���@�ˀ    @�ˀ    A�   �A ��    B���@��@    @��@    A���  �A��    B�j@��     @��     A�@  �A2    B�4r@���    @���    A��  �Ai�    B���@�ڀ    @�ڀ    A��  �Ao�    B���@��@    @��@    AϣL  �A�    B���@��     @��     Aý  �Ah    B�\�@���    @���    A��-  �A�    B�&�@��    @��    A�a�  �A    B��@��@    @��@    A�ڤ  �A[    B��Z@��     @��     Aċ:  �A ��    B���@���    @���    A��|  �@���    B�Mx@���    @���    A��5  �@�.    B��@��@    @��@    A� �  �@��    B��:@�      @�      A؁C  �@�t    B��y@��    @��    A�d  �@�ؿ    B�r�@��    @��    A�  �@�&}    B�;�@�@    @�@    A��  �@�]�    B��@�     @�     A�ʦ  �A V    B�͛@��    @��    A�:   �A�_    B��n@��    @��    B�O  �A{�    B�_-@�@    @�@    B�  �A%E    B�'�@�     @�     B)g  �Av�    B��o@�!�    @�!�    B{H  �An    B���@�%�    @�%�    Bh�  �A@I    B��b@�)@    @�)@    B4��  �AM�    B�I�@�-     @�-     B!2�  �A=$7    B�@�0�    @�0�    B�  �A�(o    B��?@�4�    @�4�    B�  �Ao�    B��c@�8@    @�8@    B�k  �A`    B�jv@�<     @�<     B�  �AY>)    B�2w@�?�    @�?�    Bvh  �AH�N    B��f@�C�    @�C�    B��  �AF��    B��C@�G@    @�G@    B�W  �A<    B��@�K     @�K     B��  �A>�!    B�Q�@�N�    @�N�    B ��  �A:��    B�w@�R�    @�R�    A�  �@��    B��@�V@    @�V@    A�Ȣ  �@ޑ+    B���@�Z     @�Z     A�yf  �A �    B�p@�]�    @�]�    A�{�  �@�3b    B�7�@�a�    @�a�    A��2  �@�=�    B���@�e@    @�e@    A���  �@�(5    B��(@�i     @�i     A��  �@�;�    B��d@�l�    @�l�    A��  �@��"    B�T�@�p�    @�p�    A���  �@�K    B��@�t@    @�t@    A��`  �@��    B���@�x     @�x     A��  �@�h    B���@�{�    @�{�    A�|�  �@�1�    B�p�@��    @��    A��!  �@˄�    B�7�@�@    @�@    A�  �@��    B��@�     @�     A�  �@�h    B�{@��    @��    A���  �@鴊    B�@    @    Aյ�  �@�8�    B~�T@�@    @�@    A�ID  �A
�    B~2�@�     @�     A��  �A`�    B}��@��    @��    A��  �A��    B}L�@    @    A���  �A(��    B|��@�@    @�@    A�͵  �A-�V    B|f�@�     @�     A�B�  �A@8{    B{�q@��    @��    B,  �A���    B{�@    @    B:S  �A�x�    B{�@�@    @�@    Bx#  �A��     Bz�4@�     @�     Bʳ  �A��j    Bz%�@��    @��    B�<  �A~��    By��@    @    B��  �A*��    By>#@�@    @�@    BL  �A4K^    Bx�F@��     @��     B7�  �A4    BxVR@���    @���    B�  �A�c    Bw�H@�ʀ    @�ʀ    Bb�  �A.?�    Bwn(@��@    @��@    B�b �AAC    Bv��@��     @��     BL� �AK}�    Bv��@���    @���    B�� �AO4�    BvJ@�ـ    @�ـ    B#�L  �A��    Bu��@��@    @��@    B�%  �A��    Bu(N@��     @��     B94  �Aʊ    Bt��@���    @���    B  �A؋�    Bt?@��    @��    A���  �A��i    Bs�=@��@    @��@    A�y�  �A�F�    BsUe@��     @��     A�ł  �A�+�    Br�z@���    @���    A��  �Aԍz    Brk}@���    @���    A�O/ �A�Ue    Bq�l@��@    @��@    A���  �A�W#    Bq�I@��     @��     A��i  �A�$�    Bq@��    @��    A�|u  �A�.    Bp��@��    @��    A��`  �A��q    Bp!r@�
@    @�
@    Bp�  �A�    Bo�@�     @�     B�z  �A��    Bo6�@��    @��    B��  �A�9    Bn��@��    @��    B	n�  �A�R�    BnK\@�@    @�@    B��  �A���    Bmլ@�     @�     B��  �Ai��    Bm_�@� �    @� �    B� �A`�w    Bl�@�$�    @�$�    B �A]��    Blt8@�(@    @�(@    B�$ �AX�-    Bk�F@�,     @�,     B� �AUز    Bk�D@�/�    @�/�    B�� �A_��    Bk3@�3�    @�3�    B��  �A\�|    Bj�@�7@    @�7@    Bs�  �Ah�A    Bj%�@�;     @�;     B�  �Ao    Bi��@�>�    @�>�    B��  �A� M    Bi9R@�B�    @�B�    Ba�  �A��D    Bh��@�F@    @�F@    A�jS  �A�X�    BhL�@�J     @�J     A��#  �A�:1    Bg�@�M�    @�M�    A���  �A���    Bg_�@�Q�    @�Q�    Bd'  �AxO�    Bf��@�U@    @�U@    BH�  �A�    BfrD@�Y     @�Y     B3�  �A���    Be��@�\�    @�\�    B�r  �A���    Be��@�`�    @�`�    Bb[  �A� �    Be�@�d@    @�d@    Bm  �A��    Bd�"@�h     @�h     B .�  �A�_    Bd 8@�k�    @�k�    B#��  �A���    Bc�@@�o�    @�o�    B#[�  �A�iM    Bc2;@�s@    @�s@    B�'  �A�i�    Bb�*@�w     @�w     B  �AФ=    BbD@�z�    @�z�    B!��  �A�    Ba��@�~�    @�~�    BAT  �A�.�    BaU�@�@    @�@    B�.  �A�Io    B`�e@�     @�     B�X  �A��C    B`g@��    @��    B�<  �A�`U    B_�@    @    B��  �A��1    B_xP@�@    @�@    B {  �A��!    B_ �@�     @�     B"��  �A�s�    B^�]@��    @��    Bb�  �A��@    B^�@    @    BV!  �A���    B]�:@�@    @�@    B��  �A�]|    B]"�@�     @�     B��  �A�    B\��@��    @��    B.  �A��    B\31@變    @變    B�=  �A�\    B[�n@�@    @�@    B�  �AЋ�    B[C�@�     @�     B�o  �B��    BZ��@��    @��    Bܜ  �B
��    BZS�@ﺀ    @ﺀ    B#Q  �B��    BY��@�@    @�@    B�}  �B�*    BYc�@��     @��     Bߺ  �B��    BX��@���    @���    B��  �B��    BXs�@�ɀ    @�ɀ    B*�  �B�    BW��@��@    @��@    B%{  �A��    BW��@��     @��     B"6  �A��
    BW�@���    @���    B#Ĭ  �A�%_    BV�M@�؀    @�؀    B"�  �A���    BV@��@    @��@    B k�  �A�'�    BU��@��     @��     BՀ  �A���    BU*q@���    @���    B3t  �A�֏    BT�@��    @��    BBM  �A��"    BT9�@��@    @��@    B�u  �A��A    BS�A@��     @��     B�f  �A��&    BSH�@���    @���    B
�T  �A�ԟ    BR�H@���    @���    A��  �AǪ	    BRW�@��@    @��@    B�,  �A��    BQ�+@��     @��     B��  �A��V    BQf�@� �    @� �    B��  �A�t�    BP��@��    @��    BJ�  �A���    BPu=@��    @��    B�  �A�J    BO��@��    @��    BB"  �B#i    BO��@�`    @�`    B"  �B�    BO@�
@    @�
@    B�M  �B��    BN�4@�     @�     B��  �B ��    BN]@�     @�     B��  �A�     BM�~@��    @��    B�W  �A���    BM'�@��    @��    B�  �A��    BL��@��    @��    B\�  �A�CC    BL5�@��    @��    B��  �A뎫    BK��@�`    @�`    B@   �A�w�    BKC�@�@    @�@    B�V  �A�LG    BJʛ@�     @�     B�R  �A��    BJQ�@�     @�     A���  �A�)�    BI�h@��    @��    A�LI  �A�P�    BI_B@� �    @� �    A���  �A���    BH�@�"�    @�"�    A�`�  �A�Q�    BHl�@�$�    @�$�    A�/Y  �A�R�    BG�@�&`    @�&`    A�F�  �A��j    BGze@�(@    @�(@    A�n�  �A��V    BG@�*     @�*     A���  �A�D    BF��@�,     @�,     A�nc  �Aћ    BFs@�-�    @�-�    A�Z%  �A��    BE�@�/�    @�/�    A��0 �A��3    BE�@�1�    @�1�    A�i�  �A���    BD�C@�3�    @�3�    A�.�  �A�3|    BD(�@�5`    @�5`    A��  �A~)f    BC�W@�7@    @�7@    A��  �A��    BC5�@�9     @�9     A߱   �A�,E    BB�Q@�;     @�;     A��D  �A�q�    BBB�@�<�    @�<�    A�w�  �A��w    BA�0@�>�    @�>�    A�z8  �A�8�    BAO�@�@�    @�@�    A�-s  �A�g    B@��@�B�    @�B�    A��n  �A��@    B@\P@�D`    @�D`    A��&  �A�B�    B?�@�F@    @�F@    Aܧo  �Aܳu    B?h�@�H     @�H     A�&(  �A���    B>�8@�J     @�J     A�^  �A�/�    B>uy@�K�    @�K�    A�.�  �A�h@    B=��@�M�    @�M�    A��  �A�    B=��@�O�    @�O�    A耩  �A�U�    B=@�Q�    @�Q�    A��[  �A��
    B<�B@�S`    @�S`    A�j  �A�LQ    B<f@�U@    @�U@    A�)  �A�l    B;��@�W     @�W     Aل4  �A�N�    B; �@�Y     @�Y     A��$  �A�	k    B:��@�Z�    @�Z�    A�'�  �A�ar    B:,�@�\�    @�\�    A�+D  �A��u    B9��@�^�    @�^�    A�Г  �A�@�    B98�@�`�    @�`�    A��@  �A�L(    B8��@�b`    @�b`    Aݻ:  �A��K    B8D�@�d@    @�d@    A�v  �A���    B7ʯ@�f     @�f     A�E�  �A{�/    B7P�@�h     @�h     Aۤ�  �A|��    B6օ@�i�    @�i�    A�[�  �A���    B6\g@�k�    @�k�    A�)E  �A�    B5�E@�m�    @�m�    A��=  �A���    B5h@�o�    @�o�    A���  �A�qk    B4��@�q`    @�q`    A��G �A�#q    B4s�@�s@    @�s@    A��� �A~�W    B3��@�u     @�u     A��� �Ao��    B3O@�w     @�w     A��p �Ab4U    B3@�x�    @�x�    Ap�A �AX�    B2��@�z�    @�z�    A_�x  �Aat    B2�@�|�    @�|�    AY#u  �A���    B1�2@�~�    @�~�    Ac1�  �A�͈    B1�@��`    @��`    Ac��  �A��E    B0��@��@    @��@    AIw  �A���    B0',@��     @��     ABE�  �A��F    B/��@��     @��     A=Nv  �A���    B/2f@���    @���    A>��  �A��    B.��@���    @���    AJ�  �A�/)    B.=�@���    @���    AT�  �A���    B-�@���    @���    AjBB  �A��3    B-H�@��`    @��`    Ak��  �@۠�    B,�*@�@    @�@    Au��  �@ۚ�    B,S�@�     @�     Ax�  �@�*�    B+�'@�     @�     Ay�  �@�'�    B+^�@��    @��    At M  �@�    B*�@��    @��    Au,�  �@�bt    B*i�@�    @�    Ar��  �@�MZ    B)��@�    @�    Ak?u  �@���    B)tW@�`    @�`    Aj�  �A ��    B(��@�@    @�@    AlM�  �A ��    B(@�     @�     Aqz{  �Am
    B(v@�     @�     Ax��  �A�s    B'��@��    @��    A��  �A��    B'"@��    @��    A�A~  �A(?    B&�r@�    @�    A��!  �A�&    B&�@�    @�    A�Ӳ  �A��    B%�@�`    @�`    A��2  �A�    B%$L@�@    @�@    A���  �AM    B$��@�     @�     A��F  �A�Z    B$.�@�     @�     A��Y  �A#ղ    B#�@��    @��    A��(  �A)�Y    B#9:@��    @��    A�5�  �A$��    B"�m@�    @�    A��  �A!>�    B"C�@�    @�    A�Ͽ  �A��    B!��@�`    @�`    A�1g  �A%{    B!M�@�@    @�@    A��  �A)p�    B �@��     @��     A��L  �A.Y�    B X5@��     @��     A�У  �A7V    B�S@���    @���    A�;7  �ABNm    Bbm@���    @���    A�3�  �AJg    B�@�Ǡ    @�Ǡ    A��U  �AXӇ    Bl�@�ɀ    @�ɀ    A�$}  �A\�b    B�@��`    @��`    Aġ  �Ag�-    Bv�@��@    @��@    A�m�  �Ap[    B��@��     @��     A�L  �A���    B��@��     @��     A��  �A���    B�@���    @���    A��N  �A�C    B��@���    @���    A�)�  �A�թ    B�@�֠    @�֠    A��  �A��C    B��@�؀    @�؀    A�A  �A�C    B�@��`    @��`    A�hz  �A��F    B��@��@    @��@    A�TH  �A��B    B#�@��     @��     B(�  �A�Ħ    B��@��     @��     B{:  �A��r    B-x@���    @���    B�o  �A��_    B�`@���    @���    BF�  �A��B    B7E@��    @��    Bc  �A�|"    B�(@��    @��    B'�  �A���    BA@��`    @��`    B��  �Aʕ�    B��@��@    @��@    BMv  �A�y�    BJ�@��     @��     A���  �A˓�    Bϕ@��     @��     A̿�  �Aګ    BTj@���    @���    AД�  �A��    B�;@���    @���    A���  �A��    B^@���    @���    A���  �A���    B��@���    @���    A���  �A�~H    Bg�@��`    @��`    A�G�  �A�J�    B�h@��@    @��@    A�?4  �B�-    Bq,@��     @��     A�o  �Bk    B��@��     @��     A��]  �B �    Bz�@���    @���    A�'�  �A���    B�j@��    @��    A���  �A�+    B�$@��    @��    A��  �A�i    B�@��    @��    Aׁ  �A�e�    B��@�`    @�`    A�b�  �A�0c    BD@�	@    @�	@    A�55  �AȌ�    B��@�     @�     A�n�  �A��    B�@�     @�     A�~,  �A��    B�N@��    @��    A�
�  �A�"�    B$�@��    @��    Aݏ� �A�@�    B��@��    @��    A�� �A�f/    B.B@��    @��    A�� �A�}_    B
��@�`    @�`    A�X� �A���    B
7�@�@    @�@    A��> �Ax��    B	�"@�     @�     Aٜ �Aj��    B	@�@�     @�     Aӊ� �Aip    B�W@��    @��    AΝ� �A|�
    BI�@��    @��    Aѩ@ �A�Ek    B΃@�!�    @�!�    A�8e �A�ń    BS@�#�    @�#�    A̓ �A��-    Bצ@�%`    @�%`    A�i% �A�2�    B\4@�'@    @�'@    A�t� �A���    B��@�)     @�)     A�.F �A�|~    BeK@�+     @�+     A�^� �A��    B��@�,�    @�,�    A���  �A�[    BnY@�.�    @�.�    AΞ   �A��    B��@�0�    @�0�    A�U  �A��    Bw`@�2�    @�2�    A�(W �A���    B��@�4`    @�4`    A��D  �A��Z    B�^@�6@    @�6@    AƬ�  �A��    B�@�8     @�8     A��B  �A�Xf    B�T@�:     @�:     A�)T  �A��    B�@�;�    @�;�    A��  �A���    B �C@�=�    @�=�    A��  �A��R    B �@�?�    @�?�    A��  �A�α    A�6T@�A�    @�A�    A�!  �A�b    A�?6@�C`    @�C`    A��m �A�l�    A�H@�E@    @�E@    Ax	n �A��    A�P�@�G     @�G     Ac�� �A��1    A�Y�@�I     @�I     A_F� �A��}    A�b�@�J�    @�J�    AZ�� �A���    A�ki@�L�    @�L�    A]�� �A��}    A�t5@�N�    @�N�    A`s� �A���    A�|�@�P�    @�P�    A];  �Au�    A���@�R`    @�R`    ASw% �Al�    A���@�T@    @�T@    AP�2 �Af��    A��D@�V     @�V     AGmB �Ab��    A� @�X     @�X     A8�8 �A[*W    A�@�Y�    @�Y�    A.�) �ASB�    A�m@�[�    @�[�    A$Y  �ALuL    A�@�]�    @�]�    A� �ABq:    A���@�_�    @�_�    AQ! �A8�2    A��x@�a`    @�a`    A�= �A7U    A�� @�c@    @�c@    A�& �A8Wx    A���@�e     @�e     AY� �A2>�    A��h@�g     @�g     A8W �A1Q�    A��@�h�    @�h�    ADe �A.A�    A���@�j�    @�j�    AN �A)�    A��<@�l�    @�l�    Ac	 �A"�a    A��@�n�    @�n�    A � �A!�&    A�e@�p`    @�p`    AH� �A'JK    A��@�r@    @�r@    AzB �A%*�    A�!�@�t     @�t     Ap� �A$d�    A�*@�v     @�v     A�D �A&	
    A�2�@�w�    @�w�    A� �A(v    A�;@�y�    @�y�    A
� �A&a�    A�C�@�{�    @�{�    A: �A'W�    A�L@�}�    @�}�    A�� �A%o    A�T�@�`    @�`    @�� �A'��    A�]@�@    @�@    @�S� �A&h�    A�e�@�     @�     @��: �A%	n    A�n@�     @�     @��; �A!�    A�vt@��    @��    @�_ �A�=    A�~�@��    @��    @�;  �A�B    AهR@�    @�    @�y4 �A��    A؏�@�    @�    @�w �AS�    Aט&@�`    @�`    @�܂ �A��    A֠�@�@    @�@    @�� �A��    Aը�@�     @�     @�� �A�|    AԱQ@�     @�     @�g �A/�    Aӹ�@��    @��    @�� �A�5    A��@��    @��    @�� �Ad�    A��h@�    @�    @�\� �A�P    A���@�    @�    @�6� �A+�    A��@�`    @�`    @ۣ �Aa    A��j@�@    @�@    @��� �AiG    A��@�     @�     @��� �A>f    A��@�     @�     @Ԃ �A��    A��Y@��    @��    @Ԕ �A��    A��@��    @��    @� � �A�Z    A��@�    @�    @�|� �A    A�5@�    @�    @�� �A3�    A�z@�`    @�`    @�]� �A_W    A�%�@�@    @�@    @��� �A��    A�-�@�     @�     @�� �AF�    A�6>@�     @�     @�N �A�    A�>|@��    @��    @�� �A�#    A�F�@��    @��    @�� �A�A    A�N�@�    @�    @�^� �A�    A�W(@�    @�    @�� �A
!�    A�_^@�`    @�`    @�'] �AO�    A�g�@�@    @�@    @�S �A�L    A�o�@�     @�     @�m �A#k    A�w�@��     @��     @��� �A��    A��$@���    @���    @�B� �A�G    A��R@���    @���    @�P� �A��    A��}@�Ơ    @�Ơ    @��  �A��    A���@�Ȁ    @�Ȁ    @��� �A�    A���@��`    @��`    @�ӈ �A�    A���@��@    @��@    @�ߌ �A��    A��@��     @��     @��� �A��    A��>@��     @��     @�� �A	��    A��`@���    @���    @��� �A��    A�Ɂ@���    @���    @�gJ �A75    A�џ@�ՠ    @�ՠ    @�"� �A��    A�ٽ@�׀    @�׀    @�M �A4�    A���@��`    @��`    @��- �A��    A���@��@    @��@    @�W� �A	&    A��@��     @��     @�� �A     A��$@��     @��     @�i� �A�B    A�;@���    @���    @��� �Ax    A�
P@���    @���    @��7 �A#    A�d@��    @��    @�� �Ae    A�v@��    @��    @��Y �A��    A�"�@��`    @��`    @�:� �A��    A�*�@��@    @��@    @�uy �A��    A�2�@��     @��     @��� �A�    A�:�@��     @��     @��� �A��    A�B�@���    @���    @�:� �Agr    A�J�@���    @���    @�� �A�    A�R�@��    @��    @��} �A��    A�Z�@���    @���    @�� �AW�    A�b�@��`    @��`    @�u� �AA    A�j�@��@    @��@    @��^ �A�I    A�r�@��     @��     @��� �A8    A�z�@��     @��     @�< �AA�    A���@���    @���    @�S� �A	f    A���@� �    @� �    @�* �AQW    A���@��    @��    @��k �A
�
    A�� @��    @��    @�� �A	p�    A���@�`    @�`    @��) �A�8    A���@�@    @�@    @�&_ �A�    A���@�
     @�
     @��o �A,    A���@�     @�     @��6 �A	ָ    A���@��    @��    @��� �A�    A���@��    @��    @ƕ( �A��    A���@��    @��    @��p �A��    A���@��    @��    @�6 �As�    A���@�`    @�`    @�'� �A�    A���@�@    @�@    @�~e �AV�    A���@�     @�     @¿ �A�    A���@�     @�     @�l_ �A��    A��@��    @��    @ţ� �A�    A�
�@��    @��    @��} �A�    A��@� �    @� �    @�o� �A��    A��@�"�    @�"�    @�� �AE�    A�"�@�$`    @�$`    @®n �A��    A�*�@�&@    @�&@    @Ɯ� �A$    A�2t@�(     @�(     @��B �Ar3    A�:f@�*     @�*     @ŚM �A(�    A�BX@�+�    @�+�    @��� �A�    A�JI@�-�    @�-�    @�6� �A�!    A�R9@�/�    @�/�    @ą� �A�3    A�Z)@�1�    @�1�    @Ɗ� �AI�    A�b@�3`    @�3`    @�m� �A4K    A�j@�5@    @�5@    @��X �A&�    A�q�@�7     @�7     @Ю� �A��    A~��@�9     @�9     @�>� �A k    A}�@�:�    @�:�    @�h� �A ��    A{�@�<�    @�<�    @�� �A&v�    Ay#\@�>�    @�>�    @�ei �A'�q    Aw36@�@�    @�@�    @��� �A'.?    AuC@�B`    @�B`    @��H �A'OL    AsR�@�D@    @�D@    @�� �A%ȱ    Aqb�@�F     @�F     @ą� �A&�    Aor�@�H     @�H     @��> �A$W    Am�p@�I�    @�I�    @�b� �A(0�    Ak�G@�K�    @�K�    @�2� �A&v�    Ai�@�M�    @�M�    @�x� �A++�    Ag��@�O�    @�O�    @ψ< �A/Z    Ae��@�Q`    @�Q`    @�*� �A+|�    AcѠ@�S@    @�S@    @�� �A+R�    Aa�u@�U     @�U     @��n �A-9�    A_�J@�W     @�W     @��= �A+m�    A^ @�X�    @�X�    @�Z� �A+j�    A\�@�Z�    @�Z�    @�Ga �A+�    AZ �@�\�    @�\�    @ǁ �A+�    AX0�@�^�    @�^�    @�ו �A1*    AV@t@�``    @�``    @��^ �A/��    ATPI@�b@    @�b@    @�� �A.�?    AR`@�d     @�d     @�8� �A0d�    APo�@�f     @�f     @�}r �A2��    AN�@�g�    @�g�    @��� �A0�"    AL��@�i�    @�i�    @�Mc �A0y�    AJ�t@�k�    @�k�    @�"" �A1�\    AH�J@�m�    @�m�    @�� �A0��    AF�!@�o`    @�o`    @�V� �A3T    AD��@�q@    @�q@    @�&� �A4U    AB��@�s     @�s     @ԲH �A9Ȟ    A@�@�u     @�u     @џ �A9)_    A>�~@�v�    @�v�    @�>� �A=�;    A=W@�x�    @�x�    @Ж� �A>ș    A;0@�z�    @�z�    @�R" �A?^�    A9.	@�|�    @�|�    @�GN �AE�O    A7=�@�~`    @�~`    @�C� �AD�
    A5M�@�@    @�@    @�[� �AF��    A3]�@�     @�     @�+� �AII�    A1mv@�     @�     @ַ �AGJ�    A/}S@��    @��    @�� �AJR,    A-�0@��    @��    @�' �AMnf    A+�@�    @�    @�eP �AQ�$    A)��@�    @�    @�Ƭ �AT�8    A'��@�`    @�`    @�SE �AX[�    A%̰@�@    @�@    @�m� �A[k�    A#ܒ@�     @�     @�wS �AZZ�    A!�u@�     @�     @��M �A\��    A�Y@��    @��    @Ԃ6 �A^�X    A>@��    @��    @�qd �A`�    A$@�    @�    @�GS �A^�N    A,@�    @�    @ߖ� �A]v�    A;�@�`    @�`    @ׁ �A^��    AK�@�@    @�@    @�L �Ag/�    A[�@�     @�     @�,� �A`��    Ak�@�     @�     @�$U �Aa�    A{�@��    @��    @ʣ� �A`r�    A��@��    @��    @ǎ� �A^��    A��@�    @�    @� �A[}�    A
�t@�    @�    @�nV �A\�G    A�g@�`    @�`    @��� �A[�    A�[@�@    @�@    @�ت �AZ��    A�R@�     @�     @��� �A\�    A�I@�     @�     @�%� �A]��    A �C@��    @��    @�1� �Ab�    @�{@��    @��    @�LA �Ab��    @�6u@�    @�    @��s �Af��    @�Vr@�    @�    @��� �Ag�     @�vr@�`    @�`    @�D �Ap��    @�v@�@    @�@    @��� �Ak1|    @�}@�     @�     @�N� �Akt    @�ֈ@��     @��     @�/� �Aj�E    @���@���    @���    @�6� �Aj�9    @��@���    @���    @�*� �Am�x    @�6�@�Š    @�Š    @�� �Anw�    @�V�@�ǀ    @�ǀ    @֣� �Am��    @�v�@��`    @��`    @пi �Ao    @ϗ@��@    @��@    @̵� �Anw�    @˷C@��     @��     @�`l �Aq$�    @��n@��     @��     @�z� �Aq|     @���@���    @���    @��  �ArNx    @��@���    @���    @ϼ� �Ar��    @�8
@�Ԡ    @�Ԡ    @Ҭ8 �AtnZ    @�XH@�ր    @�ր    @��J �Avg.    @�x�@��`    @��`    @��� �AzO�    @���@��@    @��@    @�� �A|��    @��@��     @��     @�`� �A�%�    @��k@��     @��     @�q� �A�DN    @���@���    @���    @��a �A�\�    @�@���    @���    A	�� �A��Q    @�:y@��    @��    A�� �A��    @�Z�@��    @��    A)�0 �A��    @�{F@��`    @��`    A2	� �A��    @���@��@    @��@    AB�&  �A�I    @��)@��     @��     A&�U  �A�z    @�ܢ@��     @��     @��  �A�>    @��!@���    @���    @�K  �Aج    @��@���    @���    @�e�  �A�VO    @||_@��    @��    @�9  �A��    @t�~@��    @��    @�>v  �A��m    @l��@��`    @��`    @�m+  �A��    @e?�@��@    @��@    @���  �A��    @]�#@��     @��     @�?`  �A{�    @U�r@��     @��     @��*  �AS�r    @N�@���    @���    @�C	  �A.�2    @FE5@���    @���    ��  ����      @>��@��    @��    ��  ����      @6�)@��    @��    ��  ����      @/	�@�`    @�`    ��  ����      @'KQ@�@    @�@    ��  ����      @��@�	     @�	     ��  ����      @ά@�     @�     ��  ����      @n@��    @��    ��  ����      @R=@��    @��    ��  ����      @ �@��    @��    ��  ����      ?�
@��    @��    ��  ����      ?�/�@�`    @�`    ��  ����      ?Ҵ	@�@    @�@    ��  ����      ?�83@�     @�     ��  ����      ?��z@�     @�     ��  ����      ?�@�@��    @��    ��  ����      ?��`@��    @��    ��  ����      ?�J @��    @��    ��  ����      ?k�z@�!�    @�!�    ��  ����      ?L�1@�#`    @�#`    ��  ����      ?-�&@�%@    @�%@    ��  ����      ?�Y@�'     @�'     ��  ����      >ߋ�@�)     @�)     ��  ����      >���@�*�    @�*�    ��  ����      >Gm�@�,�    @�,�    ��  ����      =�4�