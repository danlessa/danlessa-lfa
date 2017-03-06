CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.3-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20141013160000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      63 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      241 pixels     center_y      316 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-10-13 20:49:14, using ingest-tsi-12.3-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20141013.160000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:14 2017 UTC       ANDERS_armtime_timestamp      1488830654     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-10-13 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-10-13 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-10-13 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �T; BH  �rdt�M�M@�      @�      B�� �A
/3    B�H@�#�    @�#�    B��� �@�]S    B���@�'�    @�'�    B��� �@�x    B��T@�+@    @�+@    B��� �A��    B���@�/     @�/     B��L �A3��    B�|o@�2�    @�2�    B��2 �AJSc    B�W@�6�    @�6�    B�,+ �A^��    B�0�@�:@    @�:@    B��� �Aj�,    B�	�@�>     @�>     B�� �A�[    B���@�A�    @�A�    B�9� �A�mW    B��E@�E�    @�E�    B� �A���    B���@�I@    @�I@    B�E� �A\&     B�e�@�M     @�M     B� �A]K�    B�:�@�P�    @�P�    B��� �Ar�    B�Y@�T�    @�T�    B��� �A��    B��0@�X@    @�X@    B�d� �A�]T    B��i@�\     @�\     B�P �A���    B��
@�_�    @�_�    B�8� �A���    B�[@�c�    @�c�    B��� �Aµ�    B�,�@�g@    @�g@    B��� �A��    B���@�k     @�k     B�Y� �A�R#    B��@�n�    @�n�    B�� �A��o    B��@�r�    @�r�    B��� �A�e�    B�m�@�v@    @�v@    B��� �A��    B�<�@�z     @�z     B�� �A�˔    B�:@�}�    @�}�    B�
 �A�    B��p@쁀    @쁀    B��� �A��H    B��?@�@    @�@    B�zo �A�H�    B�t�@�     @�     B�s �A���    B�A�@��    @��    B��� �A�W�    B�i@쐀    @쐀    B��� �A���    B���@�@    @�@    B�>7 �A�߉    B���@�     @�     B�� �A��6    B�ry@��    @��    B�' �A�n�    B�=�@쟀    @쟀    B�}� �Aλ�    B��@�@    @�@    B��x �A�Km    B���@�     @�     B�C �A��    B��J@��    @��    B�/� �A�m/    B�h�@쮀    @쮀    B��8 �B˹    B�2�@�@    @�@    B{\ �B�    B��S@�     @�     Bz� �B�g    B���@��    @��    Bw� �B9>    B��$@콀    @콀    Br�6 �B�|    B�X5@��@    @��@    Bn6 �B$�    B�!@��     @��     Bj�= �B%�    B��@���    @���    Bn7� �B�    B��(@�̀    @�̀    Be�u �B"$�    B�zh@��@    @��@    Be�� �B"WU    B�By@��     @��     Bd� �B!�o    B�
[@���    @���    Bb� �B$    B��@�ۀ    @�ۀ    B[$� �B+k?    B���@��@    @��@    BY�F �B-��    B�`�@��     @��     BSQ �B2�n    B�(/@���    @���    BR� �B4(4    B��?@��    @��    BP� �B7W    B��(@��@    @��@    BIc� �B>    B�|�@��     @��     BF� �BAf    B�C�@���    @���    BJ�� �B>?    B�
	@���    @���    BA�� �BG`w    B��e@��@    @��@    BDg� �BD     B���@�     @�     B;� �BLkI    B�\�@��    @��    B9�" �BL�    B�"�@��    @��    B<�  �BJ�}    B��@�@    @�@    B4�� �BS .    B��X@�     @�     B2Ո �BT��    B�s�@��    @��    B.u� �BYÊ    B�9�@��    @��    B+`� �B\�/    B���@�@    @�@    B*vM �B]�    B��L@�     @�     B)) �B^^`    B���@�"�    @�"�    B!�� �Beo�    B�N�@�&�    @�&�    B1Z �Bhl    B��@�*@    @�*@    B� �Bl��    B�خ@�.     @�.     B�+ �Bo!Q    B���@�1�    @�1�    B�: �Bo��    B�bV@�5�    @�5�    B�> �Bp�Q    B�'
@�9@    @�9@    B8� �Bs[�    B��@�=     @�=     B� �Bp�v    B��4@�@�    @�@�    B�z �Bu9    B�t�@�D�    @�D�    B�� �Bu�    B�9@�H@    @�H@    B� �Bt�E    B��b@�L     @�L     B�� �Bw�    B���@�O�    @�O�    BW� �Bw�    B���@�S�    @�S�    BԂ �Bw�6    B�I�@�W@    @�W@    B�t �B|̨    B��@�[     @�[     B} �B{     B���@�^�    @�^�    BOc �B}�    B���@�b�    @�b�    B
� �B}�E    B�Y�@�f@    @�f@    B
� �B}W    B��@�j     @�j     B�A �B�k    B��G@�m�    @�m�    B�� �B�    B���@�q�    @�q�    B`X �B��H    B�h�@�u@    @�u@    B�� �B��b    B�,,@�y     @�y     B �� �B�1    B��@�|�    @�|�    A�8 �B��R    B��,@퀀    @퀀    A�� �B�e�    B�v�@�@    @�@    A��1 �B�v�    B�9�@�     @�     A�� �B��.    B��L@��    @��    A�?& �B��    B���@폀    @폀    A쒡 �B��    B���@�@    @�@    A��9 �B�ZW    B�G @�     @�     A�w �B��    B�
&@��    @��    A�c �B� �    B��@@힀    @힀    A�4 �B�[)    B��P@��@    @��@    A�Qm �B��o    B�SU@��     @��     A�} �B���    B�Q@���    @���    A߮K �B���    B��B@���    @���    A��& �B��[    B��)@��@    @��@    A��� �B��    B�_@��     @��     A��� �B�q�    B�!�@���    @���    A�LR �B�K�    B��@���    @���    A�5 �B�eM    B��i@��@    @��@    A�8" �B���    B�j#@��     @��     Aɽ} �B���    B�,�@���    @���    A�1  �B��    B��|@�ˀ    @�ˀ    A�<V �B�W?    B��@��@    @��@    A��* �B��a    B�t�@��     @��     A�� �B��0    B�7E@���    @���    A�. �B�N�    B���@�ڀ    @�ڀ    A�0! �B�Y�    B��N@��@    @��@    A��� �B�6�    B�~�@��     @��     A��@ �B�>'    B�A9@���    @���    A�G5 �B��8    B��@��    @��    A�I� �B�4�    B��@��@    @��@    A��9 �B�>&    B��e@��     @��     A� �B�Ͼ    B�J�@���    @���    A�.� �B���    B�@���    @���    A��� �B��5    B��T@��@    @��@    A��� �B�(�    B���@�      @�      A��� �B�/s    B�S�@��    @��    A�]* �B��p    B�	@��    @��    A�8L �B��@    B��:@�@    @�@    A�m� �B��    B��d@�     @�     A��~ �B�4Q    B�\�@��    @��    A��L �B��M    B��@��    @��    A��� �B�C    B���@�@    @�@    A�� �B��8    B���@�     @�     A�[� �B��/    B�d�@�!�    @�!�    A�� �B��>    B�&�@�%�    @�%�    A��� �B��*    B���@�)@    @�)@    A��R �B���    B���@�-     @�-     A� � �B�H�    B�l�@�0�    @�0�    A��� �B���    B�.�@�4�    @�4�    A�c� �B��    B���@�8@    @�8@    A��� �B�ma    B���@�<     @�<     A��U �B���    B�t�@�?�    @�?�    A�̴ �B��o    B�6�@�C�    @�C�    A�j� �B��/    B��]@�G@    @�G@    A�&� �B��S    B��6@�K     @�K     A� � �B�=�    B�|
@�N�    @�N�    A�4� �B�dC    B�=�@�R�    @�R�    A�� �B��    B���@�V@    @�V@    A��� �B�x@    B��n@�Z     @�Z     A�1� �B���    B��2@�]�    @�]�    A��� �B�p�    B�D�@�a�    @�a�    A��� �B�u�    B��@�e@    @�e@    A��� �B���    B��e@�i     @�i     A�.} �B��-    B��@�l�    @�l�    A�3� �B��}    B�K�@�p�    @�p�    A��O �B��    B�v@�t@    @�t@    A��� �B���    B��@�x     @�x     A�� �B�޿    B���@�{�    @�{�    A�O� �B���    B�Rf@��    @��    A�� �B��+    B�@�@    @�@    A��� �B��l    B�՟@�     @�     A��g �B���    B��6@��    @��    A�*
 �B���    B�X�@    @    A�k� �B�#�    B�\@�@    @�@    A��{ �B�jB    B���@�     @�     A�w~ �B�t&    B��t@��    @��    A�U$ �B��D    B�^�@    @    A��� �B�A�    B� �@�@    @�@    A�_ �B�?    B��@�     @�     A��B �B��~    B���@��    @��    A�}t �B��    B�d�@    @    A��B  �B��s    B�&t@�@    @�@    A��  �B�;    B���@�     @�     A���  �B�~P    B��\@��    @��    A���  �B���    B�j�@    @    A��  �B��    B�,:@�@    @�@    A��"  �B�J7    B���@��     @��     A�`  �B���    B��@���    @���    A�i  �B��\    B�pr@�ʀ    @�ʀ    A�#l  �B�f�    B�1�@��@    @��@    A�e;  �B���    B�l@��     @��     A��0  �B���    Bi(@���    @���    A��  �B��    B~��@�ـ    @�ـ    A�O  �B�&=    B~n�@��@    @��@    A��*  �B���    B}�?@��     @��     A�NX  �B��%    B}s�@���    @���    A��a  �B���    B|��@��    @��    A��  �B�%�    B|y+@��@    @��@    A�+T  �B�ߝ    B{��@��     @��     A�'  �B���    B{~]@���    @���    A��K  �B��    B{ �@���    @���    A�K�  �B��    Bz�}@��@    @��@    A��d  �B�	�    Bz@��     @��     A�P6  �B���    By��@��    @��    A�|  �B��    By@��    @��    A�{  �B��&    Bx��@�
@    @�
@    A���  �B���    Bx@�     @�     A���  �B�k    Bw�y@��    @��    AsNO  �B�XA    Bw�@��    @��    A�  �B�Yr    Bv�X@�@    @�@    A��  �B���    Bv�@�     @�     A{��  �B�?^    Bu�'@� �    @� �    Ap'<  �B��    Bu�@�$�    @�$�    Agf�  �B���    Bt��@�(@    @�(@    Ap�  �B�@F    Bt#B@�,     @�,     Ac��  �B��p    Bs��@�/�    @�/�    A]�  �B��D    Bs'�@�3�    @�3�    Ad3  �B��    Br�>@�7@    @�7@    A]�+  �B��u    Br,�@�;     @�;     AZ�K  �B���    Bq��@�>�    @�>�    AY�  �B�x    Bq1@�B�    @�B�    A[�  �B��    Bp�]@�F@    @�F@    AW�7  �B�T�    Bp5�@�J     @�J     AR`�  �B��    Bo��@�M�    @�M�    AN��  �B�b�    Bo:@�Q�    @�Q�    AGsT  �B�$    Bn�G@�U@    @�U@    AG��  �B�)�    Bn>z@�Y     @�Y     AM�
  �B�Y�    Bm��@�\�    @�\�    AO�g  �B��     BmB�@�`�    @�`�    AE�  �B��
    Bl�@�d@    @�d@    AT��  �B��    BlG'@�h     @�h     AIx�  �B�F3    Bk�K@�k�    @�k�    AN\�  �B��5    BkKl@�o�    @�o�    A@b�  �B�|    Bj͋@�s@    @�s@    AJ��  �B�*+    BjO�@�w     @�w     AU|  �B��L    Biѿ@�z�    @�z�    AR�  �B�Z    BiS�@�~�    @�~�    AP��  �B�QE    Bh��@�@    @�@    AL�,  �B��V    BhW�@�     @�     AR'$  �B�Ny    Bg�@��    @��    AWn  �B���    Bg\@    @    AY��  �B�T*    Bf�@�@    @�@    A`��  �B��     Bf`!@�     @�     AY��  �B���    Be�$@��    @��    Aby�  �B�dQ    Bed%@    @    A\�  �B�#�    Bd�$@�@    @�@    ATB@  �B�+�    Bdh!@�     @�     A^,�  �B��^    Bc�@��    @��    A\��  �B�>2    Bcl@變    @變    Ac�  �B�&U    Bb�@�@    @�@    Af��  �B��0    Bbo�@�     @�     A\��  �B��    Ba��@��    @��    Ad�+  �B���    Bas�@ﺀ    @ﺀ    Ae��  �B���    B`��@�@    @�@    Af7�  �B�ݐ    B`w�@��     @��     A{�f  �B�8    B_��@���    @���    Asa�  �B�D�    B_{}@�ɀ    @�ɀ    A��~  �B��q    B^�a@��@    @��@    A�,h  �B��    B^B@��     @��     A��6  �B�a�    B^!@���    @���    A��  �B�Oz    B]��@�؀    @�؀    A��;  �B��    B]�@��@    @��@    A���  �B�w�    B\��@��     @��     A�B�  �B��    B\�@���    @���    A���  �B��1    B[�`@��    @��    A�j�  �B�R�    B[3@��@    @��@    A���  �B�1�    BZ�@��     @��     A�p!  �B�d    BZ�@���    @���    A� �  �B��x    BY��@���    @���    A�x�  �B�ZT    BYn@��@    @��@    A��  �B���    BX�8@��     @��     A�L�  �B�Ə    BX @� �    @� �    A�v�  �B���    BW��@��    @��    A��J  �B�(    BW�@��    @��    A���  �B���    BV�O@��    @��    A�؈  �B��    BV@�`    @�`    A�֢  �B��    BU��@�
@    @�
@    A���  �B�I    BU!�@�     @�     A��r  �B��    BT�I@�     @�     A�Z�  �B��    BT%@��    @��    A�r8  �B�.�    BS��@��    @��    A�j>  �B��`    BS(t@��    @��    A��w  �B�gU    BR�*@��    @��    A���  �B���    BR+�@�`    @�`    A���  �B��    BQ��@�@    @�@    A�Β  �B��    BQ/A@�     @�     A���  �B���    BP��@�     @�     A�mB  �B�^�    BP2�@��    @��    A��  �B�oc    BO�K@� �    @� �    A���  �B�h    BO5�@�"�    @�"�    A�Tb  �B�=x    BN��@�$�    @�$�    A�d�  �B��    BN9G@�&`    @�&`    A�{�  �B��R    BM��@�(@    @�(@    A��  �B��z    BM<�@�*     @�*     A�Y�  �B��B    BL�7@�,     @�,     A�@  �B�!�    BL?�@�-�    @�-�    A�.  �B�7X    BK�{@�/�    @�/�    A��L  �B���    BKC@�1�    @�1�    A��{  �B��I    BJĹ@�3�    @�3�    A��8  �B�*V    BJFV@�5`    @�5`    A�-�  �B�)    BI��@�7@    @�7@    A�ת  �B��.    BII�@�9     @�9     A�k>  �B�t�    BH�%@�;     @�;     A�  �B��^    BHL�@�<�    @�<�    A�>�  �B�Z�    BG�T@�>�    @�>�    A�=-  �B���    BGO�@�@�    @�@�    A�
  �B�!�    BF�}@�B�    @�B�    A��  �B��*    BFS@�D`    @�D`    A��   �B��    BEԢ@�F@    @�F@    A��  �B�*    BEV2@�H     @�H     A�WR  �B�O|    BD��@�J     @�J     A��F  �B� <    BDYP@�K�    @�K�    A�?a  �B�T�    BC��@�M�    @�M�    A�B  �B�N�    BC\i@�O�    @�O�    A���  �B��x    BB��@�Q�    @�Q�    A��  �B��    BB_}@�S`    @�S`    A��u  �B�F;    BA�@�U@    @�U@    A�l  �B�*�    BAb�@�W     @�W     A�8v  �B��T    B@�@�Y     @�Y     A��  �B��3    B@e�@�Z�    @�Z�    A�  �B���    B?�@�\�    @�\�    A���  �B�9�    B?h�@�^�    @�^�    A��  �B��i    B>�#@�`�    @�`�    A���  �B�ԁ    B>k�@�b`    @�b`    A��o  �B��E    B=�$@�d@    @�d@    A��/  �B�c;    B=n�@�f     @�f     A��  �B��    B<�"@�h     @�h     A�%�  �B�ڼ    B<q�@�i�    @�i�    A��  �B�;�    B;�@�k�    @�k�    A�S  �B�ի    B;t�@�m�    @�m�    A��1  �B���    B:�@�o�    @�o�    A�+R  �B�d    B:w�@�q`    @�q`    A�I�  �B�0u    B9�@�s@    @�s@    A�q  �B�+�    B9zz@�u     @�u     A��<  �B�9r    B8��@�w     @�w     A�EQ  �B�k:    B8}g@�x�    @�x�    A���  �B�L�    B7��@�z�    @�z�    A�x�  �B�*    B7�P@�|�    @�|�    A�3�  �B���    B7�@�~�    @�~�    A��  �B���    B6�5@��`    @��`    A�E0  �B�/�    B6�@��@    @��@    A�n�  �B�Ӂ    B5�@��     @��     A��6  �B���    B5�@��     @��     A��  �B���    B4��@���    @���    A��'  �B���    B4
d@���    @���    A�D  �B��    B3��@���    @���    A��  �B���    B3>@���    @���    A��  �B��v    B2��@��`    @��`    A���  �B��    B2@�@    @�@    A���  �B��    B1�~@�     @�     A�H�  �B�,�    B1�@�     @�     A�_  �B�"    B0�P@��    @��    A��+  �B��w    B0�@��    @��    A�x  �B��m    B/�@�    @�    A��  �B���    B/�@�    @�    A���  �B��F    B.��@�`    @�`    A�W  �B�Z    B.Q@�@    @�@    A�,�  �B�S�    B-��@�     @�     A�I�  �B��7    B-@�     @�     A�t�  �B�6�    B,�{@��    @��    A���  �B�hv    B, �@��    @��    A���  �B�o�    B+�?@�    @�    A��"  �B�!    B+#�@�    @�    A��Y  �B�o�    B*� @�`    @�`    A�?�  �B�(#    B*&`@�@    @�@    A�O  �B��    B)��@�     @�     A���  �B���    B))@�     @�     A�j(  �B�;�    B(�z@��    @��    A���  �B��    B(+�@��    @��    A��b  �B��F    B'�4@�    @�    A�  �B�ޞ    B'.�@�    @�    A��[  �B���    B&��@�`    @�`    A���  �B��2    B&1E@�@    @�@    A�w�  �B���    B%��@��     @��     A��  �B�Y�    B%3�@��     @��     A�
	  �B�fa    B$�Q@���    @���    A��"  �B��9    B$6�@���    @���    A���  �B��    B#�@�Ǡ    @�Ǡ    A� �  �B�Z]    B#9Y@�ɀ    @�ɀ    A��r  �B�0    B"��@��`    @��`    A��	  �B�j�    B"<@��@    @��@    A���  �B�T;    B!�[@��     @��     A�1�  �B�f,    B!>�@��     @��     A���  �B��k    B �@���    @���    A���  �B��%    B AX@���    @���    A��  �B�    B«@�֠    @�֠    A��M  �B�n�    BC�@�؀    @�؀    A�	5  �B��    B�Q@��`    @��`    A��6  �B�y    BF�@��@    @��@    A��d  �B���    B��@��     @��     A���  �B�(    BIE@��     @��     A��d  �B�2    Bʕ@���    @���    A�n  �B�    BK�@���    @���    A�F�  �B�	    B�5@��    @��    A���  �B���    BN�@��    @��    A�[*  �B�!~    B��@��`    @��`    A��k  �B���    BQ @��@    @��@    A�K�  �B�=�    B�n@��     @��     A�x3  �B���    BS�@��     @��     A���  �B��v    B�@���    @���    A���  �B�2B    BVT@���    @���    A��  �B�6�    Bנ@���    @���    A���  �B��    BX�@���    @���    A��  �B�	    B�7@��`    @��`    A�   �B�%    B[�@��@    @��@    A���  �B���    B��@��     @��     A�7�  �B��    B^@��     @��     A�.�  �B�T�    B�`@���    @���    A���  �B�-    B`�@��    @��    A�=  �B�o    B��@��    @��    A��u  �B���    Bc:@��    @��    A�P  �B���    B�@�`    @�`    A��"  �B�Q'    Be�@�	@    @�	@    A��  �B�H�    B�@�     @�     A|ff  �B�B    BhX@�     @�     A�7�  �B�$c    B�@��    @��    A�L�  �B��    Bj�@��    @��    A�œ  �B��A    B�*@��    @��    A�%�  �B�Λ    Bmp@��    @��    A�m�  �B���    B�@�`    @�`    A�D�  �B���    Bo�@�@    @�@    A��  �B��<    B�?@�     @�     A���  �B���    Br�@�     @�     A�P�  �B�9@    B��@��    @��    A��  �B�k    Bu@��    @��    A�o  �B��    B�N@�!�    @�!�    A��e  �B���    Bw�@�#�    @�#�    A���  �B�ٶ    B
��@�%`    @�%`    A��  �B���    B
z@�'@    @�'@    A���  �B�l    B	�X@�)     @�)     A��  �B���    B	|�@�+     @�+     A�	  �B�8    B��@�,�    @�,�    A��`  �B�g    B@�.�    @�.�    A���  �B���    B ^@�0�    @�0�    A��!  �B���    B��@�2�    @�2�    A��2  �B��X    B�@�4`    @�4`    A��Q  �B��	    B� @�6@    @�6@    A�i  �B��    B`@�8     @�8     A���  �B�ѽ    B��@�:     @�:     A�7.  �B���    B�@�;�    @�;�    A�\�  �B�O�    B�@�=�    @�=�    A��C  �B��T    B
]@�?�    @�?�    A�a�  �B��    B��@�A�    @�A�    A��  �B�j�    B�@�C`    @�C`    A�[�  �B�c    B�@�E@    @�E@    A��X  �B��y    BW@�G     @�G     A�  �B���    B��@�I     @�I     A�<a  �B���    B�@�J�    @�J�    A��  �B���    B �@�L�    @�L�    A�'.  �B�K    B M@�N�    @�N�    A�qK  �B���    A�+@�P�    @�P�    A���  �B��,    A�-�@�R`    @�R`    A�r  �B�    A�0@�T@    @�T@    A�@�  �B��    A�2�@�V     @�V     A��  �B��    A�4�@�X     @�X     A��(  �B�sb    A�7r@�Y�    @�Y�    A��=  �B�RV    A�9�@�[�    @�[�    A���  �B��6    A�<b@�]�    @�]�    A�t|  �B�ym    A�>�@�_�    @�_�    Aә5  �B��t    A�AP@�a`    @�a`    A՛�  �B�6�    A�C�@�c@    @�c@    A�d  �B�S�    A�F=@�e     @�e     A�F  �B�g�    A�H�@�g     @�g     A�8m  �B��i    A�K)@�h�    @�h�    A�N�  �B��     A�M�@�j�    @�j�    A�7�  �B�J�    A�P@�l�    @�l�    A�P�  �B�3    A�R�@�n�    @�n�    A�8`  �B���    A�T�@�p`    @�p`    A���  �B�Z�    A�Wq@�r@    @�r@    A�uR  �B�"�    A�Y�@�t     @�t     A�8�  �B��s    A�\Y@�v     @�v     A�и  �B��8    A�^�@�w�    @�w�    A��  �B�gN    A�a?@�y�    @�y�    A�t  �B�a]    A�c�@�{�    @�{�    A�&  �B��    A�f%@�}�    @�}�    A�MV  �B��n    A�h�@�`    @�`    A��0  �B�k=    A�k
@�@    @�@    A���  �B���    A�m|@�     @�     B �-  �B�
�    A�o�@�     @�     A�9}  �B��D    A�r`@��    @��    A�m�  �B��L    A�t�@��    @��    A��  �B��    A�wC@�    @�    A�k�  �B�S!    A�y�@�    @�    A��  �B�n�    A�|%@�`    @�`    A���  �B�I�    A�~�@�@    @�@    BY�  �B��x    A܁@�     @�     B��  �B���    Aۃx@�     @�     B
�K  �B�,�    Aڅ�@��    @��    B��  �B�    AوY@��    @��    B�>  �B{@)    A؊�@�    @�    BDY  �B}��    A׍9@�    @�    BG  �Bzˢ    A֏�@�`    @�`    B#�  �B|��    AՒ@�@    @�@    B��  �B|)4    AԔ�@�     @�     BG�  �Bxz~    AӖ�@�     @�     Ba�  �B{��    Aҙh@��    @��    B�  �B{Tt    Aћ�@��    @��    Bf.  �B}�     AОG@�    @�    B0Z  �Bz�t    AϠ�@�    @�    BH  �By��    AΣ&@�`    @�`    B�  �BuV�    Aͥ�@�@    @�@    B��  �Bs�    Ą@�     @�     B}�  �Bq�    A˪s@�     @�     BT�  �Bs��    Aʬ�@��    @��    Bdn  �Bs��    AɯR@��    @��    B|
  �Bpv�    Aȱ�@�    @�    B�  �Bn�    AǴ0@�    @�    B@?  �Bo��    Aƶ�@�`    @�`    BO�  �Bm'�    AŹ@�@    @�@    B�  �Bl�z    AĻ}@�     @�     B�(  �Bp�    Aý�@��     @��     B�  �Bm:}    A��[@���    @���    B�`  �Bl�E    A���@���    @���    BVf  �Bm�    A��9@�Ơ    @�Ơ    B�.  �Bn^�    A�Ǩ@�Ȁ    @�Ȁ    B��  �Bk#:    A��@��`    @��`    B[�  �Bj�    A�̆@��@    @��@    B�#  �Bktl    A���@��     @��     B~�  �BmDB    A��e@��     @��     B��  �BmPJ    A���@���    @���    B	�  �Bt#b    A��D@���    @���    B�5  �Bq    A�س@�ՠ    @�ՠ    B��  �Br�_    A��#@�׀    @�׀    B�  �BmE�    A�ݒ@��`    @��`    B"#�  �BgoA    A��@��@    @��@    B"  �Bh)�    A��r@��     @��     B'��  �Ba�L    A���@��     @��     B*͗  �B^�c    A��R@���    @���    B1��  �BW+�    A���@���    @���    B5�"  �BSE    A��2@��    @��    B4`  �BS�=    A��@��    @��    B6��  �BP3L    A��@��`    @��`    B6�\  �BPV�    A��@��@    @��@    B<x�  �BJ�    A���@��     @��     B>��  �BH*y    A��e@��     @��     B?�C  �BF�    A���@���    @���    BC�+  �BA��    A��H@���    @���    BC��  �BA��    A���@��    @��    BG�d  �B>�3    A�+@���    @���    BI�  �B;ʋ    A��@��`    @��`    BI  �B;�w    A�@��@    @��@    BLR8  �B9j    A�	�@��     @��     BI&Q  �B=Y�    A��@��     @��     BG�  �B=n�    A�e@���    @���    BDh8  �BA�    A��@� �    @� �    BE�k  �B?��    A�J@��    @��    BD�  �B@��    A��@��    @��    B@��  �BE.�    A�1@�`    @�`    BA�t  �BC3�    A��@�@    @�@    B=�  �BG:�    A�@�
     @�
     B>[�  �BG�    A��@�     @�     B==�  �BG'    A�" @��    @��    B:��  �BIj�    A�$u@��    @��    B9�  �BL#�    A�&�@��    @��    B<3  �BHr�    A�)^@��    @��    B;��  �BH�    A�+�@�`    @�`    B<�  �BH	f    A�.I@�@    @�@    B=A�  �BG�     A�0�@�     @�     B;�'  �BH/�    A�35@�     @�     B?�  �BD;    A�5�@��    @��    B=�$  �BE��    A�8"@��    @��    B@��  �BB�    A�:�@� �    @� �    BA7  �B@�>    A�=@�"�    @�"�    BD1  �B=�$    A�?�@�$`    @�$`    BC�F  �B=*U    A�A�@�&@    @�&@    BC`;  �B<��    A�Dw@�(     @�(     BDzB  �B;��    A�F�@�*     @�*     BE��  �B9c)    A�Ih@�+�    @�+�    BE�P  �B7{    A�K�@�-�    @�-�    BG �  �B7Xx    A�N[@�/�    @�/�    BJ�  �B2:    A�P�@�1�    @�1�    BI��  �B245    A�SO@�3`    @�3`    BM=3  �B.��    A�U�@�5@    @�5@    BO��  �B(�S    A�XD@�7     @�7     BPJO  �B'-�    A�Z�@�9     @�9     BS�  �B"�    A�];@�:�    @�:�    BP�  �B"/	    A�_�@�<�    @�<�    BOp�  �B"c�    A�b3@�>�    @�>�    BK�@  �B&<%    A�d�@�@�    @�@�    BNsl  �BGx    A~�Z@�B`    @�B`    BNM�  �BnD    A|�U@�D@    @�D@    BM��  �B�    Az�Q@�F     @�F     BOIs  �Bu�    Ax�M@�H     @�H     BK�A  �B��    Av�K@�I�    @�I�    BG�4  �B(*    At�I@�K�    @�K�    BG�(  �B�    Ar�H@�M�    @�M�    BGf�  �B	�    Ap�H@�O�    @�O�    BH��  �Br�    An�I@�Q`    @�Q`    BF��  �B1R    Al�J@�S@    @�S@    BH�M  �B    Ak M@�U     @�U     BH݊  �B�E    AiQ@�W     @�W     BK=�  �B	��    Ag
U@�X�    @�X�    BL�  �B	�    AeZ@�Z�    @�Z�    BQj�  �Bn�    Aca@�\�    @�\�    BQ<�  �B��    Aah@�^�    @�^�    BT	  �A���    A_p@�``    @�``    BVw  �A�M�    A]#y@�b@    @�b@    BW�  �A�g    A[(�@�d     @�d     B]U�  �A�8    AY-�@�f     @�f     B^��  �A淕    AW2�@�g�    @�g�    Ba*�  �A��    AU7�@�i�    @�i�    Bc��  �A�}L    AS<�@�k�    @�k�    Be
�  �Aͧ&    AQA�@�m�    @�m�    Bhi�  �A� w    AOF�@�o`    @�o`    Bi��  �A��\    AMK�@�q@    @�q@    Bj�}  �A��z    AKP�@�s     @�s     Bh��  �A���    AIV@�u     @�u     Bc}  �A�rv    AG[ @�v�    @�v�    B^�  �A�l    AE`6@�x�    @�x�    B\�  �A��B    ACeL@�z�    @�z�    BX��  �A���    AAjd@�|�    @�|�    BU��  �A�k�    A?o}@�~`    @�~`    BS�  �A��z    A=t�@�@    @�@    BOہ  �A���    A;y�@�     @�     BIս  �A�U�    A9~�@�     @�     BHU<  �A�1�    A7��@��    @��    BFF  �A��T    A5�@��    @��    BE��  �A�M�    A3�+@�    @�    BF/�  �A��L    A1�L@�    @�    BBZ�  �A�:f    A/�n@�`    @�`    B>Ր  �A���    A-��@�@    @�@    B>�  �A�<f    A+��@�     @�     B<d  �A�f�    A)��@�     @�     B=��  �A�\m    A'�@��    @��    B<^�  �A�3^    A%�,@��    @��    B7��  �ApDV    A#�V@�    @�    B6t�  �AV@�    A!��@�    @�    B1��  �AE?�    A��@�`    @�`    B)zk  �A/��    A��@�@    @�@    B)��  �A*u�    A�@�     @�     BNC  �A-�    A�;@�     @�     B�+  �A.�x    A�m@��    @��    Bn  �A0�    A۠@��    @��    A���  �A6+�    A��@�    @�    A��  �A/Z�    A�
@�    @�    A؁f  �A5h�    A�B@�`    @�`    A���  �A2��    A�z@�@    @�@    A���  �A/p    A��@�     @�     A�o�  �A-Y*    A	��@�     @�     A�Wo  �A1�4    A -@��    @��    A���  �A2�o    Ak@��    @��    A��Z  �A1�    A
�@�    @�    Aj��  �A0�    A�@�    @�    Ae�a  �A0��    A /@�`    @�`    AM�t  �A.�*    @�4�@�@    @�@    ACs�  �A2(=    @�?s@�     @�     A8�q  �A2��    @�J@��     @��     A.�  �A3 �    @�T�@���    @���    A0��  �A4x    @�_)@���    @���    A-p�  �A3T�    @�i�@�Š    @�Š    A-|�  �A+P�    @�t\@�ǀ    @�ǀ    A0��  �A(O�    @�~�@��`    @��`    A7�  �A#kV    @܉�@��@    @��@    A:�
  �A�2    @ؔ@@��     @��     A;1>  �A��    @Ԟ�@��     @��     AB��  �A�'    @Щ�@���    @���    ACp�  �A�    @̴B@���    @���    AE��  �@�#    @Ⱦ�@�Ԡ    @�Ԡ    AN��  �@�    @�ɩ@�ր    @�ր    AK��  �@�?<    @��a@��`    @��`    ALF�  �@ϫ    @��@��@    @��@    AH�o  �@�
�    @���@��     @��     AG|�  �@��    @���@��     @��     ABD�  �@���    @��e@���    @���    A{ u  �@х�    @�
.@���    @���    A;�l  �@��
    @��@��    @��    A2yX  �@�|    @��@��    @��    A-��  �@��$    @�*�@��`    @��`    A#�z  �@���    @�5w@��@    @��@    A%N�  �@{Ϻ    @�@S@��     @��     A(�"  �@`ƽ    @�K1@��     @��     A#<  �@T��    @�V@���    @���    A�z  �@M�o    @�`�@���    @���    Aހ  �@Cd�    @�k�@��    @��    A�v  �@>4    @�v�@��    @��    Ai  �@1�-    @���@��`    @��`    A 4�  �@ UO    @{o@��@    @��@    @�ڿ  �@Z�    @s/a@��     @��     @�r�  �@��    @kEZ@��     @��     @  �@�"    @c[[@���    @���    @�%�  �@�t    @[qd@���    @���    @��  �@C�    @S�t@��    @��    @�X�  �?�e�    @K��@��    @��    A"j�  �@i��    @C��@�`    @�`    ��  ����      @;��@�@    @�@    ��  ����      @3�@�	     @�	     ��  ����      @+�=@�     @�     ��  ����      @$}@��    @��    ��  ����      @"�@��    @��    ��  ����      @9@��    @��    ��  ����      @Oo@��    @��    ��  ����      @e�@�`    @�`    ��  ����      ?��r@�@    @�@    ��  ����      ?�%V@�     @�     ��  ����      ?�RJ@�     @�     ��  ����      ?�P@��    @��    ��  ����      ?��f@��    @��    ��  ����      ?�َ@��    @��    ��  ����      ?��@�!�    @�!�    ��  ����      ?�4@�#`    @�#`    ��  ����      ?t��@�0�    @�0�    ��  ����      =�u