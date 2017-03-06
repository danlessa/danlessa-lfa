CDF  H   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.4-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20151013110000.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      61 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      236 pixels     center_y      320 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Zcreated by user dsmgr on machine ruby at 2015-10-13 13:49:02, using ingest-tsi-12.4-0.el6      ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20151013.110000.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:45 2017 UTC       ANDERS_armtime_timestamp      1488830685     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2015-10-13 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2015-10-13 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2015-10-13 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �VI�BH  �rdt�M�M@�V     @�V     A% �  �B���    A��<@�Y�    @�Y�    A"��  �B�k}    A���@�]�    @�]�    A'S�  �B��    A���@�a@    @�a@    A��  �B�i    A��K@�e     @�e     A?�  �B���    A���@�h�    @�h�    AXm  �B�9�    A���@�l�    @�l�    A&�  �B��:    A��\@�p@    @�p@    A �{  �B���    A��@�t     @�t     A�  �B�/�    A���@�w�    @�w�    AK3  �B���    A��o@�{�    @�{�    A�  �B���    A�� @�@    @�@    A+�  �B���    A���@�     @�     A�  �B�y    A���@��    @��    A�2  �B�H�    A�5@㊀    @㊀    A!�;  �B���    A�|�@�@    @�@    A$y�  �B�A    A�z�@�     @�     A'a�  �B�	�    A�xL@��    @��    A(�0  �B��+    A�u�@㙀    @㙀    A16  �B�w�    A�s�@�@    @�@    A=~�  �B��"    A�qd@�     @�     AD��  �B�;]    A�o@��    @��    AF�  �B��    A�l�@㨀    @㨀    AWx  �B�=    A�j}@�@    @�@    AU��  �B��C    A�h1@�     @�     AZ�  �B�Q�    A�e�@��    @��    A`�  �B��     A�c�@㷀    @㷀    Aj��  �B��Q    A�aK@�@    @�@    Ag�  �B�(    A�^�@�     @�     Aht  �B�"J    A�\�@���    @���    A`Lf  �B�.    A�Zg@�ƀ    @�ƀ    Ao��  �B�\�    A�X@��@    @��@    AfO  �B��[    A�U�@��     @��     Al�2  �B��o    A�S�@���    @���    Acy9  �B�1    A�Q8@�Հ    @�Հ    Ap�V  �B�f�    A�N�@��@    @��@    AP��  �B�x�    A�L�@��     @��     AG[�  �B��P    A�JV@���    @���    AH��  �B�i�    A�H
@��    @��    A+,�  �B�*f    A�E�@��@    @��@    A��  �B�3S    A�Cs@��     @��     A�	  �B�6-    A�A(@���    @���    A�2  �B�$5    A�>�@��    @��    @��J  �B�Q�    A�<�@��@    @��@    @��  �B�h�    A�:F@��     @��     @��>  �B��    A�7�@���    @���    @�*�  �B�R}    A�5�@��    @��    @��H  �B��~    A�3e@�@    @�@    @�5_  �B�;*    A�1@�
     @�
     @�\  �B���    A�.�@��    @��    A	�  �B��#    A�,�@��    @��    A7�  �B��    A�*8@�@    @�@    Av�  �B���    A�'�@�     @�     A��  �B��i    A�%�@��    @��    A��  �B��<    A�#V@� �    @� �    A��  �B���    A�!
@�$@    @�$@    A
~  �B�cj    A��@�(     @�(     A7�  �B���    A�s@�+�    @�+�    A�  �B��    A�(@�/�    @�/�    AM�  �B���    A��@�3@    @�3@    Ad�  �B���    A��@�7     @�7     @�в  �B�Љ    A�E@�:�    @�:�    @�f�  �B�M�    A��@�>�    @�>�    @�]+  �B��R    A��@�B@    @�B@    @Բ�  �B�n    A�a@�F     @�F     @̇�  �B��    A�
@�I�    @�I�    @�>�  �B��w    A��@�M�    @�M�    @�5>  �B��b    A�|@�Q@    @�Q@    @��%  �B�	�    A�/@�U     @�U     @�ɪ  �B�R    A� �@�X�    @�X�    @��  �B��n    A���@�\�    @�\�    @��u  �B�?b    A��I@�`@    @�`@    @��  �B�]�    A���@�d     @�d     @��q  �B��6    A���@�g�    @�g�    @��  �B���    A��b@�k�    @�k�    @��6  �B�:W    A��@�o@    @�o@    @��  �B��i    A���@�s     @�s     @�۴  �B��/    A��y@�v�    @�v�    @�Y�  �B�Ex    A��+@�z�    @�z�    @�'  �B���    A���@�~@    @�~@    @�`  �B�;�    A��@�     @�     @�b  �B�j    A��?@��    @��    @��  �B�s    A���@䉀    @䉀    @�{8  �B���    A��@�@    @�@    @�x�  �B�S�    A��R@�     @�     @�}�  �B���    A��@��    @��    @�%I  �B�A�    A�ٳ@䘀    @䘀    @��'  �B���    A��c@�@    @�@    @�[z  �B�1�    A��@�     @�     @ǒ$  �B�C�    A���@��    @��    @�4}  �B�{�    A��r@䧀    @䧀    @��  �B���    A��!@�@    @�@    @�g�  �B�s#    A���@�     @�     @��L  �B���    A��~@��    @��    @�*�  �B��b    A��-@䶀    @䶀    @��  �B��    A���@�@    @�@    @�b�  �B�&S    A�@�     @�     @��  �B�~N    A��6@���    @���    AUw  �B�S�    A��@�ŀ    @�ŀ    Atr  �B�    A�@��@    @��@    A�h  �B���    A�<@��     @��     AR\  �B�2�    A��@���    @���    @��]  �B�e    A���@�Ԁ    @�Ԁ    AC  �B�v�    A��?@��@    @��@    @�)d  �B��@    A���@��     @��     @�s�  �B�g    A���@���    @���    @��  �B�n�    A��?@��    @��    @�ڊ  �B�3}    A���@��@    @��@    @�Z:  �B�'�    A���@��     @��     @��   �B�b�    A��<@���    @���    @�z�  �B��    A���@��    @��    @��}  �B�~�    A���@��@    @��@    @�K5  �B�$�    A��5@��     @��     @�b  �B��    A���@���    @���    @�\  �B�$0    B LB@��    @��    @�ۺ  �B�#�    B �@�@    @�@    @�·  �B�O>    BI�@�	     @�	     @��  �B��    BȻ@��    @��    @}��  �B�Ղ    BG�@��    @��    @��  �B���    B�`@�@    @�@    @��  �B���    BE2@�     @�     @�:  �B�]>    B�@��    @��    @�<q  �B�m�    BB�@��    @��    @x��  �B��i    B��@�#@    @�#@    @_�  �B��@    B@y@�'     @�'     @G"  �B��k    B�J@�*�    @�*�    @30Z  �B�$�    B>@�.�    @�.�    @1B\  �B�5�    B��@�2@    @�2@    @!k  �B´	    B;�@�6     @�6     @>�  �B���    B��@�9�    @�9�    @B�  �Bæ    B9[@�=�    @�=�    @��  �B� �    B�*@�A@    @�A@    @ �  �Bº    B	6�@�E     @�E     @=  �BåJ    B	��@�H�    @�H�    @	��  �B�p'    B
4�@�L�    @�L�    ?���  �B��^    B
�e@�P@    @�P@    @e�  �B��    B23@�T     @�T     @�?  �BÄ�    B�@�W�    @�W�    @��  �BÌ�    B/�@�[�    @�[�    ?�9�  �B�    B��@�_@    @�_@    ?���  �B�8    B-h@�c     @�c     @ ��  �Bü�    B�4@�f�    @�f�    @��  �B�{�    B+ @�j�    @�j�    @��  �B�JO    B��@�n@    @�n@    @u  �Bé�    B(�@�r     @�r     @��  �BÖ�    B�c@�u�    @�u�    @��  �B��;    B&-@�y�    @�y�    @?O�  �B��     B��@�}@    @�}@    @"��  �B¤�    B#�@�     @�     @C�i  �B��    B��@��    @��    @M�  �B�H�    B!U@刀    @刀    @l �  �B�X�    B�@�@    @�@    @[��  �B��M    B�@�     @�     @f�k  �B���    B��@��    @��    @���  �B���    Bv@嗀    @嗀    @x(  �B���    B�>@�@    @�@    @�05  �B���    B@�     @�     @���  �B�2�    B��@��    @��    @�%j  �B���    B�@妀    @妀    @�}�  �B�N�    B�X@�@    @�@    @�U�  �B�ނ    B@�     @�     @�(_  �B�p�    B��@��    @��    @��+  �B��    B�@嵀    @嵀    @��I  �B���    B�k@�@    @�@    @�j�  �B��(    B/@�     @�     @�K�  �B�    B��@���    @���    @�r  �B�~�    B�@�Ā    @�Ā    @���  �B���    B�w@��@    @��@    @��  �B���    B9@��     @��     @���  �B��    B��@���    @���    @�\�  �B��N    B�@�Ӏ    @�Ӏ    @��?  �B��    B�}@��@    @��@    @�B  �B�m�    B=@��     @��     @���  �B��2    B��@���    @���    @��  �B�\d    B�@��    @��    @켉  �B���    B�z@��@    @��@    @���  �B��x    B9@��     @��     @ؠ�  �B�     B�@���    @���    @�U�  �B�=�    B��@��    @��    @ݔ�  �B���    B }p@��@    @��@    @�'�  �B�'�    B �-@��     @��     @��  �B�X�    B!z�@���    @���    @��.  �B���    B!��@� �    @� �    A B$  �B��    B"x^@�@    @�@    A�  �B�,�    B"�@�     @�     AI1  �B�y�    B#u�@��    @��    A(7  �B��,    B#�@��    @��    A�  �B���    B$sC@�@    @�@    A�;  �B�w�    B$��@�     @�     Az  �B�6    B%p�@��    @��    A"��  �B�`�    B%�i@��    @��    A*��  �B�c�    B&n@�"@    @�"@    A)�  �B��#    B&��@�&     @�&     A-��  �B���    B'k�@�)�    @�)�    A#=  �B�Q�    B'�>@�-�    @�-�    A6�  �B���    B(h�@�1@    @�1@    A6�   �B��F    B(�@�5     @�5     AN��  �B���    B)fW@�8�    @�8�    AJ�  �B�|�    B)�	@�<�    @�<�    AQB   �B���    B*c�@�@@    @�@@    ASc�  �B�W�    B*�k@�D     @�D     AUճ  �B��    B+a@�G�    @�G�    AP`?  �B�Ì    B+��@�K�    @�K�    A]O�  �B��    B,^x@�O@    @�O@    AY�?  �B��    B,�&@�S     @�S     AT�  �B�H�    B-[�@�V�    @�V�    AN�  �B��    B-ڀ@�Z�    @�Z�    AN�   �B���    B.Y,@�^@    @�^@    AF�  �B��    B.��@�b     @�b     A?�  �B���    B/V�@�e�    @�e�    A2�*  �B�i7    B/�+@�i�    @�i�    A9  �B���    B0S�@�m@    @�m@    A0�C  �B���    B0�}@�q     @�q     A(E  �B��O    B1Q$@�t�    @�t�    A�  �B��    B1��@�x�    @�x�    Akp  �B��    B2Nq@�|@    @�|@    A��  �B��    B2�@�     @�     A�  �B�     B3K�@��    @��    A)D�  �B���    B3�^@懀    @懀    A)�*  �B���    B4I@�@    @�@    A)f�  �B���    B4ǣ@�     @�     A#_'  �B�cS    B5FE@��    @��    A"_�  �B��Y    B5��@斀    @斀    A+y�  �B�b%    B6C�@�@    @�@    A4��  �B�=t    B6�$@�     @�     A7�  �B��    B7@�@��    @��    AK,m  �B�f�    B7�_@楀    @楀    Ac�#  �B�Dj    B8=�@�@    @�@    Avp�  �B���    B8��@�     @�     Au��  �B�	a    B9;1@��    @��    A��  �B��    B9��@洀    @洀    A�nf  �B�d�    B:8d@�@    @�@    A��  �B���    B:��@�     @�     A��  �B�Y)    B;5�@��    @��    A��^ �B�Ҡ    B;�(@�À    @�À    A�9s  �B�x>    B<2�@��@    @��@    A�1  �B�{�    B<�R@��     @��     A�Q�  �B�h�    B=/�@���    @���    A��\ �B���    B=�w@�Ҁ    @�Ҁ    A�M�  �B�+S    B>-@��@    @��@    A�E  �B�e3    B>��@��     @��     A�p  �B�+�    B?*(@���    @���    A�]�  �B���    B?��@��    @��    AР�  �B��o    B@'C@��@    @��@    A�1V  �B�*}    B@��@��     @��     A���  �B���    BA$[@���    @���    A�
�  �B���    BA��@���    @���    A���  �B���    BB!n@��@    @��@    A�  �B�j�    BB��@��     @��     A��  �B�e�    BC}@���    @���    A�7U  �B�    BC�@���    @���    A�%�  �B��|    BD�@�@    @�@    A��  �B��    BD�@�     @�     Aۇ�  �B��    BE�@�
�    @�
�    A�Z�  �B��x    BE�@��    @��    A�uj  �B���    BF�@�@    @�@    A��j  �B���    BF�@�     @�     A�c�  �B�%�    BG�@��    @��    A��U  �B��B    BG�	@��    @��    A�8\  �B�'�    BH�@�!@    @�!@    A��]  �B�9    BH��@�%     @�%     Bi �Bvw�    BIw@�(�    @�(�    B*�h �Bax    BI��@�,�    @�,�    B+�) �Ba�    BJ	e@�0@    @�0@    B3�<  �BY��    BJ��@�4     @�4     B=�V  �BO�    BKN@�7�    @�7�    BJ\- �B=�    BK��@�;�    @�;�    BI]r �B;��    BL2@�?@    @�?@    BC[ �B@Ԃ    BL��@�C     @�C     B/� �BN��    BM @�F�    @�F�    B�� �Bb@    BM~~@�J�    @�J�    B �  �Bb:�    BM��@�N@    @�N@    B$�x �B^�F    BN{U@�R     @�R     B#�� �Bf�    BN��@�U�    @�U�    B�� �Bk��    BOx&@�Y�    @�Y�    B�4 �Bw��    BO��@�]@    @�]@    B� �B{h�    BPt�@�a     @�a     B(� �Bwq^    BP�T@�d�    @�d�    B��  �B��    BQq�@�h�    @�h�    B��  �B~�    BQ�@�l@    @�l@    B	C�  �B��    BRnu@�p     @�p     B �  �B�J�    BR��@�s�    @�s�    A�t!  �B�Z}    BSk.@�w�    @�w�    A��� �B��    BS�@�{@    @�{@    A��� �B�^�    BTg�@�     @�     A��� �B�cc    BT�8@��    @��    A��� �B��8    BUd�@熀    @熀    A��L  �B�ϐ    BU��@�@    @�@    A�$  �B�˞    BVa4@�     @�     A�W0  �B�0�    BV߄@��    @��    A�}�  �B�fp    BW]�@畀    @畀    A�6�  �B�p�    BW� @�@    @�@    A���  �B�K�    BXZl@�     @�     A`�i  �B���    BXص@��    @��    AQ��  �B���    BYV�@礀    @礀    AI�?  �B���    BY�C@�@    @�@    A@%Q  �B�Ǭ    BZS�@�     @�     A:V#  �B�t�    BZ��@��    @��    A7�  �B��    B[P@糀    @糀    A@�  �B�Έ    B[�I@�@    @�@    A87�  �B���    B\L�@�     @�     A6`�  �B���    B\��@��    @��    A;�  �B�"�    B]H�@�    @�    A��  �B��    B]�1@��@    @��@    @��  �B��b    B^Ef@��     @��     @�f  �B�    B^Ù@���    @���    A͉  �B��&    B_A�@�р    @�р    A,ؒ  �B�)\    B_��@��@    @��@    Asvq  �B�[�    B`>&@��     @��     A��`  �B��:    B`�Q@���    @���    A�5C  �B��    Ba:z@���    @���    A���  �B�H2    Ba��@��@    @��@    B�
  �B�D�    Bb6�@��     @��     B��  �B}��    Bb��@���    @���    B{�  �Bt��    Bc3@��    @��    B',  �Bh^�    Bc�$@��@    @��@    B-�w  �Ba��    Bd/@@��     @��     B&��  �Bf��    Bd�Y@���    @���    B�`  �BxZW    Be+p@���    @���    B�  �B}��    Be��@�@    @�@    B$T  �B~�    Bf'�@�     @�     B̏  �Bg�    Bf��@�	�    @�	�    B�  �B��    Bg#�@��    @��    B:B  �B�d    Bg��@�@    @�@    B~�  �B{ӫ    Bh�@�     @�     B#  �Bp��    Bh��@��    @��    B!#�  �Bmö    Bi�@��    @��    B$�A  �Bi��    Bi��@� @    @� @    B&H�  �Bf�S    Bj�@�$     @�$     B"�{  �Bd%�    Bj��@�'�    @�'�    B!� �Bbϕ    Bk�@�+�    @�+�    B-�R �BZ�-    Bk��@�/@    @�/@    B>�y �BI�}    Bl�@�3     @�3     BX<P �B2�    Bl��@�6�    @�6�    Bmm� �B��    Bm�@�:�    @�:�    B~!� �BY�    Bm�q@�>@    @�>@    B�� �B/�    BnY@�B     @�B     B|�1 �B(�    Bn�=@�E�    @�E�    Bz4) �A��k    Bo@�I�    @�I�    B{3f �A�5�    Bo��@�M@    @�M@    Bon� �A���    Bo��@�Q     @�Q     BpTQ �A���    Bp|�@�T�    @�T�    Br$� �A�@    Bp��@�X�    @�X�    Bu$ �A��    BqxY@�\@    @�\@    B�� �Aԟ�    Bq�(@�`     @�`     B�ݙ �A��    Brs�@�c�    @�c�    B�px �A�~�    Br�@�g�    @�g�    B�n� �A�/M    Bso�@�k@    @�k@    B�Z� �A�    Bs�C@�o     @�o     B�� �A�͠    Btk@�r�    @�r�    B��J �A�p�    Bt�@�v�    @�v�    B��� �A���    Bufs@�z@    @�z@    B�p� �A�32    Bu�&@�~     @�~     B� �A�Ā    Bva�@��    @��    B��/ �A� �    Bv߃@腀    @腀    B�9�  �Aђ�    Bw]+@�@    @�@    B� � �A�y�    Bw��@�     @�     B�i�  �A� /    BxXq@��    @��    B�ý  �A�;8    Bx�@蔀    @蔀    B�nL  �A���    ByS�@�@    @�@    B�.�  �A��    By�=@�     @�     B��S  �A�{�    BzN�@��    @��    B�zS  �A���    Bz�[@裀    @裀    B�w  �A���    B{I�@�@    @�@    B��o  �A��    B{�i@�     @�     B��  �A��$    B|D�@��    @��    B�$�  �A�-    B|�f@貀    @貀    B�?]  �A��f    B}?�@�@    @�@    B���  �A�k�    B}�Q@�     @�     B�+  �A�    B~:�@��    @��    B�ʜ  �A��    B~�*@���    @���    B�)  �A��    B5�@��@    @��@    B��  �A�d    B��@��     @��     B�(�  �A���    B�&@���    @���    B���  �A��    B�V�@�Ѐ    @�Ѐ    B�Ӡ  �A�L�    B��{@��@    @��@    B��_  �A�ޟ    B��"@��     @��     B�<.  �A�(    B��@���    @���    B���  �A�YZ    B�Qg@�߀    @�߀    B�vz  �A��    B��@��@    @��@    B�\�  �A��!    B�Σ@��     @��     B��5  �A��0    B�<@���    @���    B�ݍ  �A�p�    B�K�@��    @��    B��!  �A�|�    B��g@��@    @��@    B��R  �A�9�    B���@��     @��     B�!�  �A�ʚ    B��@���    @���    B�|  �A�I    B�F@���    @���    B���  �A���    B���@�@    @�@    B�>  �A�    B��!@�     @�     B��H  �A�L    B��@��    @��    B�3�  �A�z�    B�@#@��    @��    B��=  �Aȹ�    B�~�@�@    @�@    B�m�  �A�EE    B��@�     @�     B�P�  �A��Q    B���@��    @��    B��1  �A���    B�:@��    @��    B���  �A���    B�xq@�@    @�@    B��x  �A�    B���@�#     @�#     B���  �A���    B��F@�&�    @�&�    B���  �A���    B�3�@�*�    @�*�    B��o  �A�q    B�r@�.@    @�.@    B���  �A���    B��k@�2     @�2     B��[  �A��    B���@�5�    @�5�    B���  �A�N�    B�-@�9�    @�9�    B�G�  �A��h    B�kn@�=@    @�=@    B�a<  �A���    B���@�A     @�A     B��T  �A��n    B��@�D�    @�D�    B�2  �A�      B�&O@�H�    @�H�    B��!  �A��    B�d�@�L@    @�L@    B���  �A�x7    B���@�P     @�P     B���  �A���    B��@�S�    @�S�    B���  �A��U    B�B@�W�    @�W�    B�b �A��    B�]t@�[@    @�[@    B�^ �A��-    B���@�_     @�_     B�b5 �A�`j    B���@�b�    @�b�    B�u� �A��    B��@�f�    @�f�    B�9N �A��u    B�V@�j@    @�j@    B�f� �Aķ    B��+@�n     @�n     B��� �A�Hl    B��B@�q�    @�q�    B��� �A�p�    B�S@�u�    @�u�    B�b2 �Aԧ�    B�N`@�y@    @�y@    B�"h �A��    B��g@�}     @�}     B�I� �A�9m    B��i@��    @��    B�Zo  �A��.    B�g@鄀    @鄀    B���  �A�|8    B�F^@�@    @�@    B���  �A땫    B��P@�     @�     Br #  �BW    B��=@��    @��    Bad�  �B-�    B� $@铀    @铀    BZ}  �B4�    B�>@�@    @�@    Bb�5  �B,��    B�{�@�     @�     BO@  �B@<�    B���@��    @��    BWt  �B7�    B���@颀    @颀    B\i�  �B2�Z    B�5N@�@    @�@    Bb� �B,�    B�s@�     @�     BpHF �B�    B���@��    @��    B�ʨ �B�_    B��@鱀    @鱀    B�_� �B��    B�,2@�@    @�@    B��� �B��    B�i�@�     @�     B��� �A��    B��{@��    @��    B��d �Bz    B��@���    @���    B� �B2)    B�"�@��@    @��@    Bi{�  �B%�    B�`2@��     @��     B�i>  �B��    B���@���    @���    By8�  �B1,    B��2@�π    @�π    B}�  �B��    B��@��@    @��@    Bu��  �B�Q    B�V@��     @��     Br}�  �B�P    B��u@���    @���    BS�  �B;h�    B���@�ހ    @�ހ    BSր  �B;pl    B�#@��@    @��@    BY�f  �B5_    B�Km@��     @��     Bd�  �B*yM    B���@���    @���    B[	F  �B4W{    B���@��    @��    BL�  �BBj     B�@��@    @��@    BQ=�  �B>I    B�@9@��     @��     Bnh� �B!m    B�}T@���    @���    Bji� �B%t    B��e@���    @���    Bc�] �B+u�    B��l@� @    @� @    B?��  �BOL�    B�4i@�     @�     B='�  �BR�    B�qZ@��    @��    B]��  �B1��    B��A@��    @��    Bf��  �B(r�    B��@�@    @�@    Br~?  �B��    B�'�@�     @�     Bp��  �B+�    B�d�@��    @��    Bh\  �B'<�    B��j@��    @��    BqR� �B k    B��@�@    @�@    Ba��  �B-�e    B��@�"     @�"     BU��  �B9��    B�WI@�%�    @�%�    BK�  �BCa�    B���@�)�    @�)�    BX-Z �B6�    B��G@�-@    @�-@    BI�g �BEaP    B��@�1     @�1     BCQ�  �BLI    B�I@�4�    @�4�    B:�  �BT�%    B��\@�8�    @�8�    B6�K  �BX��    B���@�<@    @�<@    B5� �BZ>�    B���@�@     @�@     B%� �Bi{z    B�9�@�C�    @�C�    Bu  �B���    B�u�@�G�    @�G�    Bk�  �B���    B���@�K@    @�K@    B��  �B�S�    B���@�O     @�O     B�  �B���    B�)�@�R�    @�R�    A��a  �B��k    B�e�@�V�    @�V�    A�v  �B�*�    B��K@�Z@    @�Z@    A�f[  �B�
)    B���@�^     @�^     A�&�  �B���    B��@�a�    @�a�    A�;�  �B��    B�T@�e�    @�e�    A�  �B�|6    B��p@�i@    @�i@    A�9I  �B���    B���@�m     @�m     A��  �B�Jk    B�@�p�    @�p�    A���  �B�v    B�A2@�t�    @�t�    A���  �B��    B�|F@�x@    @�x@    A�mO  �B�ہ    B��C@�|     @�|     A� �  �B��h    B��(@��    @��    B,�  �BcK�    B�,�@ꃀ    @ꃀ    B5�� �BY˰    B�g�@�@    @�@    BFr] �BI�    B��?@�     @�     B,A�  �Bc$    B�ܽ@��    @��    B6N%  �BY3x    B�@ꒀ    @ꒀ    B:�]  �BT�v    B�Qd@�@    @�@    B?JL  �BP�    B���@�     @�     B<ŭ  �BR��    B�Ŗ@��    @��    B4��  �BZ��    B���@ꡀ    @ꡀ    B2��  �B\˖    B�9J@�@    @�@    B+^"  �Bc�4    B�r�@�     @�     B$0]  �Bk`    B��y@��    @��    B��  �Brq
    B���@가    @가    B�R  �Bs�    B�@�@    @�@    B��  �Bva�    B�X3@�     @�     B�9  �Bt��    B��$@��    @��    B"R  �Bt`�    B���@꿀    @꿀    Bz�  �B}�    B��@��@    @��@    Bi  �B�}j    B�;@��     @��     B�O  �B��)    B�sI@���    @���    A�vF  �B�u&    B��c@�΀    @�΀    A���  �B�7~    B��N@��@    @��@    Aޖ�  �B��    B�@��     @��     B+�  �BRw    B�R�@���    @���    A��P  �B�1�    B���@�݀    @�݀    B	��  �B��n    B���@��@    @��@    B��  �B�co    B���@��     @��     B,G  �Bc0    B�.�@���    @���    B5\  �BZ4    B�d�@��    @��    B@�W  �BN�D    B��@��@    @��@    BAZ  �BM�    B��@��     @��     BI�;  �BEQ    B��@���    @���    B?μ  �BO]    B�<@���    @���    BG�  �BH&     B�q"@��@    @��@    BM�r �BA�/    B���@�     @�     BN�a  �B@2�    B��k@��    @��    BJ�  �BD��    B��@�
�    @�
�    BK�/  �BC��    B�Bm@�@    @�@    BF�G  �BH    B�u�@�     @�     B8a  �BV��    B��@��    @��    B7�q  �BW>    B���@��    @��    B
�  �B��    B�1@�@    @�@    B @�  �Bn�n    B�@*@�!     @�!     B�{  �Bq�F    B�q�@�$�    @�$�    B#&  �Bz+    B���@�(�    @�(�    A��  �B�sI    B�ӂ@�,@    @�,@    A��  �B�jx    B��@�0     @�0     B  P  �B���    B�3h@�3�    @�3�    A�v  �B���    B�b�@�7�    @�7�    A�\  �B��~    B��B@�;@    @�;@    A���  �B��    B��]@�?     @�?     A��j  �B�@u    B���@�B�    @�B�    A���  �B�e�    B��@�F�    @�F�    A��M  �B�J    B�F@�J@    @�J@    A���  �B�    B�q�@�N     @�N     A���  �B���    B���@�Q�    @�Q�    Aۮ�  �B��     B��@�U�    @�U�    B��  �By��    B���@�Y@    @�Y@    BK�  �Bu��    B�M@�]     @�]     B'V�  �Bg�l    B�AD@�`�    @�`�    B-�I  �BaA�    B�hg@�d�    @�d�    B,��  �Bbf�    B���@�h@    @�h@    B��  �B�q�    B��@�l     @�l     BLw  �Br��    B�؁@�o�    @�o�    B�  �B��P    B���@�s�    @�s�    A��a  �B��    B�u@�w@    @�w@    A��  �B��-    B�?�@�{     @�{     A��6  �B��    B�`>@�~�    @�~�    A�a�  �B��b    B�y@낀    @낀    A��  �B���    B���@�@    @�@    A��3  �B��g    B��j@�     @�     A�
<  �B�6�    B��@��    @��    A�&m  �B�     B��f@둀    @둀    A�V  �B��    B�	o@�@    @�@    A�  �B���    B�!@�     @�     Aw��  �B���    B�7d@��    @��    Ae��  �B��i    B�L<@렀    @렀    ATz  �B��V    B�_�@�@    @�@    Abcx  �B�@�    B�q{@�     @�     Al9�  �B� R    B���@��    @��    Ary  �B�"K    B���@므    @므    A��  �B�zm    B���@�@    @�@    A��y  �B��g    B��F@�     @�     B+  �B�u�    B��*@��    @��    B��  �Bx	y    B��b@뾀    @뾀    A�0  �B���    B���@��@    @��@    A�e  �B��-    B�ƿ@��     @��     B	M  �Bpb    B���@���    @���    B�a  �B��@    B��E@�̀    @�̀    A�M/  �B�p�    B���@��@    @��@    A�Z�  �B�0A    B���@��     @��     B>s  �B���    B��(@���    @���    B-b~  �B_*+    B���@�܀    @�܀    B(�x  �Bdi�    B���@��@    @��@    B$%(  �Bh�g    B���@��     @��     B)Tl  �Bd_V    B��1@���    @���    B#��  �BjD    B��@��    @��    B'z   �BfG�    B��J@��@    @��@    B ��  �BmY    B�{�@��     @��     A�tw  �B��    B�k
@���    @���    B�C  �B�~    B�X�@���    @���    B�  �B~�T    B�D�@��@    @��@    BX�  �B�    B�/M@�     @�     A�q  �B�+N    B��@��    @��    A�i�  �B��    B� Q@�	�    @�	�    A�  �B�,    B���@�@    @�@    A���  �B�M    B���@�     @�     A��8  �B���    B���@��    @��    Au�a  �B���    B���@��    @��    A�p�  �B�X    B�t	@�@    @�@    A�v  �B��5    B�Tc@�      @�      A� �  �B�5�    B�3�@�#�    @�#�    AF��  �B��     B��@�'�    @�'�    Ag��  �B�yX    B���@�+@    @�+@    A:�J  �B��F    B��@�/     @�/     A/�}  �B�xU    B��S@�2�    @�2�    A5��  �B���    B���@�6�    @�6�    AJ}�  �B��q    B�Z@�:@    @�:@    A_{  �B�K    B�2�@�>     @�>     Ae��  �B�6�    B�
L@�A�    @�A�    A��S  �B���    B��@@�E�    @�E�    Av"�  �B�\|    B��r@�I@    @�I@    A�kp  �B��    B���@�M     @�M     A���  �B���    B�a�@�P�    @�P�    A�@   �B��j    B�5�@�T�    @�T�    A��  �B�x�    B�	C@�X@    @�X@    A�LT  �B�5    B��@�\     @�\     A���  �B��    B��]@�_�    @�_�    Aגb  �B��    B��@�c�    @�c�    A�J�  �B�{�    B�Q1@�g@    @�g@    Bd  �B�=�    B�!�@�k     @�k     B"~  �B�B>    B���@�n�    @�n�    A�Y�  �B��    B���@�r�    @�r�    A��L  �B��    B���@�v@    @�v@    A�  �B���    B�_p@�z     @�z     A�ͱ  �B��L    B�-�@�}�    @�}�    A��  �B���    B���@쁀    @쁀    A��  �B��1    B��@�@    @�@    A���  �B�I    B��.@�     @�     A�i  �B�z    B�b�@��    @��    A��:  �B��    B�/K@쐀    @쐀    A�k  �B���    B��U@�@    @�@    A�u  �B���    B��@�     @�     A�(�  �B�4�    B��r@��    @��    AĿ�  �B�9�    B�]�@쟀    @쟀    A�Q�  �B���    B�(X@�@    @�@    A�ܦ  �B�#=    B���@�     @�     A��  �B�w6    B��@��    @��    A��  �B��K    B��@쮀    @쮀    A��  �B���    B�P�@�@    @�@    A���  �B��    B�X@�     @�     A���  �B��    B��@��    @��    A��  �B��&    B���@콀    @콀    A�|w  �B�    B�uu@��@    @��@    A�n�  �B�1�    B�>@��     @��     A��  �B���    B�x@���    @���    A�E  �B�`�    B�ά@�̀    @�̀    A�P6  �B�WM    B���@��@    @��@    Aպ#  �B�h�    B�^�@��     @��     A�R	  �B��"    B�&.@���    @���    AԪW  �B���    B���@�ۀ    @�ۀ    A�M  �B�-    B���@��@    @��@    A�4�  �B�.�    B�|$@��     @��     A�~  �B��.    B�C$@���    @���    A֯�  �B�g�    B�	�@��    @��    A˫�  �B�o|    B�г@��@    @��@    A��  �B�V    B��C@��     @��     Aʑo  �B�V�    B�]�@���    @���    AɌ�  �B�x�    B�#�@���    @���    A�9a  �B��    B��'@��@    @��@    A�d�  �B�X�    B��2@�     @�     B�; �B��G    B�v@��    @��    B  �B��.    B�;�@��    @��    A˭1  �B�H�    B��@�@    @�@    A��  �B�K     B��0@�     @�     A�V=  �B�<�    B���@��    @��    A���  �B�k=    B�R@��    @��    A���  �B���    B�L@�@    @�@    A�p�  �B��"    B��w@�     @�     A�M  �B�=�    B���@�"�    @�"�    A�[�  �B���    B�f�@�&�    @�&�    A�G2  �B�u9    B�+i@�*@    @�*@    Aƽc  �B�V+    B��7@�.     @�.     A·�  �B�Q�    B���@�1�    @�1�    A��  �B���    B�y�@�5�    @�5�    A��  �B�&$    B�>@�9@    @�9@    A��R  �B���    B��@�=     @�=     B��  �B{=�    B���@�@�    @�@�    B��  �Bt'>    B��S@�D�    @�D�    B �k  �Bl%0    B�O�@�H@    @�H@    B �  �Bl�.    B��@�L     @�L     B+U�  �Bb�S    B���@�O�    @�O�    B,�  �BaVC    B���@�S�    @�S�    B/�&  �B^U�    B�_�@�W@    @�W@    B,_
  �Ba[v    B�#�@�[     @�[     B&�  �Bg��    B��@�^�    @�^�    B��  �Bxz�    B��w@�b�    @�b�    B�2  �By�V    B�o4@�f@    @�f@    BO�  �By)�    B�2�@�j     @�j     B*  �Br<�    B���@�m�    @�m�    BW�  �Bz2�    B��@�q�    @�q�    B��  �B��}    B�}�@�u@    @�u@    A�,  �B�x�    B�A@�y     @�y     A��  �B�LW    B�u@�|�    @�|�    A��"  �B��    B���@퀀    @퀀    A�7�  �B��    B��"@�@    @�@    A�d�  �B�?    B�Nf@�     @�     A�}  �B�J�    B��@��    @��    A�'�  �B�KV    B���@폀    @폀    B�D  �B{�    B���@�@    @�@    B\  �B���    B�[@�     @�     B�  �B��!    B�@��    @��    B�  �B��    B��@힀    @힀    B�d  �Bo��    B��@��@    @��@    B
��  �B���    B�f�@��     @��     B*�>  �BbX/    B�)�@���    @���    B�^  �BoT%    B��@���    @���    B  �By_�    B��z@��@    @��@    B��  �B}p�    B�r@@��     @��     B�  �B�h    B�4�@���    @���    B�K  �B{��    B���@���    @���    B,�  �B`��    B��\@��@    @��@    B$�W  �Bi�R    B�|�@��     @��     B*�?  �Bc\L    B�?�@���    @���    B$g  �Bi�    B�+@�ˀ    @�ˀ    B'_�  �BfA�    B�Ķ@��@    @��@    B%�  �Bhf    B��8@��     @��     B%ء  �Bg�T    B�I�@���    @���    B-��  �B_��    B�&@�ڀ    @�ڀ    B2݄  �BZ�    B�Β@��@    @��@    B7-�  �BVgt    B���@��     @��     B7^�  �BUq�    B�ST@���    @���    B0��  �B[�G    B��@��    @��    B3#�  �BX�    B���@��@    @��@    B6a  �BU
Q    B��D@��     @��     B; G  �BO �    B�\�@���    @���    BE�  �BE/�    B��@���    @���    BL'�  �B>�    B���@��@    @��@    Br5% �B�    B��(@�      @�      Bu�� �Bl�    B�eR@��    @��    By��  �B��    B�'u@��    @��    Bw�/  �BtY    B��@�@    @�@    By��  �B��    B���@�     @�     B�Q   �B.e    B�m�@��    @��    B�v�  �B o
    B�/�@��    @��    B�(�  �A�m     B���@�@    @�@    B��W  �A�z=    B���@�     @�     B�i�  �A�C    B�u�@�!�    @�!�    B��'  �A�x�    B�7�@�%�    @�%�    B���  �A��L    B���@�)@    @�)@    B��
  �A�R�    B���@�-     @�-     B��  �A�8~    B�}�@�0�    @�0�    B��%  �A��!    B�?x@�4�    @�4�    B�� �A��"    B�X@�8@    @�8@    B��& �A�Os    B��3@�<     @�<     B��=  �B
"~    B��	@�?�    @�?�    B�,s  �B[2    B�F�@�C�    @�C�    B�g  �AډU    B��@�G@    @�G@    B��� �A�Xp    B��s@�K     @�K     B��h �A�     B��8@�N�    @�N�    B��� �AŤ�    B�M�@�R�    @�R�    B��g �A��!    B��@�V@    @�V@    B�S �A��n    B��o@�Z     @�Z     B�J �A�p�    B��$@�]�    @�]�    B�� �A�3�    B�T�@�a�    @�a�    B�å �A��'    B��@�e@    @�e@    B�n� �A�G    B��,@�i     @�i     B�<C  �A�R=    B���@�l�    @�l�    B�@�  �A��K    B�[t@�p�    @�p�    B��  �A�o�    B�@�t@    @�t@    B�x  �A�_!    B�ޭ@�x     @�x     B�}�  �A���    B��E@�{�    @�{�    B���  �A�|    B�a�@��    @��    B��Y  �A���    B�#j@�@    @�@    B��S  �A�    B���@�     @�     B�t�  �A�gT    B���@��    @��    B�TB  �A�B    B�h	@    @    B�I  �A��~    B�)�@�@    @�@    B���  �B�Z    B��@�     @�     B},~  �B��    B���@��    @��    B|9�  �Bg�    B�n@    @    Bz��  �B�W    B�/~@�@    @�@    B}�F  �B�    B���@�     @�     B~ �  �BP4    B��e@��    @��    B���  �Bf    B�s�@    @    Bz'�  �B�s    B�5A@�@    @�@    Bs?�  �B"    B���@�     @�     Bi��  �B$    B��@��    @��    BZt�  �B2�    B�yv@    @    BN �  �B?N�    B�:�@�@    @�@    BE�8  �BG�    B��7@��     @��     B;��  �BQK�    B���@���    @���    B)�N  �Bb��    B�~�@�ʀ    @�ʀ    B}�  �Bu�    B�@F@��@    @��@    B�n  �Bs{    B��@��     @��     B��  �Brw�    B��@���    @���    BD  �Br�_    B|@�ـ    @�ـ    B&D  �Bms}    B~�@��@    @��@    B"��  �BfR�    B~�@��     @��     B�  �Bp�S    B}�C@���    @���    B(�  �Bp]�    B}�@��    @��    B
��  �B}    B|�\@��@    @��@    B�  �Bw�D    B|�@��     @��     A��B  �B��A    B{�c@���    @���    A�Z�  �B���    B{�@���    @���    A�MA  �B���    Bz�Z@��@    @��@    A�+J  �B�_    Bz!�@��     @��     B �B  �B���    By�@@��    @��    A�W� �B��i    By&�@��    @��    A�V  �B�-�    Bx�@�
@    @�
@    A˔Q  �B���    Bx+{@�     @�     A�j4  �B��    Bw��@��    @��    AՓm  �B�v    Bw0:@��    @��    A��.  �B��    Bv��@�@    @�@    A�k�  �B��    Bv4�@�     @�     A��   �B��:    Bu�<@� �    @� �    A��  �B��    Bu9�@�$�    @�$�    A�  �B��W    Bt��@�(@    @�(@    A�U  �B��^    Bt>@�,     @�,     A�+A  �B�n    Bs�c@�/�    @�/�    A�D?  �B�    BsB�@�3�    @�3�    A��/  �B���    Br��@�7@    @�7@    A�6F  �B�^0    BrG@�;     @�;     Aޑl  �B��o    Bq�S@�>�    @�>�    A��  �B���    BqK�@�B�    @�B�    A�Y�  �B�t�    Bp͸@�F@    @�F@    A��  �B� �    BpO�@�J     @�J     A�z:  �B�$�    Bo�@�M�    @�M�    A��)  �B��Y    BoT7@�Q�    @�Q�    A  �B�0    Bn�[@�U@    @�U@    A��Y  �B��    BnX}@�Y     @�Y     A�"I  �B���    Bmڛ@�\�    @�\�    A���  �B�@Z    Bm\�@�`�    @�`�    A��  �B�J    Bl��@�d@    @�d@    A���  �B���    Bl`�@�h     @�h     A�y�  �B���    Bk��@�k�    @�k�    A�&�  �B�T    Bke@�o�    @�o�    A�M!  �B��l    Bj�@�s@    @�s@    A�̹  �B�cV    Bji!@�w     @�w     A뾓  �B��#    Bi�)@�z�    @�z�    A��  �B���    Bim.@�~�    @�~�    A싑  �B���    Bh�1@�@    @�@    A�d.  �B���    Bhq1@�     @�     A��Y  �B��9    Bg�/@��    @��    A�9  �B��    Bgu*@    @    A��Q  �B��    Bf�#@�@    @�@    A�g�  �B��    Bfy@�     @�     A��x  �B�;0    Be�@��    @��    A�GO  �B���    Be|�@    @    A�A�  �B��o    Bd��@�@    @�@    A���  �B���    Bd��