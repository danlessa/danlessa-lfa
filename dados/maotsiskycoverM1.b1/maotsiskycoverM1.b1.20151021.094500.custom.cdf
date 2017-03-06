CDF  �   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.4-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20151021094500.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      62 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      237 pixels     center_y      320 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Zcreated by user dsmgr on machine ruby at 2015-10-21 11:49:00, using ingest-tsi-12.4-0.el6      ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20151021.094500.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:46 2017 UTC       ANDERS_armtime_timestamp      1488830686     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2015-10-21 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2015-10-21 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2015-10-21 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �V&ՀBH  �rdt�M�M@�#�    @�#�    ��  ����      =�$�@�'@    @�'@    ��  ����      >M/�@�+     @�+     ��  ����      >�g.@�.�    @�.�    ��  ����      >�6�@�2�    @�2�    ��  ����      ?�O@�6@    @�6@    ��  ����      ?0�m@�:     @�:     ��  ����      ?PS�@�=�    @�=�    ��  ����      ?o�,@�A�    @�A�    ��  ����      ?��f@�E@    @�E@    ��  ����      ?�F�@�I     @�I     ��  ����      ?��F@�L�    @�L�    ��  ����      ?���@�P�    @�P�    ��  ����      ?�dy@�T@    @�T@    ��  ����      ?�2@�X     @�X     ��  ����      ?�� @�[�    @�[�    ��  ����      ?���@�_�    @�_�    ��  ����      @��@�c@    @�c@    ��  ����      @
vp@�g     @�g     ��  ����      @P�@�j�    @�j�    ��  ����      @+�@�n�    @�n�    ��  ����      @"9@�r@    @�r@    ��  ����      @)��@�v     @�v     ��  ����      @1��@�y�    @�y�    ��  ����      @9�Y@�}�    @�}�    A��  �A޲�    @Aq @�@    @�@    A���  �A�)    @IK�@�     @�     A�$�  �A�M    @Q&�@��    @��    A���  �A�n    @Y�@ጀ    @ጀ    A�W  �A��t    @`ܜ@�@    @�@    A�sO  �A��y    @h��@�     @�     A� �  �B�    @p��@��    @��    A��\  �B%    @xm�@ᛀ    @ᛀ    A�w�  �B�    @�$S@�@    @�@    A�~T  �BE�    @��@�     @�     A�1�  �B)�    @��p@��    @��    A�PG  �B:�    @��@᪀    @᪀    A� Q  �B"��    @�ڝ@�@    @�@    A�US  �B'b!    @��:@�     @�     A���  �B'}J    @���@��    @��    A�>y  �B*�    @���@Ṁ    @Ṁ    A�J]  �B+;F    @��)@�@    @�@    A�7"  �B-F�    @�~�@��     @��     A�ٵ  �B+�#    @�l�@���    @���    A�sU  �B-Y�    @�Z;@�Ȁ    @�Ȁ    A��N  �B-l�    @�G�@��@    @��@    A�P�  �B,ob    @�5�@��     @��     A̳�  �B+�    @�#p@���    @���    A�Ր  �B,ls    @�4@�׀    @�׀    AǱ  �B,n�    @���@��@    @��@    A� ,  �B+�n    @���@��     @��     A��[  �B,)�    @�ږ@���    @���    A���  �B+�(    @��h@��    @��    A��  �B*�    @ζ>@��@    @��@    A�l�  �B*��    @Ҥ@��     @��     A�z�  �B)Y    @֑�@���    @���    A�A�  �B,&e    @��@���    @���    A�%  �B)�    @�m�@��@    @��@    A�	o  �B+��    @�[�@��     @��     A�#T  �B,��    @�I�@� �    @� �    A��  �B,^�    @�7u@��    @��    A�j6  �B+�     @�%e@�@    @�@    A׽~  �B.�G    @�X@�     @�     A�_@  �B299    @�O@��    @��    A��(  �B4Ƀ    @��I@��    @��    A�D  �B9�)    @��E@�@    @�@    AȒ  �B=��    A �@�     @�     A�@�  �BC%*    Aܤ@��    @��    Aüd  �BF��    Aӧ@�"�    @�"�    A�9�  �BG�Q    Aʫ@�&@    @�&@    A�\�  �BMw�    A��@�*     @�*     A�m�  �BL�m    A
��@�-�    @�-�    A�Yj  �BOU    A��@�1�    @�1�    A�:  �BN�    A��@�5@    @�5@    A�R�  �BOX�    A��@�9     @�9     A��s  �BI�    A��@�<�    @�<�    A�L  �BKEP    A��@�@�    @�@�    A�M�  �BJ     A��@�D@    @�D@    A��  �BHh�    Az@�H     @�H     B;.  �BJ�c    Aq@�K�    @�K�    B|  �BO?�    Ah1@�O�    @�O�    B�M  �BV	L    A_E@�S@    @�S@    Bwx  �B]�    A VZ@�W     @�W     B�  �Bh��    A"Mp@�Z�    @�Z�    B�V  �Bu     A$D�@�^�    @�^�    B .  �B���    A&;�@�b@    @�b@    A�   �B��
    A(2�@�f     @�f     A�s�  �B�f�    A*)�@�i�    @�i�    A�C  �B�b2    A, �@�m�    @�m�    A�5�  �B��    A.@�q@    @�q@    A�  �B���    A0)@�u     @�u     Ab��  �B�#�    A2H@�x�    @�x�    AU�_  �B��    A3�g@�|�    @�|�    AI�  �B�dx    A5�@�@    @�@    AJv�  �B�9�    A7�@�     @�     AYŉ  �B�X!    A9��@��    @��    Ax��  �B�w9    A;��@⋀    @⋀    A�)�  �B���    A=�@�@    @�@    A��=  �B�L$    A?�7@�     @�     A���  �B���    AA�]@��    @��    A�	p  �B���    AC��@⚀    @⚀    A��  �B��J    AE��@�@    @�@    A��  �B��i    AG��@�     @�     B�d  �B��    AI��@��    @��    B�  �Bw&�    AK�&@⩀    @⩀    B4|  �Bt�~    AM�P@�@    @�@    B��  �Bs�    AO�{@�     @�     BH�  �Bm��    AQx�@��    @��    B
�* �Br�.    ASo�@⸀    @⸀    B�� �Bd��    AUg @�@    @�@    BW! �B\0-    AW^.@��     @��     B#��  �BP�L    AYU\@���    @���    B�u �BKD;    A[L�@�ǀ    @�ǀ    B�Z �BC�$    A]C�@��@    @��@    B��  �B;
�    A_:�@��     @��     B�  �B/�    Aa2@���    @���    BS1  �B&ݻ    Ac)K@�ր    @�ր    B;  �B&�    Ae |@��@    @��@    B 3}  �B�    Ag�@��     @��     BN�  �B�+    Ai�@���    @���    A�!�  �Bm    Ak@��    @��    A�Y�  �B�    Al�F@��@    @��@    A�|  �B�    An�y@��     @��     A���  �B!    Ap�@���    @���    A��X  �BEq    Ar��@��    @��    A���  �B��    At�@��@    @��@    B�w  �B�7    Av�J@��     @��     B�(  �B<?    Ax�@���    @���    B,j  �B��    Az��@��    @��    B1*  �B�e    A|��@�@    @�@    Bs�  �Bm%    A~�@�     @�     B~i  �BI0    A�R�@��    @��    Byy  �B"�    A�NE@��    @��    B�2  �B${�    A�I�@�@    @�@    B ��  �B%�^    A�E{@�     @�     A�M  �B&c�    A�A@��    @��    A��t  �B)ߌ    A�<�@�!�    @�!�    A��{  �B(3�    A�8M@�%@    @�%@    A։�  �B"Xq    A�3�@�)     @�)     Aԣ�  �BIM    A�/�@�,�    @�,�    Aޭ�  �B~.    A�+@�0�    @�0�    A�e]  �BM�    A�&�@�4@    @�4@    A�{�  �BEq    A�"U@�8     @�8     A��  �B*    A��@�;�    @�;�    A���  �B��    A��@�?�    @�?�    A��  �B׫    A�&@�C@    @�C@    A��  �B#�    A��@�G     @�G     A��$  �B'�    A�\@�J�    @�J�    A�-  �B+�=    A��@�N�    @�N�    A�w`  �B0��    A��@�R@    @�R@    A�q�  �B4�f    A��-@�V     @�V     A��8  �B8K�    A���@�Y�    @�Y�    A�f&  �B;��    A��b@�]�    @�]�    A���  �B@�G    A���@�a@    @�a@    A���  �BD�    A��@�e     @�e     A���  �BI�    A��0@�h�    @�h�    A�$�  �BJ�    A���@�l�    @�l�    A�<  �BL��    A��c@�p@    @�p@    A���  �BO��    A���@�t     @�t     A�`0  �BU��    A�ז@�w�    @�w�    A���  �B]��    A��/@�{�    @�{�    A���  �BhT�    A���@�@    @�@    A�+  �Blq�    A��`@�     @�     A�:  �Bvn�    A���@��    @��    A�7�  �B~�!    A���@㊀    @㊀    A�L�  �B���    A��(@�@    @�@    AЦi  �B�:    A���@�     @�     A�.�  �B��i    A��W@��    @��    A��u  �B�}�    A���@㙀    @㙀    A��Z  �B��V    A���@�@    @�@    A���  �B���    A��@�     @�     A�J� �B��    A���@��    @��    A�!z �B��    A��F@㨀    @㨀    A�F �B�9    A���@�@    @�@    A�  �B���    A��p@�     @�     A��; �B��6    A��@��    @��    A�s� �B�
�    A���@㷀    @㷀    A�4  �B��    A��,@�@    @�@    A�wK  �B�ܝ    A���@�     @�     A̷�  �B��    A�R@���    @���    A�~w  �B��j    A�z�@�ƀ    @�ƀ    AޜJ  �B�0�    A�vv@��@    @��@    A��  �B��    A�r@��     @��     A��  �B�j�    A�m�@���    @���    B;� �B|��    A�i)@�Հ    @�Հ    B� �Bw{�    A�d�@��@    @��@    BB� �Bt��    A�`I@��     @��     B� �Bw��    A�[�@���    @���    B�x �Byr�    A�Wf@��    @��    B
�` �By��    A�R�@��@    @��@    B�= �Bz��    A�N�@��     @��     B�i �B{�    A�J@���    @���    B�� �By^�    A�E�@��    @��    Big �Bxm_    A�A&@��@    @��@    B�m �Bz�    A�<�@��     @��     B�B �Bz`/    A�8;@���    @���    B޼ �B{�q    A�3�@��    @��    B*F �BP�    A�/N@�@    @�@    B�� �B��    A�*�@�
     @�
     B�s �B�.L    A�&^@��    @��    B[p �B��    A�!�@��    @��    BZ� �B��    A�k@�@    @�@    A�� �B�F�    A��@�     @�     A�M �B�k�    A�v@��    @��    A�  �B��\    A��@� �    @� �    A�H� �B�VM    A�}@�$@    @�$@    A��@ �B��k    A� @�(     @�(     A�XU �B�7L    A��@�+�    @�+�    A�� �B�>    A��@�/�    @�/�    A�c �B��Y    A���@�3@    @�3@    A�� �B���    A��@�7     @�7     A�c�  �B���    A���@�:�    @�:�    BV�  �B���    A���@�>�    @�>�    B	�z  �B��i    A��|@�B@    @�B@    B:�  �B~V�    A���@�F     @�F     BK� �B)�    A��t@�I�    @�I�    BC �B���    A���@�M�    @�M�    B�l �B�En    A��h@�Q@    @�Q@    A�~H �B�a2    A���@�U     @�U     A�]�  �B���    A��X@�X�    @�X�    A���  �B���    A���@�\�    @�\�    A���  �B�Y�    A��D@�`@    @�`@    A�j�  �B�4�    Aپ�@�d     @�d     A�g� �B�G�    Aں-@�g�    @�g�    A��O  �B��    A۵�@�k�    @�k�    A���  �B��    Aܱ@�o@    @�o@    Aؙ  �B�W{    Aݬ�@�s     @�s     A�'q  �B�)�    Aާ�@�v�    @�v�    A֩	  �B���    Aߣ`@�z�    @�z�    A��J  �B���    A���@�~@    @�~@    A�('  �B���    A�:@�     @�     A���  �B�X/    A╥@��    @��    A�G�  �B���    A�@䉀    @䉀    A�O  �B�}�    A�x@�@    @�@    A�,  �B��    A��@�     @�     A��y  �B��    A�G@��    @��    Asy  �B�:&    A�~�@䘀    @䘀    AL�  �B���    A�z@�@    @�@    A>�  �B��l    A�ut@�     @�     A w�  �B�y
    A�p�@��    @��    A"�  �B�H�    A�l7@䧀    @䧀    A�   �B�<�    A�g�@�@    @�@    A�0  �B���    A�b�@�     @�     A#  �B�Q    A�^R@��    @��    AƖ  �B��G    A�Y�@䶀    @䶀    Az�  �B���    A�U@�@    @�@    A&��  �B���    A�P`@�     @�     A.	B  �B��    A�K�@���    @���    A>�  �B�ý    A�G@�ŀ    @�ŀ    AT�  �B���    A�Bc@��@    @��@    A_5�  �B��    A�=�@��     @��     Am�  �B���    A�9	@���    @���    A`�9 �B���    A�4Y@�Ԁ    @�Ԁ    Agb� �B���    A�/�@��@    @��@    Asg� �B�;�    A�*�@��     @��     AsY �B�D�    A�&B@���    @���    Az | �B�n!    A�!�@��    @��    A�*� �B�o�    A��@��@    @��@    A�F� �B�(2    A�@��     @��     A�� �B�į    A�d@���    @���    A���  �B���    A��@��    @��    A�"�  �B�	6    B �@��@    @��@    A֍w  �B��Y    B ��@��     @��     A�  �B�}F    B 7@���    @���    B �B��i    B}�@��    @��    B�7 �B|U8    B�t@�@    @�@    B� �Bp��    By@�	     @�	     B)w� �Be��    B��@��    @��    B2� �B\�    BtJ@��    @��    B9?  �BU:     B��@�@    @�@    B?��  �BN��    Bo~@�     @�     BB��  �BJH    B�@��    @��    B?�b  �BLE�    Bj�@��    @��    B=x  �BM!�    B�G@�#@    @�#@    B9�~  �BP>�    Be�@�'     @�'     B8/�  �BR��    B�s@�*�    @�*�    B9�Y  �BQ8    Ba@�6     @�6     B3�a  �BW��    Bٿ@�9�    @�9�    B0�0  �BZ�"    B	WP@�=�    @�=�    B-Ү �B]��    B	��@�A@    @�A@    B+s,  �B`�/    B
Ro@�E     @�E     B&m[  �Bg��    B
��@�H�    @�H�    B�� �BtP�    BM�@�L�    @�L�    B�� �BvT8    B�@�P@    @�P@    B�� �Byp    BH�@�T     @�T     BV� �Bz@9    B�)@�W�    @�W�    B{� �BtF    BC�@�[�    @�[�    B�� �Bu��    B�9@�_@    @�_@    BdO �Bv$�    B>�@�c     @�c     B�x �B}�    B�E@�f�    @�f�    B�� �B��    B9�@�j�    @�j�    B	� �B��2    B�M@�n@    @�n@    B�L �B�    B4�@�r     @�r     B �  �B�l    B�P@�u�    @�u�    A�$� �B�y	    B/�@�y�    @�y�    B�@  �B��\    B�N@�}@    @�}@    B�T �B�ڂ    B*�@�     @�     B� �B���    B�H@��    @��    Bu �B�;�    B%�@刀    @刀    B#> �B��z    B�=@�@    @�@    B;� �B�$N    B �@�     @�     Bo  �B}r�    B�-@��    @��    BRD �Bs(    B�@嗀    @嗀    B?Q �BxT�    B�@�@    @�@    BJN �By8    B�@�     @�     B:� �ByJr    B��@��    @��    B�' �Bs�s    Bp@妀    @妀    B��  �Br��    B��@�@    @�@    B�\ �Bo��    BO@�     @�     B$t�  �Bj��    B��@��    @��    B$�t �Bjt1    B(@嵀    @嵀    B�z �Bp��    B��@�@    @�@    B�9 �Brǫ    B�@�     @�     BI� �Br    Bd@���    @���    B
� �BtR�    B��@�Ā    @�Ā    B��  �Bt��    Bz0@��@    @��@    B�  �Bt�0    B��@��     @��     Be�  �By�\    Bt�@���    @���    B�  �B|�    B�W@�Ӏ    @�Ӏ    BM�  �Bw��    Bo�@��@    @��@    B�% �B{��    B�@��     @��     B��  �BM    Bjq@���    @���    Bj�  �B{��    B��@��    @��    B2F  �ByVF    Be%@��@    @��@    B0�  �Bx�]    B�}@��     @��     B  �B��    B _�@���    @���    B]�  �B�_�    B �(@��    @��    B�&  �B�D�    B!Z|@��@    @��@    A�>  �B���    B!��@��     @��     A�c  �B�$�    B"U@���    @���    A��d  �B�r^    B"�l@� �    @� �    A�|V  �B��=    B#O�@�@    @�@    A���  �B��-    B#�@�     @�     A��  �B�5�    B$JN@��    @��    A�t�  �B��=    B$ǖ@��    @��    A�F�  �B�d�    B%D�@�@    @�@    A��  �B��!    B%�!@�     @�     A��  �B�8�    B&?d@��    @��    Aт  �B�W    B&��@��    @��    A�)�  �B���    B'9�@�"@    @�"@    A��W  �B��p    B'�"@�&     @�&     A��  �B���    B(4^@�)�    @�)�    A�3�  �B��c    B(��@�-�    @�-�    Aٞ�  �B�F�    B).�@�1@    @�1@    A�om �B��    B)�@�5     @�5     A��`  �B�e�    B*);@�8�    @�8�    A���  �B��g    B*�m@�<�    @�<�    B��  �B���    B+#�@�@@    @�@@    B�'  �B���    B+��@�D     @�D     B�K  �B�j�    B,�@�G�    @�G�    B��  �B��    B,�%@�K�    @�K�    B�|  �B}�`    B-N@�O@    @�O@    B�� �B��y    B-�u@�S     @�S     B�  �B~��    B.�@�V�    @�V�    BgQ �Bx,,    B.��@�Z�    @�Z�    B\ �B|!X    B/�@�^@    @�^@    B8� �Bv;I    B/��@�b     @�b     B#c  �Bh�9    B0@�e�    @�e�    B"|<  �BkM    B0�5@�i�    @�i�    B `�  �Bku    B1N@�m@    @�m@    BI�  �Bq5v    B1~d@�q     @�q     B��  �Br�    B1�y@�t�    @�t�    B�G  �Btx    B2x�@�x�    @�x�    B Z �Br�    B2��@�|@    @�|@    B��  �BsL9    B3r�@�     @�     Bԑ  �BtG�    B3�@��    @��    B�_  �Bp�#    B4l�@懀    @懀    BX�  �Bp    B4��@�@    @�@    B�y  �Bs�    B5f�@�     @�     B��  �Bw    B5��@��    @��    B�� �B~��    B6`�@斀    @斀    B �B���    B6��@�@    @�@    B�, �B��$    B7Z�@�     @�     B��  �B�    B7׾@��    @��    B��  �B~�    B8T�@楀    @楀    B%�  �B��&    B8ѩ@�@    @�@    A��~  �B�`�    B9N�@�     @�     A�a�  �B���    B9ˉ@��    @��    A���  �B�@�    B:Hv@洀    @洀    A��b  �B�!�    B:�`@�@    @�@    A���  �B���    B;BG@�     @�     B-  �B�    B;�,@��    @��    Bn  �B�t.    B<<@�À    @�À    BGn  �B�xw    B<��@��@    @��@    B�  �B�$�    B=5�@��     @��     A�  �B�b�    B=��@���    @���    A���  �B�d�    B>/|@�Ҁ    @�Ҁ    A���  �B�b�    B>�P@��@    @��@    B c  �B�ƶ    B?)"