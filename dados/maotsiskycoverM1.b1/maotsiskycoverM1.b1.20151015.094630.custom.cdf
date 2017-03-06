CDF  s   
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.4-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20151015094630.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      62 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.013625       ccd_width_factor      	0.013312       center_x      240 pixels     center_y      317 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Zcreated by user dsmgr on machine ruby at 2015-10-15 11:49:00, using ingest-tsi-12.4-0.el6      ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20151015.094630.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:45 2017 UTC       ANDERS_armtime_timestamp      1488830685     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2015-10-15 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2015-10-15 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2015-10-15 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �V�BH  �rdt�M�M@�.�    @�.�    ��  ����      =�@�2�    @�2�    ��  ����      >"�e@�6@    @�6@    ��  ����      >��@�:     @�:     ��  ����      >���@�=�    @�=�    ��  ����      ?�a@�A�    @�A�    ��  ����      ?'+�@�E@    @�E@    ��  ����      ?F��@�I     @�I     ��  ����      ?fd&@�L�    @�L�    ��  ����      ?� W@�P�    @�P�    ��  ����      ?�ά@�T@    @�T@    ��  ����      ?��@�X     @�X     ��  ����      ?�k�@�[�    @�[�    ��  ����      ?�:@�_�    @�_�    ��  ����      ?��@�c@    @�c@    ��  ����      ?��R@�g     @�g     ��  ����      ?�@�j�    @�j�    ��  ����      @ �k@�n�    @�n�    ��  ����      @��@�r@    @�r@    ��  ����      @�M@�v     @�v     ��  ����      @p�@�y�    @�y�    ��  ����      @ XP@�}�    @�}�    ��  ����      @(?�@�@    @�@    ��  ����      @0'r@�     @�     ��  ����      @8@��    @��    ��  ����      @?��@ጀ    @ጀ    @�x�  �A��    @G�_@�@    @�@    @���  �A��    @O�@�     @�     @�  �A�
    @W��@��    @��    @哌  �A��@    @_��@ᛀ    @ᛀ    @�̟  �A��    @g}[@�@    @�@    @�l�  �A��    @oe,@�     @�     @��  �A��    @wM@��    @��    @ف  �A�pZ    @4�@᪀    @᪀    @���  �A��    @��g@�@    @�@    @ׇ�  �A�~j    @��^@�     @�     @�x�  �A�_�    @�vY@��    @��    @��r  �A�u    @�jW@Ṁ    @Ṁ    @�\'  �A��    @�^Y@�@    @�@    @�|�  �A��    @�R_@��     @��     @�z_  �A���    @�Fh@���    @���    @�|  �A��`    @�:t@�Ȁ    @�Ȁ    @އ�  �A�J�    @�.�@��@    @��@    @ۧ�  �A���    @�"�@��     @��     @��  �A�~=    @��@���    @���    @�Oe  �A���    @�
�@�׀    @�׀    @��  �A�̢    @���@��@    @��@    @޺�  �A�ו    @��@��     @��     @�r�  �A��    @��'@���    @���    @��O  �A�N�    @��N@��    @��    @�s�  �A�eh    @��x@��@    @��@    @���  �A���    @�å@��     @��     @�g�  �A�x�    @ʷ�@���    @���    @���  �A�w    @ά	@���    @���    @�{O  �A��4    @Ҡ?@��@    @��@    @�8#  �A�&    @֔x@��     @��     @��  �A��m    @ڈ�@� �    @� �    @�3  �A�;5    @�|�@��    @��    @��  �A�b�    @�q7@�@    @�@    @��  �A�@H    @�e|@�     @�     @�*�  �A�
�    @�Y�@��    @��    @�  �A�T�    @�N@��    @��    @���  �A��3    @�B^@�@    @�@    @�p �A�$^    @�6�@�     @�     @뉘 �A�S�    @�+@��    @��    @��� �A�o�    @�Y@�"�    @�"�    @��= �A���    A	�@�&@    @�&@    @� �A�6     A@�*     @�*     @�� �A�sL    A�6@�-�    @�-�    @�� �A��    A�g@�1�    @�1�    @�&B �A���    A�@�5@    @�5@    @�C �A�$�    A
��@�9     @�9     @�/  �A��z    A�@�<�    @�<�    @�  �A�ώ    A�7@�@�    @�@�    @��� �A�ʉ    A�n@�D@    @�D@    @�~ �A���    Aէ@�H     @�H     @��� �A�ŧ    A��@�K�    @�K�    Ai� �A���    A�@�O�    @�O�    A �� �A��    A�X@�S@    @�S@    @�[, �A�w�    A��@�W     @�W     A� �A���    A��@�Z�    @�Z�    A�v �A���    A�@�^�    @�^�    A8� �A��    A �U@�b@    @�b@    A|� �A�ǥ    A"��@�f     @�f     A[� �A�`k    A$��@�i�    @�i�    A�L �A��     A&�@�m�    @�m�    A.( �A��    A(�e@�q@    @�q@    A�� �A��*    A*��@�u     @�u     A�! �A�D�    A,��@�x�    @�x�    A� �A�%�    A.�<@�|�    @�|�    A� �A���    A0�@�@    @�@    A	�� �A��U    A2y�@�     @�     A
 �A�z�    A4t@��    @��    A_� �A�k     A6nk@⋀    @⋀    A� �A��K    A8h�@�@    @�@    A
�� �A�R�    A:c@�     @�     A� �A���    A<]X@��    @��    AH� �A���    A>W�@⚀    @⚀    A�� �A�U    A@Q�@�@    @�@    An �A�V    ABLO@�     @�     A	�� �A�M�    ADF�@��    @��    A�� �A�Xv    AF@�@⩀    @⩀    A� �A��f    AH;N@�@    @�@    AJL �A��-    AJ5�@�     @�     AV� �A���    AL/�@��    @��    A�� �A�#    AN*V@⸀    @⸀    AT �A�e    AP$�@�@    @�@    A2� �A��    AR
@��     @��     A�� �A��    ATe@���    @���    A� �A��3    AV�@�ǀ    @�ǀ    A�B �A��    AX@��@    @��@    A� �A�    AZ|@��     @��     A	 �A�p6    A\�@���    @���    ACM �A��    A]�:@�ր    @�ր    A�v �A�}n    A_��@��@    @��@    AҜ �A��
    Aa��@��     @��     Av� �A��    Ac�]@���    @���    AXT �A���    Ae��@��    @��    A�� �A�8    Ag�#@��@    @��@    A� �A���    Aiۇ@��     @��     A( �A��    Ak��@���    @���    AN �A��B    Am�Q@��    @��    Aѯ �A�H�    Aoʷ@��@    @��@    A�- �A��\    Aq�@��     @��     A�  �A��l    As��@���    @���    Aղ �A�o    Au��@��    @��    A�\ �A���    Aw�U@�@    @�@    A�  �A�:W    Ay��@�     @�     AS �A���    A{�)@��    @��    A<	 �A��    A}��@��    @��    A�� �A���    A��@�@    @�@    AYh �A���    A��5@�     @�     A�� �A��J    A��k@��    @��    Aq �A���    A�ơ@�!�    @�!�    A�7 �A���    A���@�%@    @�%@    A � �A�p6    A��@�)     @�)     A�� �A�2    A��F@�,�    @�,�    A�� �A���    A��}@�0�    @�0�    A�, �A�2    A���@�4@    @�4@    A�� �A�Ra    A���@�8     @�8     A";� �A�~p    A��%@�;�    @�;�    A"λ �A�!�    A��]@�?�    @�?�    A$@� �A��D    A���@�C@    @�C@    A"�+ �A��    A���@�G     @�G     A'm �A�    A��@�J�    @�J�    A%� �A�E�    A��A@�N�    @�N�    A(к �A���    A��z@�R@    @�R@    A(�+ �A��K    A���@�V     @�V     A(�� �A�V�    A���@�Y�    @�Y�    A)W �A�e�    A��'@�]�    @�]�    A%� �A�N�    A��a@�a@    @�a@    A$Te �A���    A���@�e     @�e     A(5 �A�s5    A���@�h�    @�h�    A(�� �A�<b    A��@�l�    @�l�    A)� �A�d    A��K@�p@    @�p@    A/� �A��    A���@�t     @�t     A+�v �A�G�    A���@�w�    @�w�    A.�? �A�    A���@�{�    @�{�    A1
 �A�[�    A��7@�@    @�@    A,$o �A���    A�~s@�     @�     A+H* �A���    A�{�@��    @��    A.м �A�E0    A�x�@㊀    @㊀    A.�; �A�Sb    A�v%@�@    @�@    A2�7 �A�],    A�sa@�     @�     A/W� �A���    A�p�@��    @��    A3g� �A�pB    A�m�@㙀    @㙀    A,�h �A��    A�k@�@    @�@    A+�A �A�G%    A�hP@�     @�     A2�� �A�v�    A�e�@��    @��    A3�� �A���    A�b�@㨀    @㨀    A1� �A��%    A�`@�@    @�@    A+�C �A�e}    A�]@@�     @�     A3M �A���    A�Z|@��    @��    A/�A �A��P    A�W�@㷀    @㷀    A-� �A���    A�T�@�@    @�@    A/�= �A�֎    A�R1@�     @�     A/�� �A�D�    A�Om@���    @���    A0Ti �A�t+    A�L�@�ƀ    @�ƀ    A4)V �A�s�    A�I�@��@    @��@    A.ً �A���    A�G!@��     @��     A1 �A���    A�D]@���    @���    A3� �A��.    A�A�@�Հ    @�Հ    A1� �A��R    A�>�@��@    @��@    A.�� �A��r    A�<@��     @��     A-� �A���    A�9M@���    @���    A5lf �A��j    A�6�@��    @��    A4 9 �A��    A�3�@��@    @��@    A7�u �A��o    A�1@��     @��     A96 �A�-�    A�.<@���    @���    A:�d �A�0�    A�+x@��    @��    A>�� �A�?�    A�(�@��@    @��@    A?( �A��8    A�%�@��     @��     A>�< �A��g    A�#*@���    @���    AB2 �A��'    A� e@��    @��    ADJ� �A��l    A��@�@    @�@    AC$� �A�V�    A��@�
     @�
     AA  �A�	    A�@��    @��    ACN �A�[�    A�P@��    @��    A<+N �A�e{    A��@�@    @�@    A@�� �A��%    A��@�     @�     A@O� �A���    A��@��    @��    A?�O �A��i    A�
8@� �    @� �    A?�� �A�Z�    A�r@�$@    @�$@    AA�T �A�9o    A��@�(     @�(     AB<� �A���    A��@�+�    @�+�    AE�s �A�WK    A��@�/�    @�/�    A@8	 �A�$b    A��W@�3@    @�3@    AG:~ �A�d�    A���@�7     @�7     AC�: �A��&    A���@�:�    @�:�    AP�  �A��    A�� @�>�    @�>�    AD�  �A���    A��8@�B@    @�B@    AO�� �A��>    A��p@�F     @�F     AN1 �A�4�    A��@�I�    @�I�    AM3� �A��
    A���@�M�    @�M�    AJ� �A��Y    A��@�Q@    @�Q@    AO� �A�yC    A��L@�U     @�U     AN6 �A���    A���@�X�    @�X�    AM�� �A���    A�ݹ@�\�    @�\�    AO5� �A�y0    A���@�`@    @�`@    AD`� �A�	�    A��$@�d     @�d     AHJ� �A}�6    A��Y@�g�    @�g�    A?�� �A~��    A�Ҏ@�k�    @�k�    A<, �A��?    A���@�o@    @�o@    AA
: �A�    A���@�s     @�s     AE�X �A%<    A��+@�v�    @�v�    AD�` �A���    A��^@�z�    @�z�    A>rl �Ay�C    A�đ@�~@    @�~@    A9! �A{��    A���@�     @�     A?H� �A}]�    A��@��    @��    AAbW �A~޿    A�(@䉀    @䉀    AE�� �A}{X    A�Z@�@    @�@    AD�  �A|��    A㶋@�     @�     AFE� �A{�t    A䳼@��    @��    AD�� �Ax��    A��@䘀    @䘀    A?B� �Axz�    A�@�@    @�@    A?i% �Ay�    A�L@�     @�     AC �A|x�    A�{@��    @��    AF�� �A|�    A饩@䧀    @䧀    AI � �A|�N    A��@�@    @�@    AP&� �A{;�    A�@�     @�     AJI� �Ax�    A�2@��    @��    AI�� �Ax��    A�_@䶀    @䶀    AOj� �Az�    A@�@    @�@    AI�^ �A|�    A@�     @�     AG�| �A{��    A��@���    @���    AP�� �Az<�    A�@�ŀ    @�ŀ    AS�� �A{��    A�7@��@    @��@    AW'( �Ay�    A�`@��     @��     APM� �Av�    A�@���    @���    ANش �Ao�v    A���@�Ԁ    @�Ԁ    AP*y �Aq�     A���@��@    @��@    AO� �Ap�    A�~@��     @��     AJ�Y �As�>    A�{(@���    @���    ANt� �Au�L    A�xN@��    @��    AM�� �AvȲ    A�us@��@    @��@    AL/ �Aw�*    A�r�@��     @��     AM�P �At��    A�o�@���    @���    ARm �Aw�    A�l�@��    @��    AW�� �AvΓ    A�j@��@    @��@    AX2n �Av��    A�g&@��     @��     AZ� �Avу    B 2$@���    @���    AVn �Avyf    B ��@��    @��    AX�� �Ao�r    B/D@�@    @�@    A`(� �Aq�p    B��@�	     @�	     AeL� �Ap��    B,d@��    @��    A_�� �Aro�    B��@��    @��    A[̭ �Ap�K    B)�@�@    @�@    Aa�F �Ar�    B�@�     @�     A_ �Ar�)    B&�@��    @��    AdG. �At̖    B�,@��    @��    Af> �As��    B#�@�#@    @�#@    Ac;� �As�m    B�F@�'     @�'     A`�� �As�    B �@�*�    @�*�    Acg� �At�    B�_@�.�    @�.�    Aim �Ap�[    B�@�2@    @�2@    Ae�n �An�-    B�v@�6     @�6     AZ�J �Aj�    B@�9�    @�9�    A\4� �Al��    B��@�=�    @�=�    Aa�� �An�-    B	@�A@    @�A@    Abf� �Ao`�    B	��@�E     @�E     Ag� �An��    B
(@�H�    @�H�    Ab� �AiC�    B
��@�L�    @�L�    Agז �Am�    B9@�P@    @�P@    Ab� �Agڇ    B��@�T     @�T     A`� �Af��    BH@�W�    @�W�    Aa�_ �Af�|    B��@�[�    @�[�    A`�� �Af-�    BV@�_@    @�_@    A`�� �Adt�    B��@�c     @�c     A]4 �Af0�    B	b@�f�    @�f�    Aeq� �Aa�N    B��@�j�    @�j�    Ab�x �AZu�    Bk@�n@    @�n@    Ac�1 �AX�    B��@�r     @�r     Aa �AV�    Bs@�u�    @�u�    A]]H �AT�    B��@�y�    @�y�    A]t� �AS��    B y@�}@    @�}@    AZA �A5~�    B~�@�     @�     AmT� �AS    B�|@��    @��    Agp� �AS��    B{�@刀    @刀    A`Ԉ �AQg�    B�~@�@    @�@    Aaj\ �ARX�    Bx�@�     @�     Ad�� �AR�I    B�}@��    @��    Ai�� �ARʉ    Bu�@嗀    @嗀    Aof6 �AM��    B�z@�@    @�@    Aq� �AO�    Br�@�     @�     AtL� �APpO    B�u@��    @��    Ay�d �AL+    Bo�@妀    @妀    Av�e �AOy    B�n@�@    @�@    Avw� �AI[�    Bl�@�     @�     A{#y �AF�#    B�d@��    @��    A��� �AHd�    Bi�@嵀    @嵀    A�e� �AG��    B�W@�@    @�@    A��g �AG��    Bf�@�     @�     A~�x �AC�6    B�I@���    @���    A�V� �AD��    Bc�@�Ā    @�Ā    A� k �AD~+    B�7@��@    @��@    A��� �AE:.    B`�@��     @��     A�Ғ �AF��    B�$@���    @���    A�� �AG3v    B]�@�Ӏ    @�Ӏ    A��W �AH׍    B�@��@    @��@    A��� �AI �    BZ�@��     @��     A��� �AM%0    B��@���    @���    A��p �AM��    BWf@��    @��    A�I� �AOѬ    B��@��@    @��@    A��m �AO�,    BTI@��     @��     A�Q� �AN��    Bҹ@���    @���    A��� �AH\+    B Q)@��    @��    A�jD �AGn8    B ϗ@��@    @��@    A�L� �AIjp    B!N@��     @��     A��� �AJ�<    B!�s@���    @���    A�\+ �AKjY    B"J�@� �    @� �    A��! �AG�    B"�K@�@    @�@    A��� �AK#�    B#G�@�     @�     A�uh �AI��    B#� @��    @��    A�� �AJ2�    B$D�@��    @��    A�ځ �AL	    B$��@�@    @�@    A��n �AOq�    B%AZ@�     @�     A��� �AO��    B%��@��    @��    A�l �AM $    B&>(@��    @��    A�O3 �AL�    B&��@�"@    @�"@    A�J! �AG�%    B':�@�&     @�&     A��� �AE�    B'�V@�)�    @�)�    A�� �AC��    B(7�@�-�    @�-�    A��� �AF�    B(�@�1@    @�1@    A� �ABu    B)4|@�5     @�5     A��� �A@6�    B)��@�8�    @�8�    A�fY �ACg`    B*1;@�<�    @�<�    A��� �ABU    B*��@�@@    @�@@    A�� �AAݳ    B+-�@�D     @�D     A��� �AC>    B+�T@�G�    @�G�    A�n[ �A=�    B,*�@�K�    @�K�    A�)Q �A?�p    B,�@�O@    @�O@    A�Y! �A:v�    B-'e@�S     @�S     A�G� �A:��    B-��@�V�    @�V�    A�6� �A:�    B.$@�Z�    @�Z�    A�;� �A9�h    B.�m@�^@    @�^@    A�� �A8�.    B/ �@�b     @�b     A�Vj �A8b�    B/�@�e�    @�e�    A��� �A6��    B0l@�i�    @�i�    A�Da �A8_�    B0��@�m@    @�m@    A�� �A;��    B1@�q     @�q     A��� �A>T{    B1�b@�t�    @�t�    A�� �A8��    B2�@�x�    @�x�    A�0R �A7�g    B2� @�|@    @�|@    A�R� �A7�    B3N@�     @�     A�zl �A3��    B3��@��    @��    A�d� �A3	�    B4�@懀    @懀    A��� �A0w    B4�1@�@    @�@    A��f �A.+�    B5z@�     @�     A��= �A0t�    B5��@��    @��    A�nD �A,��    B6	
@斀    @斀    A��� �A*��    B6�P@�@    @�@    A��o �A)��    B7�@�     @�     A��� �A(    B7��@��    @��    A�09 �A"kL    B8@楀    @楀    A�89 �A#b    B8�\@�@    @�@    A��� �A#�h    B8��@�     @�     A��� �A#��    B9|�@��    @��    AÉ� �A&�.    B9�@洀    @洀    A�4 �A'�G    B:yU@�@    @�@    A�Ͻ �A*�    B:��@�     @�     A�h� �A+^�    B;u�@��    @��    A�8K �A&�    B;�@�À    @�À    AĽ �A&��    B<r:@��@    @��@    AÀ' �A*jL    B<�p@��     @��     Aʃ� �A+��    B=n�@���    @���    A�� �A+=�    B=��@�Ҁ    @�Ҁ    A� �A/�C    B>k
@��@    @��@    A�u �A2�    B>�;@��     @��     Aȴ� �A1��    B?gj@���    @���    A�&� �A3
�    B?�@��    @��    A֣� �A2gk    B@c�@��@    @��@    A�U  �A'*�    B@��@��     @��     A�rk �A!��    BA`@���    @���    A�Ƚ �AI�    BA�A@���    @���    A�ұ �A��    BB\h@��@    @��@    A� �A��    BBڍ@��     @��     A�Y" �A yj    BCX�@���    @���    A�{? �A7�    BC��@���    @���    A��Q �A�    BDT�@�@    @�@    A��  �@�[�    BD�@�     @�     A�� �@���    BEQ1@�
�    @�
�    A�d� �@��    BE�M@��    @��    A��m �@�PK    BFMg@�@    @�@    A��p �@�=8    BFˀ@�     @�     A�i� �@똲    BGI�@��    @��    A�A� �@��    BGǭ@��    @��    A�� �@��    BHE�@�!@    @�!@    A��2 �@���    BH��@�%     @�%     A�* �@� �    BIA�@�(�    @�(�    A�w �@�    BI��@�,�    @�,�    A��� �@�!C    BJ> @�0@    @�0@    A�� �@�8�    BJ�@�4     @�4     A��g �@�[    BK:@�7�    @�7�    A�� �@��j    BK�@�;�    @�;�    A�� �@�j�    BL6"@�?@    @�?@    A�*S �@˯�    BL�'@�C     @�C     A�� �@��s    BM2)@�F�    @�F�    A��- �@��N    BM�*@�J�    @�J�    A�� �@�X    BN.(@�N@    @�N@    A�
3 �@��`    BN�%@�R     @�R     A�� �@���    BO* @�U�    @�U�    A�b� �@��	    BO�@�Y�    @�Y�    A�pP �@�s    BP&@�]@    @�]@    A�:e �@���    BP�@�a     @�a     A��� �@��^    BQ!�@�d�    @�d�    A��/ �@��n    BQ��@�h�    @�h�    B p0 �@���    BR�@�l@    @�l@    A�j� �@�R    BR��@�p     @�p     B � �@��i    BS�@�s�    @�s�    A�g) �@�[B    BS��@�w�    @�w�    A�!z �@�w    BT~@�{@    @�{@    A�Í �@�H    BT�b@�     @�     A��� �@���    BUD@��    @��    A��� �@���    BU�$@熀    @熀    A��U �@�(    BV@�@    @�@    A��� �@��-    BV��@�     @�     A�� �@�8�    BW�@��    @��    A��� �@��e    BW��@畀    @畀    A��P �@�    BXb@�@    @�@    A�[� �@��H    BX�4@�     @�     B^� �@��    BY @��    @��    B 0� �@�.w    BY}�@礀    @礀    B�� �@�t�    BY��@�@    @�@    B�s �@�f    BZyd@�     @�     B5 �@� h    BZ�*@��    @��    B�G �@�    B[t�@糀    @糀    Bq� �@���    B[�@�@    @�@    B �@�y�    B\pm@�     @�     Br �@ۨF    B\�)@��    @��    B�� �@�
@    B]k�@�    @�    B�Z �@�7�    B]�@��@    @��@    B�� �@�=1    B^gL@��     @��     B	*R �@��c    B^��@���    @���    B�X �@��    B_b�@�р    @�р    B�� �@���    B_�W@��@    @��@    B�= �@���    B`]�@��     @��     B� �@���    B`ۥ@���    @���    B	r� �@�o�    BaYH@���    @���    BqA �@��2    Ba��@��@    @��@    BQ� �@��    BbT�@��     @��     B�X �@���    Bb�@���    @���    B� �@��    BcO�@��    @��    A��� �@�$�    Bc�J@��@    @��@    Bk� �@��J    BdJ�@��     @��     A�"� �@��    Bd�i@���    @���    A�P �@���    BeE�@���    @���    A��� �@�]�    Be�z@�@    @�@    B � �@��w    Bf@�@�     @�     BC� �@�G|    Bf�@�	�    @�	�    B�� �@���    Bg;�@��    @��    B�? �@�(G    Bg�w@�@    @�@    BB� �@��    Bh6�@�     @�     B�� �@��    Bh�a@��    @��    Bߚ �@��)    Bi1�@��    @��    B	�� �@��    Bi�=@� @    @� @    B�� �@�!>    Bj,�@�$     @�$     B
�� �@�    Bj�@�'�    @�'�    B�� �@�^    Bk'l@�+�    @�+�    B
Zg �@���    Bk��@�/@    @�/@    B
۸ �@�p�    Bl"$@�3     @�3     B�j �@�}�    Bl�z@�6�    @�6�    B>W �@�w_    Bm�@�:�    @�:�    BEp �@�WC    Bm�@�>@    @�>@    B�� �@�B    Bne@�B     @�B     B�u �@��    Bn��@�E�    @�E�    B�� �@�J�    Bo�@�I�    @�I�    B�� �@�}�    Bo�-@�M@    @�M@    B� �@���    Bpg@�Q     @�Q     B� �@��    Bp��@�T�    @�T�    B| �@���    Bq�@�X�    @�X�    B�] �@���    Bq��@�\@    @�\@    B5h �@��    Br%@�`     @�`     BU~ �@�l�    Br~J@�c�    @�c�    Bu� �@���    Br�j@�g�    @�g�    B4� �@���    Bsx�@�k@    @�k@    B~� �@�P�    Bs��@�o     @�o     B� �@���    Btr�@�r�    @�r�    B� �@��    Bt�@�v�    @�v�    Ba� �@�(�    Bul�@�z@    @�z@    B�� �@��y    Bu��@�~     @�~     B� �@�;    Bvf�@��    @��    BH� �@���    Bv��@腀    @腀    B( �@��    Bw`�@�@    @�@    B�� �@�B.    Bwݪ@�     @�     B?� �@�(�    BxZ�@��    @��    B�� �@�ch    Bx�z@蔀    @蔀    B�� �@���    ByT[@�@    @�@    B#�� �@�J\    By�6@�     @�     B$�� �@�F    BzN@��    @��    B$2$ �@��    Bz��@裀    @裀    B# �@ɍv    B{G�@�@    @�@    B)�J �@�#�    B{�k@�     @�     B) �@ɍR    B|A*@��    @��    B$4 �@�b�    B|��@貀    @貀    B$O �@͋�    B}:�@�@    @�@    B%�� �@���    B}�D@�     @�     B(�� �@د    B~3�@��    @��    B$� �@��    B~��@���    @���    B*�( �@�l_    B-'@��@    @��@    B-d �@�W�    B��@��     @��     B+�; �@��X    B�%@���    @���    B-�� �@�Ü    B�Qh@�Ѐ    @�Ѐ    B-� �@���    B���@��@    @��@    B,�� �A�    B���@��     @��     B*t� �A
�    B� @���    @���    B(� �AjP    B�JV@�߀    @�߀    B()a �AH�    B���@��@    @��@    B(�G �A�X    B�Ƹ@��     @��     B+ �A 0�    B��@���    @���    B'�� �A�e    B�C@��    @��    B.�h �A�p    B��/@��@    @��@    B1� �AV_    B��O@��     @��     B5�� �A 8�    B��l@���    @���    B2�& �A|m    B�;�@���    @���    B-;� �A    B�y�@�@    @�@    B)h �AWq    B���@�     @�     B's+ �AQ�    B���@��    @��    B%� �A/(    B�3�@��    @��    B!{� �A�    B�q�@�@    @�@    Bռ �A΃    B���@�     @�     B!�  �Alx�    B���@��    @��    B  �A$\�    B�+�@��    @��    Bd�  �A*1�    B�i�@�@    @�@    B�p  �A0s    B���@�#     @�#     B$�  �AF�    B��@�&�    @�&�    B"+�  �AM�?    B�#e@�*�    @�*�    B(5�  �AB��    B�aE@�.@    @�.@    B*p�  �A2��    B�� @�2     @�2     B+_�  �A"i�    B���@�5�    @�5�    B$9  �At+    B��@�9�    @�9�    B$2]  �AI�    B�X�@�=@    @�=@    B&.C �A~5    B��[@�A     @�A     B&X� �A�    B��@�D�    @�D�    B#� �A��    B��@�H�    @�H�    B �� �A�W    B�O�@�L@    @�L@    B"A �A"�    B��B@�P     @�P     B)z  �AR�    B���@�S�    @�S�    B�Z �A;��    B��@�W�    @�W�    B�� �AJ7�    B�F4@�[@    @�[@    B�� �A`i$    B���@�_     @�_     B�2 �An�    B��c@�b�    @�b�    B� �At�    B���@�f�    @�f�    B� �Aye`    B�<z@�j@    @�j@    Br �Alޙ    B�y�@�n     @�n     B%;W �AcT�    B��w@�q�    @�q�    B;�A  �A|~d    B���@�u�    @�u�    B?��  �Ao�:    B�2Z@�y@    @�y@    B �!  �A/bb    B�o�@�}     @�}     B,.  �A=�    B��!@��    @��    B6�Y  �AF`�    B��{@鄀    @鄀    B.H}  �AU    B�'�@�@    @�@    B' D �ANn�    B�e@�     @�     B#]� �Aa[�    B��\@��    @��    B$q� �Ak̄    B�ߘ@铀    @铀    B � �A}p�    B��@�@    @�@    B'��  �A+    B�Y�@�     @�     B>�%  �A�Y�    B��@��    @��    B@�  �A�X�    B��9@颀    @颀    B/��  �A���    B�M@�@    @�@    B � �A�    B�NY@�     @�     B*u9  �A���    B��[@��    @��    B,}L �A�NU    B��V@鱀    @鱀    B0� �A�    B�G@�@    @�@    B1�� �A��    B�B/@�     @�     B7� �A�X    B�@��    @��    B8-�  �A�/6    B���@���    @���    B>�  �A��    B���@��@    @��@    BB�� �A��_    B�5p@��     @��     BG�� �A��.    B�r'@���    @���    BK� �A���    B���@�π    @�π    BQJ� �A���    B��x@��@    @��@    BRm� �A�1�    B�(@��     @��     BU�� �A��[    B�d�@���    @���    BV� �A���    B��@�ހ    @�ހ    B\2� �Ao}�    B�ݖ@��@    @��@    Bc�1 �Ad��    B�@��     @��     Bd� �ATD�    B�V`@���    @���    Bl_2 �A:-�    B���@��    @��    BnO� �A6s    B���@��@    @��@    Br� �A+}    B�3@��     @��     Br��  �A    B�Ga@���    @���    Bpg  �A�G    B���@���    @���    Bp�  �A��    B���@� @    @� @    Bn�  �@�$�    B���@�     @�     BnO� �@ԴU    B�7�@��    @��    Bp' �@���    B�sv@��    @��    Bp� �@�.F    B��P@�@    @�@    By[- �@��    B��@�     @�     Bz$  �@���    B�&�@��    @��    Br-� �@�ƅ    B�b�@��    @��    Bn�| �@�Ζ    B��@�@    @�@    Bp�� �@ȱ�    B�٦@�"     @�"     Bi�� �@�2Q    B�@�%�    @�%�    Bl} �A �9    B�P�@�)�    @�)�    B��   �AJߔ    B���@�-@    @�-@    Bw`�  �An��    B��#@�1     @�1     BO�R  �A��    B�U@�4�    @�4�    BB�#  �A���    B�=t@�8�    @�8�    B,&�  �A���    B�x�@�<@    @�<@    B+�Z  �A�{A    B��x@�@     @�@     B(��  �A�    B��[@�C�    @�C�    B'ca  �A�(�    B�))@�G�    @�G�    B/�  �A��5    B�c�@�K@    @�K@    B5��  �A���    B���@�O     @�O     B:z"  �A�    B��@�R�    @�R�    B2�x  �A���    B��@�V�    @�V�    B6�  �A�    B�M�@�Z@    @�Z@    B='o  �A�S�    B��)