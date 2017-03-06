CDF     
      time          I   command_line      tsi_ingest -s mao -f M1    process_version       ingest-tsi-12.2-0.el6      dod_version       tsiskycover-b1-3.2     site_id       mao    facility_id       !M1: Manacapuru, Amazonia, Brazil       
data_level        b1     input_source      4/data/collection/mao/maotsiM1.00/20140822100800.dat    resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     sampling_interval         30 seconds     averaging_interval        None       tsi_name      TSI440/660     serial_number         102    band_hw       
30 pixels      band_top      -20 pixels     
brightness        7      cam_mir_dist      69 mm      camera_arm_width      
15 pixels      camera_head_height        100 pixels     camera_head_width         
40 pixels      ccd_height_factor         	0.012922       ccd_width_factor      	0.013000       center_x      240 pixels     center_y      325 pixels     image_display_height      427 pixels     image_display_width       320 pixels     image_height      640 pixels     image_small_height        160 pixels     image_small_width         120 pixels     image_width       480 pixels     max_proc_zen      80 degrees     orientation       north      reference_fitCoef0        	0.316565       reference_fitCoefMag0         214.187698     reference_fitCoef1        	0.000240       reference_fitCoefMag1         
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
datastream        maotsiskycoverM1.b1    history       Ycreated by user dsmgr on machine tin at 2014-08-22 12:49:03, using ingest-tsi-12.2-0.el6       ANDERS_input_file         V/http/ftp/armguest_user/bernardinelid1/184047/maotsiskycoverM1.b1.20140822.100800.cdf      ANDERS_processing_timestamp       Mon Mar  6 20:04:10 2017 UTC       ANDERS_armtime_timestamp      1488830650     ANDERS_version        ?ANDERS hg id:055c83de3fe4 rev:141+ Wed Sep 2 20:32:41 PDT 2015           	base_time                string        2014-08-22 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         �   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2014-08-22 00:00:00 0:00          �   time                	long_name         Time offset from midnight      units         'seconds since 2014-08-22 00:00:00 0:00          �   percent_thin                	long_name         Percent thin cloud     units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   sunny                   	long_name         Sunshine meter     units         	unitless       	valid_min                	valid_max               missing_value         ��     flag_values             flag_meanings         .sun_blocked_by_cloud sun_not_blocked_by_cloud      comment       70 = sun blocked by cloud, 1 = sun not blocked by cloud          �   percent_opaque                  	long_name         Percent opaque cloud       units         %      	valid_min                	valid_max         B�     
resolution        ?�     missing_value         �<         �   alt              	long_name         Altitude above mean sea level      units         m      standard_name         	altitude            �   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4     standard_name         
longitude           �   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�     standard_name         	latitude            �   qc_solar_altitude                   	long_name         ;Quality check results on field: Sun altitude above horizon     units         	unitless       description       7See global attributes for individual bit descriptions.          �   solar_altitude                  	long_name         Sun altitude above horizon     units         degree     	valid_min         ´     	valid_max         B�     
resolution        ?�     missing_value         �<         �S�� BH  �rdt�M�M@��     @��     ��  ����      <�@���    @���    ��  ����      >��@�׀    @�׀    ��  ����      >�*@��@    @��@    ��  ����      >���@��     @��     ��  ����      >�'m@���    @���    ��  ����      ?�9@��    @��    ��  ����      ?>$�@��@    @��@    ��  ����      ?]l�@��     @��     ��  ����      ?|��@���    @���    ��  ����      ?��@���    @���    ��  ����      ?���@��@    @��@    ��  ����      ?�Es@��     @��     ��  ����      ?���@� �    @� �    ��  ����      ?̌`@��    @��    ��  ����      ?�/�@�@    @�@    ��  ����      ?���@�     @�     ��  ����      ?�u�@��    @��    ��  ����      @�s@��    @��    ��  ����      @]�@�@    @�@    ��  ����      @/>@�     @�     ��  ����      @ �@��    @��    ��  ����      @$��@�"�    @�"�    ��  ����      @,� @�&@    @�&@    ��  ����      @4t#@�*     @�*     ��  ����      @<E6@�-�    @�-�    @���  �A<�%    @D:@�1�    @�1�    @�h;  �A<��    @K�/@�5@    @�5@    @���  �A<�    @S�@�9     @�9     @�A�  �A=Ud    @[��@�<�    @�<�    @�?f  �A>|�    @cY�@�@�    @�@�    @�z�  �AC t    @k*e@�D@    @�D@    @��:  �AEo@    @r�@�H     @�H     @��5  �AI9,    @zˠ@�K�    @�K�    @���  �AJxC    @�N@�O�    @�O�    @��   �AQj�    @�6M@�S@    @�S@    @��  �ASG    @��@�W     @�W     @�$5  �AVت    @��@�Z�    @�Z�    @�ˎ  �AX�>    @���@�^�    @�^�    @�3�  �A]*=    @���@�b@    @�b@    @ܕ�  �A_��    @���@�f     @�f     @�fD  �Aa��    @���@�i�    @�i�    @�|  �Ab��    @���@�m�    @�m�    @�Ԑ  �Ad�y    @�v�@�q@    @�q@    @���  �Adݩ    @�^�@�u     @�u     @�  �Ad�c    @�F�@�x�    @�x�    @��  �Af_    @�.�@�|�    @�|�    @��  �Af?�    @�w@�@    @�@    @�d�  �Ad�
    @��A@�     @�     @즺  �Ae�    @��@��    @��    @���  �Ae5�    @�ͻ@⋀    @⋀    @�  �Abk�    @õj@�@    @�@    @  �Ag(�    @ǝ@�     @�     @��{  �Ah;�    @˄�@��    @��    @���  �Agä    @�l@@⚀    @⚀    @�$S  �Ai8p    @�S�@�@    @�@    @�M  �AiKe    @�;J@�     @�     @��  �AjD�    @�"�@��    @��    @�ip  �Am�1    @�
.@⩀    @⩀    @���  �Aqբ    @��@�@    @�@    A �  �AqJ�    @���@�     @�     A??  �Am�K    @��<@��    @��    A�  �Anp�    @@⸀    @⸀    A�  �AqM�    @�@�@    @�@    Ad�  �At=�    @�u�@��     @��     A�  �Am�    @�]@���    @���    A�	  �At    @�D5@�ǀ    @�ǀ    A�m  �Av�    A�@��@    @��@    A�  �Aum    A	)@��     @��     A��  �As�.    A��@���    @���    A��  �Av��    A�#@�ր    @�ր    A
�}  �Aw�    A�@��@    @��@    A	��  �Aw�    A
�@��     @��     A
py  �Av��    A�r@���    @���    A@�  �As�D    A��@��    @��    AŬ  �Av/�    A�6@��@    @��@    A:�  �AuU�    A��@��     @��     Ay�  �Asu�    A��@���    @���    A@�  �AtO�    A�4@��    @��    A��  �Ar��    A~~@��@    @��@    A؇  �Asq    Aq�@��     @��     A	��  �At�     Ae@���    @���    A�  �Av �    AX:@��    @��    A�%  �Av�.    A Km@�@    @�@    A�  �AsFO    A">�@�     @�     A��  �Av     A$1�@��    @��    A]  �At    A&$�@��    @��    As!  �As�2    A( @�@    @�@    A�
  �Ak    A*@�     @�     A�S  �Ap#�    A+�'@��    @��    A�+  �Ao�,    A-�1@�!�    @�!�    A�  �Aj`}    A/�6@�%@    @�%@    AV~  �Ak�_    A1�4@�)     @�)     A+G  �Ak�    A3�,@�,�    @�,�    AP2  �AiJ    A5�@�0�    @�0�    AI�  �Ai�y    A7�@�4@    @�4@    A�Z  �Aio�    A9��@�8     @�8     A��  �Ag��    A;��@�;�    @�;�    AS  �AjI�    A=��@�?�    @�?�    A�?  �Ad��    A?{}@�C@    @�C@    A�  �Af�     AAnJ@�G     @�G     A��  �Ah    ACa@�J�    @�J�    Aĭ  �Ad#    AES�@�N�    @�N�    AO  �Ae6    AGF�@�R@    @�R@    Aʡ  �Ab\�    AI9>@�V     @�V     AE�  �Ah�    AK+�@�Y�    @�Y�    A?�  �Af�M    AM�@�]�    @�]�    A�o  �Ac_�    AO0@�a@    @�a@    A��  �Adl`    AQ�@�e     @�e     A�z  �A`�    AR�[@�h�    @�h�    A�  �AcRF    AT��@�l�    @�l�    A��  �Ab��    AV�k@�p@    @�p@    Ah  �AQ��    AX��@�t     @�t     A��  �A^c�    AZ�`@�w�    @�w�    A�|  �A]��    A\��@�{�    @�{�    A]�  �A^p�    A^�9@�@    @�@    A�<  �A\g;    A`��@�     @�     A+'  �A[��    Ab��@��    @��    A  �AY�T    Ad|K@㊀    @㊀    A��  �AY>o    Afn�@�@    @�@    A��  �AZgj    Ah`�@�     @�     A�  �AV��    AjS@��    @��    AK  �AW"�    AlES@㙀    @㙀    A��  �AZ��    An7�@�@    @�@    AW�  �AYX#    Ap)�@�     @�     A�  �A\�?    Ar�@��    @��    A��  �AV�k    At�@㨀    @㨀    A>  �AU�    Au��@�@    @�@    A�  �AWCj    Aw�@�     @�     Aq  �AS�    Ay�@��    @��    A0<  �ASp    A{�@㷀    @㷀    A�   �AR��    A}��@�@    @�@    A��  �AN��    A��@�     @�     A��  �AP�(    A���@���    @���    A�  �AR�    A���@�ƀ    @�ƀ    A&  �AK6T    A�ǽ@��@    @��@    A ��  �AN��    A���@��     @��     A�5  �AP7�    A���@���    @���    A �	  �AL��    A��d@�Հ    @�Հ    A$�J  �AL�c    A��?@��@    @��@    A#�  �AJ�9    A��@��     @��     A&��  �AI1W    A���@���    @���    A%RJ  �AG$�    A���@��    @��    A&k�  �AE��    A���@��@    @��@    A%��  �AC(    A��F@��     @��     A&B~  �AC}�    A��@���    @���    A'�H  �AB�    A�x�@��    @��    A(�  �AA�f    A�q@��@    @��@    A&��  �A?��    A�j4@��     @��     A$�  �AA!�    A�b�@���    @���    A&u  �A?m�    A�[�@��    @��    A(u  �A>N1    A�T:@�@    @�@    A'W&  �A=y�    A�L�@�
     @�
     A'W  �A=�    A�E~@��    @��    A&Bt  �A?�	    A�>@��    @��    A%�  �A<�    A�6�@�@    @�@    A(�8  �A:��    A�/B@�     @�     A&,U  �A;�%    A�'�@��    @��    A#�^  �A9�    A� Y@� �    @� �    A(�d  �A6g�    A��@�$@    @�$@    A+3�  �A7+{    A�^@�(     @�(     A(�T  �A4�^    A�	�@�+�    @�+�    A&]�  �A-q    A�P@�/�    @�/�    A#��  �A(�8    A���@�3@    @�3@    A!��  �A(�    A��0@�7     @�7     A!M  �A)�\    A��@�:�    @�:�    A�D  �A"�    A���@�>�    @�>�    A�C  �A#P}    A��\@�B@    @�B@    A	  �A"V�    A�Է@�F     @�F     A�Z  �A}�    A��@�I�    @�I�    Ag�  �Aq.    A��]@�M�    @�M�    A!�R  �A�-    A���@�Q@    @�Q@    A��  �A1d    A���@�U     @�U     A"͢  �AI    A��3@�X�    @�X�    A�q  �A�3    A��p@�\�    @�\�    A�  �A �w    A���@�`@    @�`@    A�p  �A M`    A���@�d     @�d     Ag�  �A G    A��	@�g�    @�g�    A�  �A �    A��2@�k�    @�k�    A��  �A!˄    A�V@�o@    @�o@    A�  �A ��    A�wt@�s     @�s     A
3  �A �    A�o�@�v�    @�v�    A3  �A    A�g�@�z�    @�z�    A��  �A�     A�_�@�~@    @�~@    A3�  �Ad�    A�W�@�     @�     Asz  �Ak�    A�O�@��    @��    A��  �A"Sh    A�G�@䉀    @䉀    A�  �A!S�    A�?�@�@    @�@    A�  �A"�j    A�7�@�     @�     A��  �A"}�    A�/�@��    @��    A��  �A#�9    A�'�@䘀    @䘀    A��  �A��    A�h@�@    @�@    A�<  �A#�    A�F@�     @�     AU  �A%h�    A�@��    @��    A�  �A#\�    A��@䧀    @䧀    A"   �A!�y    A���@�@    @�@    A ��  �A ��    A���@�     @�     A!��  �A �V    A��I@��    @��    A#8�  �A!�    A��@䶀    @䶀    A"_  �A�    A�ݻ@�@    @�@    A!ID  �A�=    A��l@�     @�     A#0�  �A��    A��@���    @���    A��  �Aj�    A�Ļ@�ŀ    @�ŀ    A��  �A��    AļZ@��@    @��@    A��  �As5    Aų�@��     @��     A�  �A�l    Aƫ�@���    @���    A!c  �A�    Aǣ@�Ԁ    @�Ԁ    A�!  �Aa�    AȚ�@��@    @��@    A�N  �A�s    Aɒ@��     @��     AA�  �A��    Aʉ�@���    @���    AU  �A~    Aˁ@��    @��    AD  �A	�    A�xv@��@    @��@    A�  �A�s    A�o�@��     @��     A^  �A"    A�g?@���    @���    A��  �At�    A�^�@��    @��    A-z  �A�    A�U�@��@    @��@    A�R  �A)o    A�M=@��     @��     A�0  �A�    A�D�@���    @���    A$G  �Ad�    A�;�@��    @��    A�v  �A�o    A�3 @�@    @�@    A��  �AD�    A�*4@�	     @�	     A�%  �AW�    A�!b@��    @��    ALj  �Aٶ    A��@��    @��    A2`  �A�A    A��@�@    @�@    A  �ASR    A��@�     @�     AC�  �A��    A���@��    @��    AH�  �A	w    A���@��    @��    A[  �Ai    A���@�#@    @�#@    A��  �Aq�    A���@�'     @�'     Aq�  �A��    A���@�*�    @�*�    At�  �A�    A���@�.�    @�.�    A  �A	�    A�ǰ@�2@    @�2@    A�  �A�    Aྑ@�6     @�6     A��  �AJl    A�l@�9�    @�9�    A؄  �Az    A�?@�=�    @�=�    A ��  �A��    A�@�A@    @�A@    A �  �A�G    A��@�E     @�E     A��  �A �    A同@�H�    @�H�    A��  �A!��    A�B@�L�    @�L�    A
�  �A"�A    A�}�@�P@    @�P@    A
�  �A%{�    A�t�@�T     @�T     A�  �A'y�    A�k6@�W�    @�W�    @��%  �@��v    A�a�@�[�    @�[�    A?�  �A%dK    A�X^@�_@    @�_@    A	�  �A&ND    A�N�@�c     @�c     A
�  �A*l�    A�Eg@�f�    @�f�    A	�  �A*:*    A�;�@�j�    @�j�    A�
  �A*��    A�2P@�n@    @�n@    AU�  �A)��    A�(�@�r     @�r     A�  �A'F    A�@�u�    @�u�    A6>  �A%Ȱ    A�t@�y�    @�y�    A"�  �A$��    A��@�}@    @�}@    A%��  �A#��    A�@�     @�     A'2w  �A�y    A��O@��    @��    A/	�  �A!73    A��@刀    @刀    A6�#  �A#�    A��@�@    @�@    A<��  �A"�;    A���@�     @�     A@��  �A&R�    A�� @��    @��    AA_�  �A'�|    A��@嗀    @嗀    ABE�  �A&�'    A��&@�@    @�@    AC�i  �A*�    A��,@�     @�     AAh�  �A-!]    A��*@��    @��    AA�  �A0��    A��@妀    @妀    A?�  �A3�     A��@�@    @�@    A9Y�  �A7�r    A���@�     @�     A<��  �A:n�    B @e@��    @��    A:��  �A?r<    B �N@嵀    @嵀    A<�/  �AA'�    B63@�@    @�@    A7aX  �AA#�    B�@�     @�     A;M�  �A?r�    B+�@���    @���    A:jJ  �AA�<    B��@�Ā    @�Ā    A=��  �AF�    B!�@��@    @��@    A>\�  �AD��    B�e@��     @��     AG��  �AE`    B/@���    @���    AK�  �AC4    B��@�Ӏ    @�Ӏ    AN��  �ACЭ    B�@��@    @��@    AP;�  �A;�    B�n@��     @��     AQjs  �A=H    B$@���    @���    AZsW  �A< c    B|�@��    @��    A^ �  �A7�    B��@��@    @��@    A`  �A5��    Br*@��     @��     Ak"�  �A `    B��@���    @���    An��  �A ��    Bgk@��    @��    Aq�[  �A!��    B�@��@    @��@    Av�U  �A"t;    B	\�@��     @��     A~S�  �A 7K    B	�&@���    @���    A~2  �A�    B
Q�@� �    @� �    A}N�  �A2B    B
�5@�@    @�@    A���  �AD�    BF�@�     @�     A�d�  �A"�    B�/@��    @��    A�+�  �A�K    B;�@��    @��    A��  �A�f    B�@�@    @�@    A���  �A�<    B0@�     @�     A�2  �A��    B��@��    @��    A���  �A�p    B%E@��    @��    A�5)  �A�M    B��@�"@    @�"@    A���  �A��    B�@�&     @�&     A�*  �A2    B�D@�)�    @�)�    A���  �A(�    B�@�-�    @�-�    A��Y  �Aw�    B��@�1@    @�1@    A�
�  �A��    B@�5     @�5     A�5,  �A^    B}L@�8�    @�8�    A���  �AA    B��@�<�    @�<�    A���  �AG    Bq�@�@@    @�@@    A�5�  �A�    B��@�D     @�D     A��'  �A�    Be�@�G�    @�G�    A���  �A��    B�@�K�    @�K�    A�1�  �A    BZ-@�O@    @�O@    A�0�  �A�T    B�>@�S     @�S     A~Gd  �A�    BNI@�V�    @�V�    A|=�  �A"qC    B�O@�Z�    @�Z�    AzH  �A~    BBN@�^@    @�^@    Au�{  �A�J    B�G@�b     @�b     Ax��  �AN�    B6:@�e�    @�e�    Av�2  �A    B�'@�i�    @�i�    A}�  �Ag4    B*@�m@    @�m@    A~��  �A�    B��@�q     @�q     A�|�  �AK�    B�@�t�    @�t�    A�A�  �AE�    B��@�x�    @�x�    A��e  �AO9    Bi@�|@    @�|@    A���  �Ab0    B�0@�     @�     A��  �A�6    B�@��    @��    A�Ɠ  �A�z    B~�@懀    @懀    A��K  �A�T    B�^@�@    @�@    A�5  �A�    Br@�     @�     A��\  �A�O    B�@��    @��    A���  �A      BeQ@斀    @斀    A�_�  �A    B��@�@    @�@    A�x�  �Aɹ    BX|@�     @�     A� _  �At;    B�@��    @��    A�l  �A��    BK�@楀    @楀    A�~8  �A�     B�@�@    @�@    A�!  �A��    B >@�     @�     Af�  �A�    B ��@��    @��    A~:'  �A��    B!1V@洀    @洀    A~�  �A۸    B!��@�@    @�@    A��d  �A��    B"$@�     @�     A�;�  �A ��    B"�c@��    @��    A�V  �@�K_    B#�@�À    @�À    A~��  �@��#    B#��@��@    @��@    A�b�  �A7    B$	.@��     @��     A~��  �A�S    B$�c@���    @���    A�  �A��    B$��@�Ҁ    @�Ҁ    A�)�  �A7D    B%t�@��@    @��@    A�،  �@�@m    B%��@��     @��     A��  �AO    B&f�@���    @���    A�(W  �@�./    B&��@��    @��    A�"^  �@�%    B'Y @��@    @��@    A�[d  �@�EJ    B'��@��     @��     A�*�  �@�v    B(J�@���    @���    A���  �@�"    B(��@���    @���    A�\�  �@�-�    B)<�@��@    @��@    A���  �@� �    B)��@��     @��     A�{�  �@���    B*.�@���    @���    A�:  �@���    B*�O@���    @���    A��  �@�ŵ    B+ @�@    @�@    A�tC  �@�V    B+��@�     @�     A���  �@��2    B,�@�
�    @�
�    A�k�  �@鮥    B,�5@��    @��    A���  �@��    B-�@�@    @�@    A�&_  �@��    B-{t@�     @�     A���  �@�'    B-�@��    @��    A�{  �@��v    B.l�@��    @��    A�v  �@�R    B.�@�!@    @�!@    A�@�  �@��\    B/]�@�%     @�%     A�w�  �@��f    B/��@�(�    @�(�    A�h4  �@��R    B0N`@�,�    @�,�    A�H9  �@�S�    B0ƽ@�0@    @�0@    A�o�  �@���    B1?@�4     @�4     A���  �@��A    B1�[@�7�    @�7�    A�}�  �@��    B2/�@�;�    @�;�    A��g  �@��    B2��@�?@    @�?@    A�7  �@��3    B3 @�C     @�C     A��  �@�rc    B3�'@�F�    @�F�    A���  �@�Q�    B4C@�J�    @�J�    A�+�  �@���    B4�T@�N@    @�N@    A��  �@��    B5 \@�R     @�R     A�tN  �@��    B5xY@�U�    @�U�    A���  �@�4�    B5�M@�Y�    @�Y�    A���  �@���    B6h7@�]@    @�]@    A���  �@�p�    B6�@�a     @�a     A��O  �@�dn    B7W�@�d�    @�d�    A�Kx  �@��    B7Ϸ@�h�    @�h�    A�0�  �@���    B8Gx@�l@    @�l@    A���  �@�e    B8�.@�p     @�p     A���  �@�xK    B96�@�s�    @�s�    A�+  �@���    B9�{@�w�    @�w�    A��	  �@�o'    B:&@�{@    @�{@    A��  �@�DV    B:��@�     @�     A�tf  �@��    B;@��    @��    A�v@  �@��6    B;��@熀    @熀    A�S  �@���    B<�@�@    @�@    A���  �@�b�    B<{^@�     @�     A�W�  �@�f�    B<�@��    @��    A��<  �@���    B=i�@畀    @畀    A�-&  �@���    B=�9@�@    @�@    A���  �@��x    B>Xk@�     @�     A�4�  �@��e    B>ϒ@��    @��    A���  �@x�    B?F�@礀    @礀    A�U�  �@jPD    B?��@�@    @�@    A�/~  �@g��    B@4�@�     @�     A���  �@_]    B@��@��    @��    A��<  �@T�$    BA"�@糀    @糀    A���  �@U A    BA�@�@    @�@    A�[  �@\�1    BBQ@�     @�     A�z%  �@KF�    BB�@��    @��    A���  �@E�Z    BB��@�    @�    A���  �@J��    BCt}@��@    @��@    A��D  �@P��    BC�@��     @��     A���  �@>��    BDa�@���    @���    A��  �@H\Q    BD�6@�р    @�р    A�<A  �@F�!    BEN�@��@    @��@    A�Ͱ  �@L�    BE�@��     @��     Aɚc  �@C^�    BF;y@���    @���    A�7  �@L�    BF��@���    @���    A��_  �@U�Z    BG(@��@    @��@    A�/  �@Y�-    BG�C@��     @��     A�d�  �@[�    BHk@���    @���    Aǂ�  �@c�x    BH��@��    @��    A���  �@Y.�    BI �@��@    @��@    A�H�  �@Wg�    BIv�@��     @��     AƸ�  �@X��    BI�~@���    @���    A�Ӏ  �@Z    BJb_@���    @���    A�N�  �@\"�    BJ�2@�@    @�@    A���  �@e�j    BKM�@�     @�     Aƻ�  �@w��    BKë@�	�    @�	�    A�l�  �@�W�    BL9Q@��    @��    A��  �@�c1    BL��@�@    @�@    A�(�  �@���    BM$q@�     @�     A��"  �@�2�    BM��@��    @��    A���  �@�q    BNR@��    @��    A�:}  �@���    BN��@� @    @� @    A��z  �@�:�    BN��@�$     @�$     A�3@  �@��    BOo0@�'�    @�'�    A�`�  �@�7     BO�Z@�+�    @�+�    A���  �@|�x    BPYs@�/@    @�/@    A��  �@r8i    BP�}@�3     @�3     A�A�  �@Z��    BQCu@�6�    @�6�    A��  �@T�<    BQ�^@�:�    @�:�    A�d�  �@X�    BR-5@�>@    @�>@    A�^  �@�?�    BR��@�B     @�B     A��3  �@���    BS�@�E�    @�E�    Aۺ�  �@���    BS�U@�I�    @�I�    A��$  �@{��    BS��@�M@    @�M@    A�4'  �@~�C    BTti@�Q     @�Q     Aκ�  �@��k    BT��@�T�    @�T�    A��5  �@��F    BU]5@�X�    @�X�    AʣB  �@�(q    BUс@�\@    @�\@    A���  �@Ҫ�    BVE�@�`     @�`     Aɥ`  �@�G�    BV��@�c�    @�c�    A�9�  �ADI    BW-�@�g�    @�g�    A�M  �A��    BW��@�k@    @�k@    A�1�  �A7F�    BX�@�o     @�o     A�Z  �AU��    BX��@�r�    @�r�    A�(�  �Ac�'    BX��@�v�    @�v�    A�-  �AbMB    BYq<@�z@    @�z@    A�l  �A��    BY��@�~     @�~     A�e  �A��t    BZXi@��    @��    A���  �A��    BZ��@腀    @腀    A�t�  �A�)�    B[?F@�@    @�@    A�0�  �A��N    B[��@�     @�     A��#  �A�c�    B\%�@��    @��    BË  �A�    B\��@蔀    @蔀    B��  �A�u�    B]	@�@    @�@    B�  �A��6    B]@�     @�     B�l  �A�#�    B]��@��    @��    B�#  �A�J�    B^d�@裀    @裀    BXk  �A��    B^�x@�@    @�@    A��  �B
m�    B_J@�     @�     A��  �B    B_��@��    @��    A��  �B"�    B`/#@貀    @貀    Aӯ�  �B(�Z    B`��@�@    @�@    A�Fr  �B3w�    Ba�@�     @�     AƠ�  �B?�    Ba�@��    @��    A�Ʀ  �BK�P    Ba�@���    @���    A��  �BPʪ    Bbj @��@    @��@    A�,�  �B]b    Bb�@��     @��     A��M  �B`f|    BcM�@���    @���    Bi�  �BE��    Bc��@�Ѐ    @�Ѐ    B
H  �BYi�    Bd1=@��@    @��@    B/�[  �B-    Bd��@��     @��     B%��  �BB&�    Be7@���    @���    B9��  �B4ٙ    Be��@�߀    @�߀    BlD�  �B 9�    Be��@��@    @��@    Bf3  �B�+    Bfg�@��     @��     B{r�  �A���    Bf��@���    @���    Bn��  �Bcm    BgI�@��    @��    B�W  �A��    Bg��@��@    @��@    B��#  �A�x?    Bh+q@��     @��     B��  �A�Ǎ    Bh�@���    @���    B��e  �A�}5    Bi�@���    @���    B��  �A��    Bi|�@�@    @�@    B��I  �A�     Bi�>@�     @�     B�[�  �A��6    Bj]j@��    @��    B�  �A�7�    Bj�z@��    @��    B��?  �A�}�    Bk=m@�@    @�@    B�1�  �A�yY    Bk�B@�     @�     B�CF  �Ax�T    Bl�@��    @��    B�Κ  �A���    Bl��@��    @��    Bw��  �A��p    Bl�@�@    @�@    Bm��  �A�KU    Bmkk@�#     @�#     Bmh�  �A�'�    Bmڨ@�&�    @�&�    Bg�6  �A�(�    BnI�@�*�    @�*�    Be��  �A�2    Bn��@�.@    @�.@    Bf�  �A�c    Bo'�@�2     @�2     Bcxg  �A��f    Bo�a@�5�    @�5�    B]�u  �A���    Bp�@�9�    @�9�    B\/�  �A���    Bpsz@�=@    @�=@    B[�U  �A�cq    Bp��@�A     @�A     B\:�  �A���    BqP@�D�    @�D�    B\7�  �A�6    Bq�#@�H�    @�H�    BV��  �A�X\    Br,@�L@    @�L@    BW?�  �A�0    Br��@�P     @�P     BY��  �A�O�    Bs�@�S�    @�S�    B_fk  �Aܰ�    Bsu