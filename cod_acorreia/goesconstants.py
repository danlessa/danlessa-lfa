#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

# These are satellite calibration constants for GOES-8, -10, -12, -13, and -14
# Usage in a separate code:
#
# from goesconstants import *
#

# VIS channel definitions
RadVISUnits = 'Radiance  [W/(m2 sr um)]'
RefVISUnits = 'Reflectance'
X0 = 29.0

mVIS8 = 0.5502
kVIS8 = 0.001062
C8 = 1.711  # Valid for GOES-10 from Jun2005-Dec2007

mVIS10 = 0.5582154
kVIS10 = 0.001110
C10 = 1.711  # Valid for GOES-10 from Jun2005-Dec2007

mVIS12 = 0.5771
kVIS12 = 0.001141
C12 = 1.0  # Valid for GOES-12 until Aug2007

mVIS13 = 0.6118208
kVIS13 = 0.001160
C13 = 1.255  # Valid for GOES-13 for FEB2013

mVIS14 = 0.5860814
kVIS14 = 0.00110636
C14 = 1.1171


# IR channel definitions
RadIRUnits = 'Radiance  [mW/(m2 sr cm-1)]'
TUnits = 'Brightness Temperature  [K]'
C1 = 1.191066e-5  # mW/(m2 sr cm-4)
C2 = 1.438833  # K/(cm-1)

# GOES 8 IR Channels
m28 = 227.3889
m38 = 38.8383
m48 = 5.2285
m68 = 5.0273
bb28 = 68.2167
bb38 = 29.1287
bb48 = 15.6854
bb68 = 15.3332
n28 = 2556.71
n38 = 1481.91
n48 = 934.30
n68 = 837.06
a28 = -0.618007
a38 = -0.656443
a48 = -0.519333
a68 = -0.383077
b28 = 1.001825
b38 = 1.001914
b48 = 1.002834
b68 = 1.000856
g28 = -6.021442e-7
g38 = -9.535221e-7
g48 = -3.005194e-6
g68 = 6.026892e-7
factor8 = 165.5 / 251.821294  # RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
# Solar constant averaged over GOES-8 Ch2 (W/m2) (Platnick and Fontenla, 2008)
F0_8 = 9.66116872950445
Bfunction8 = [8.00197E-6, 2.31711E-49, 1.99042E+1, -3.65106E-8, -1.93902E-10, 1.09415, -4.18331E-2, 6.63342E-4, -5.59083E-6, 2.64532E-8, -6.67303E-11, 7.02400E-14,
              3.41519E-7, 7.13808E-52, 5.86638E-4, 3.08368E-9, 1.41409E-11, 4.32652E-6, 3.26338E-8, 1.36474E-10, 5.48652E-13, 2.13147E-15, 7.51904E-18, 2.18539E-20]
satlat8 = 0.00  # Satellite subpoint latitude
satlon8 = -74.999954  # Satellite subpoint longitude


# GOES 10 IR Channels
m210 = 227.3889
m310 = 38.8383
m410 = 5.2285
m610 = 5.0273
bb210 = 68.2167
bb310 = 29.1287
bb410 = 15.6854
bb610 = 15.3332
n210 = 2552.9845
n310 = 1486.2212
n410 = 936.10260
n610 = 830.88473
a210 = -0.63343694
a310 = -0.66500842
a410 = -0.36939128
a610 = -0.32763317
b210 = 1.0013206
b310 = 1.0017857
b410 = 1.0017466
b610 = 1.0014057
g210 = -4.2038547e-7
g310 = -7.3885254e-7
g410 = -1.4981835e-6
g610 = -9.5563444e-7
# RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
factor10 = 164.6 / 252.267139805
# Solar constant averaged over GOES-10 Ch2 (W/m2) (Platnick and Fontenla, 2008)
F0_10 = 9.56249106144777
Bfunction10 = [1.23081E-04, 2.33024E-49, 1.99044E+01, -1.67954E-06, 5.61210E-09, 1.05730E+00, -4.06959E-02, 6.49225E-04, -5.50192E-06, 2.61627E-08, -6.62992E-11,
               7.00790E-14, 3.62089E-07, 8.98442E-52, 7.34503E-04, 3.69227E-09, 1.60863E-11, 8.18432E-03, 2.45021E-04, 3.01127E-06, 1.94516E-08, 6.96821E-11, 1.31320E-13, 1.01769E-16]
satlat10 = 0.00  # Satellite subpoint latitude
satlon10 = -59.999935  # Satellite subpoint longitude


# GOES 12 IR Channels
m212 = 227.3889
m312 = 38.8383
m412 = 5.2285
m612 = 5.5297
bb212 = 68.2167
bb312 = 29.1287
bb412 = 15.6854
bb612 = 16.5892
n212 = 2562.45
n312 = 1536.43
n412 = 933.21
n612 = 751.91
a212 = -0.727744
a312 = -5.278975
a412 = -0.534982
a612 = -0.177244
b212 = 1.002131
b312 = 1.016476
b412 = 1.002693
b612 = 1.000138
g212 = -1.173898e-6
g312 = -7.754348e-6
g412 = -2.667092e-6
g612 = 1.163496e-6
# RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
factor12 = 174.6 / 264.8402490292
# Solar constant averaged over GOES-12 Ch2 (W/m2) (Platnick and Fontenla, 2008)
F0_12 = 9.72421593651545
Bfunction12 = [1.15088E-004, 2.22766E-049, 1.99043E+001, -1.48400E-006, 4.67326E-009, 9.42467E-001, -3.75665E-002, 6.13815E-004, -5.28598E-006, 2.53975E-008, -6.47620E-011,
               6.86758E-014, 2.80977E-007, 5.85508E-052, 5.00534E-004, 2.54033E-009, 1.16359E-011, 4.28505E-006, 3.04522E-008, 1.26965E-010, 5.06002E-013, 1.95785E-015, 6.91205E-018, 2.01590E-020]
satlat12 = 0.00  # Satellite subpoint latitude
satlon12 = -74.999954  # Satellite subpoint longitude


# GOES 13 IR Channels
m213 = 227.3889
m313 = 38.8383
m413 = 5.2285
m613 = 5.5297
bb213 = 68.2167
bb313 = 29.1287
bb413 = 15.6854
bb613 = 16.5892
n213 = 2561.7421
n313 = 1522.5182
n413 = 937.23449
n613 = 749.82589
a213 = -1.4755462
a313 = -4.1556932
a413 = -0.52227011
a613 = -0.16089410
b213 = 1.0028656
b313 = 1.0142082
b413 = 1.0023802
b613 = 1.0006896
g213 = -5.8203946e-7
g313 = -8.0255086e-6
g413 = -2.0798856e-6
g613 = -3.9853774e-7
# RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
factor13 = 232.4 / 354.1531901968
# Solar constant averaged over GOES-13 Ch2 (W/m2) (Platnick and Fontenla, 2008)
F0_13 = 9.74384716655927
Bfunction13 = [3.22252E-004, 2.29538E-049, 1.99043E+001, -4.15215E-006, 1.32596E-008, 8.85164E-001, -3.57508E-002, 5.90154E-004, -5.12428E-006, 2.47901E-008, -6.35867E-011,
               6.77820E-014, 3.25625E-007, 5.56472E-052, 4.61646E-004, 2.56534E-009, 1.26165E-011, 4.35725E-006, 3.09653E-008, 1.29105E-010, 5.14528E-013, 1.99084E-015, 7.02852E-018, 2.04987E-020]
satlat13 = -0.125907  # Satellite subpoint latitude
satlon13 = -74.574043  # Satellite subpoint longitude


# GOES 14 IR Channels
m214 = 227.3889
m314 = 38.8383
m414 = 5.2285
m614 = 5.5297
bb214 = 68.2167
bb314 = 29.1287
bb414 = 15.6854
bb614 = 16.5892
n214 = 2577.3518
n314 = 1519.3488
n414 = 933.98541
n614 = 752.88143
a214 = -1.5565294
a314 = -3.9656363
a414 = -0.50944128
a614 = -0.16549136
b214 = 1.0027731
b314 = 1.0133248
b414 = 1.0029288
b614 = 1.0001953
g214 = -4.0683469e-7
g314 = -7.5834376e-6
g414 = -3.3213255e-6
g614 = 9.0998038e-7
# RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
factor14 = 244.1 / 367.8086681308
# Solar constant averaged over GOES-14 Ch2 (W/m2) (Platnick and Fontenla, 2008)
F0_14 = 9.89394846310663
Bfunction14 = [6.53738E-004, 2.13558E-049, 1.99043E+001, -8.20827E-006, 2.56581E-008, 9.69643E-001, -3.82718E-002, 6.20658E-004, -5.31370E-006, 2.54131E-008, -6.45656E-011,
               6.82709E-014, 2.69752E-007, 3.67356E-052, 3.27507E-004, 1.85944E-009, 9.82179E-012, 4.09287E-006, 2.90865E-008, 1.21271E-010, 4.83309E-013, 1.87004E-015, 6.60205E-018, 1.92549E-020]
satlat14 = -0.248697  # Satellite subpoint latitude
satlon14 = -105.926392  # Satellite subpoint longitude
