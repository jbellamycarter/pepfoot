#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright 2019 Jeddidiah Bellamy-Carter, University of Nottingham
#
#    This file is part of PepFoot.
#
#    PepFoot is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    PepFoot is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with PepFoot.  If not, see <https://www.gnu.org/licenses/>.


import os
import re
import numpy as np
from scipy import signal
from datetime import datetime
from pyteomics import mass
from collections import OrderedDict


# Default settings values

default_enzymes = {
        'ArgC':'[R]',
        'AspN':'(?=[D])',
        'AspN + N-term Glu':'(?=[DE])',
        'BNPS-Skatole':'[W]',
        'Chymotrypsin':'[FLMYW](?=[^P])',
        'Chymotrypsin High Specificity':'[FYW](?=[^P])',
        'Clostripain':'[R]',
        'Cyanogen bromide':'[M]',
        'Elastase':'[AVSLI]',
        'Formic acid':'[D]',
        'GluC':'[E]',
        'Glutamyl endopeptidase':'[E]',
        'Hydroxylamine':'[N](?=[G])',
        'Iodosobenzoic acid':'[W]',
        'LysC':'[K]',
        'LysN':'(?=[K])',
        'NTCB + Ni':'(?=[C])',
        'Neutrophil elastase':'[VA]',
        'Pepsin':'[FL](?=[^AGV])',
        'Proline endopeptidase':'[P](?=[^P])',
        'Proteinase K':'[AEFILTVWY]',
        'Thermolysin':'[^DE](?=[AFILMV])',
        'Trypsin':'[KR](?=[^P])'
        }

default_mods = {
        'Acetyl': {
            'ID': 'ac',
            'Gain': 'C2H3O',
            'Loss': 'H',
            'Mass': 42.0106,
            'Residues': 'K'
        },
        'Acrylamide adduct': {
            'ID': 'acr',
            'Gain': 'C3H6ON',
            'Loss': 'H',
            'Mass': 71.0371,
            'Residues': 'C'
        },
        'Amidation': {
            'ID': 'am',
            'Gain': 'NH2',
            'Loss': 'OH',
            'Mass': -0.984,
            'Residues': 'DE'
        },
        'Biotinylation': {
            'ID': 'bt',
            'Gain': 'C10H15N2O2S',
            'Loss': 'H',
            'Mass': 226.0776,
            'Residues': 'K'
        },
        'Carbamylation': {
            'ID': 'ca',
            'Gain': 'H2CNO',
            'Loss': 'H',
            'Mass': 43.0058,
            'Residues': 'KRCM'
        },
        'Carbamidomethyl': {
            'ID': 'cam',
            'Gain': 'C2H4NO',
            'Loss': 'H',
            'Mass': 57.0215,
            'Residues': 'C'
        },
        'Deamidation': {
            'ID': 'deam',
            'Gain': 'OH',
            'Loss': 'NH2',
            'Mass': 0.984,
            'Residues': 'NQ'
        },
        'Nitrosylation': {
            'ID': 'n',
            'Gain': 'NO2',
            'Loss': 'H',
            'Mass': 44.9851,
            'Residues': 'WY'
        },
        'Oxidation': {
            'ID': 'ox',
            'Gain': 'O',
            'Loss': '',
            'Mass': 15.994915,
            'Residues': 'ACDEFGHIKLMNPQRSTVWY'
        },
        'Oxidation of Met': {
            'ID': 'oxm',
            'Gain': 'O',
            'Loss': '',
            'Mass': 15.994915,
            'Residues': 'M'
        },
        'Phosphate': {
            'ID': 'p',
            'Gain': 'H2PO3',
            'Loss': 'H',
            'Mass': 79.9663,
            'Residues': 'STY'
        },
        'Photoleucine': {
            'ID': 'pleu',
            'Gain': 'C5H9NO2',
            'Loss': '',
            'Mass': 115.0633,
            'Residues': 'ACDEFGHIKLMNPQRSTVWY'
        },
        'Sulfation': {
            'ID': 's',
            'Gain': 'HO3S',
            'Loss': 'H',
            'Mass': 79.9568,
            'Residues': 'STY'
        },
        'Aryldiazarine-TDBA': {
            'ID': 'tdba',
            'Gain': 'C9H5F3O2',
            'Loss': '',
            'Mass': 202.0242,
            'Residues': 'ACDEFGHIKLMNPQRSTVWY'
        }
}


# Dictionary containing masses and relative abundances of elements.
# Adapted from pyteomics mass.nist_mass to remove zero abundances
# and only relative to most abundant isotope
_rel_mass = {
        'Ac': {227: (227, 1.0)},
        'Ag': {107: (106.905097, 1.0), 109: (108.904752, 0.9290495572831265)},
        'Al': {27: (26.98153863, 1.0)},
        'Am': {243: (243, 1.0)},
        'Ar': {36: (35.967545106, 0.0033785038800083936), 38: (37.9627324, 0.0006345362413567027), 40: (39.9623831225, 1.0)},
        'As': {75: (74.9215965, 1.0)},
        'At': {210: (210, 1.0)},
        'Au': {197: (196.9665687, 1.0)},
        'B': {10: (10.012937, 0.2484394506866417), 11: (11.0093054, 1.0)},
        'Ba': {130: (129.9063208, 0.0014784233869842953), 132: (131.9050613, 0.001408686434768055), 134: (133.9045084, 0.03371084270133058), 135: (134.9056886, 0.09194119780189128), 136: (135.9045759, 0.10954280454127033), 137: (136.9058274, 0.15665708945856233), 138: (137.9052472, 1.0)},
        'Be': {9: (9.0121822, 1.0)},
        'Bh': {272: (272, 1.0)},
        'Bi': {209: (208.9803987, 1.0)},
        'Bk': {247: (247, 1.0)},
        'Br': {79: (78.9183371, 1.0), 81: (80.9162906, 0.9727756954034326)},
        'C': {12: (12.0, 1.0), 13: (13.0033548378, 0.010815728292732234)},
        'Ca': {40: (39.96259098, 1.0), 42: (41.95861801, 0.006674162635004797), 43: (42.9587666, 0.001392599622450769), 44: (43.9554818, 0.021518243055054107), 46: (45.9536926, 4.126221103557835e-05), 48: (47.952534, 0.0019290083659132874)},
        'Cd': {106: (105.906459, 0.0435085276714236), 108: (107.904184, 0.030978071702053602), 110: (109.9030021, 0.43473720849286457), 111: (110.9041781, 0.44552732335537765), 112: (111.9027578, 0.8398886181691612), 113: (112.9044017, 0.4253393665158371), 114: (113.9033585, 1.0), 116: (115.904756, 0.26070309780717016)},
        'Ce': {136: (135.907172, 0.002091577162238553), 138: (137.905991, 0.00283776144714528), 140: (139.9054387, 1.0), 142: (141.909244, 0.1256529112492934)},
        'Cf': {251: (251, 1.0)},
        'Cl': {35: (34.96885268, 1.0), 37: (36.96590259, 0.3199577613516367)},
        'Cm': {247: (247, 1.0)},
        'Cn': {285: (285, 1.0)},
        'Co': {59: (58.933195, 1.0)},
        'Cr': {50: (49.9460442, 0.05185644893721133), 52: (51.9405075, 1.0), 53: (52.9406494, 0.11339197269331296), 54: (53.9388804, 0.02822566207974794)},
        'Cs': {133: (132.905451933, 1.0)},
        'Cu': {63: (62.9295975, 1.0), 65: (64.9277895, 0.44613159797541574)},
        'Db': {268: (268, 1.0)},
        'Ds': {281: (281, 1.0)},
        'Dy': {156: (155.924283, 0.001981599433828733), 158: (157.924409, 0.003361641896673744), 160: (159.9251975, 0.08241330502477), 161: (160.9269334, 0.6684005661712668), 162: (161.9267984, 0.9014508138711959), 163: (162.9287312, 0.8809624911535738), 164: (163.9291748, 1.0)},
        'Er': {162: (161.928778, 0.004148882189654658), 164: (163.9292, 0.047786765364295734), 166: (165.9302931, 1.0), 167: (166.9320482, 0.6825955884547653), 168: (167.9323702, 0.8052413216726861), 170: (169.9354643, 0.44503477300540256)},
        'Es': {252: (252, 1.0)},
        'Eu': {151: (150.9198502, 0.9160758766047136), 153: (152.9212303, 1.0)},
        'F': {19: (18.99840322, 1.0)},
        'Fe': {54: (53.9396105, 0.06370294483074307), 56: (55.9349375, 1.0), 57: (56.935394, 0.02309436100878436), 58: (57.9332756, 0.0030734354905508207)},
        'Fm': {257: (257, 1.0)},
        'Fr': {223: (223, 1.0)},
        'Ga': {69: (68.9255736, 1.0), 71: (70.9247013, 0.6636720569641312)},
        'Gd': {152: (151.919791, 0.008051529790660225), 154: (153.9208656, 0.08776167471819646), 155: (154.922622, 0.5958132045088567), 156: (155.9221227, 0.824074074074074), 157: (156.9239601, 0.6300322061191627), 158: (157.9241039, 1.0), 160: (159.9270541, 0.8800322061191626)},
        'Ge': {70: (69.9242474, 0.5550108932461874), 72: (71.9220758, 0.7437363834422658), 73: (72.9234589, 0.21132897603485837), 74: (73.9211778, 1.0), 76: (75.9214026, 0.21323529411764702)},
        'H': {1: (1.00782503207, 1.0), 2: (2.0141017778, 0.00011501322652104991)},
        'He': {3: (3.0160293191, 1.3400017956024062e-06), 4: (4.00260325415, 1.0)},
        'Hf': {174: (173.940046, 0.004561003420752566), 176: (175.9414086, 0.1499429874572406), 177: (176.9432207, 0.5302166476624858), 178: (177.9436988, 0.7776510832383123), 179: (178.9458161, 0.3882554161915621), 180: (179.94655, 1.0)},
        'Hg': {196: (195.965833, 0.0050234427327528475), 198: (197.966769, 0.33389149363697257), 199: (198.9682799, 0.5649698593436034), 200: (199.968326, 0.7736101808439385), 201: (200.9703023, 0.4413931681178835), 202: (201.970643, 1.0), 204: (203.9734939, 0.23007367716008037)},
        'Ho': {165: (164.9303221, 1.0)},
        'Hs': {270: (270, 1.0)},
        'I': {127: (126.904473, 1.0)},
        'In': {113: (112.904058, 0.0448229025180232), 115: (114.903878, 1.0)},
        'Ir': {191: (190.960594, 0.594896331738437), 193: (192.9629264, 1.0)},
        'K': {39: (38.96370668, 1.0), 40: (39.96399848, 0.00012545827118502306), 41: (40.96182576, 0.07216745784012327)},
        'Kr': {78: (77.9203648, 0.006229490936529384), 80: (79.916379, 0.040114412058890624), 82: (81.9134836, 0.20343236176671875), 83: (82.914136, 0.20180041061996598), 84: (83.911507, 1.0), 86: (85.91061073, 0.3032095039219471)},
        'La': {138: (137.907112, 0.000900810729656691), 139: (138.9063533, 1.0)},
        'Li': {6: (6.015122795, 0.08213396818526132), 7: (7.01600455, 1.0)},
        'Lr': {262: (262, 1.0)},
        'Lu': {175: (174.9407718, 1.0), 176: (175.9426863, 0.026588645929576018)},
        'Md': {258: (258, 1.0)},
        'Mg': {24: (23.9850417, 1.0), 25: (24.98583692, 0.126598303582732), 26: (25.982592929, 0.13938473224458792)},
        'Mn': {55: (54.9380451, 1.0)},
        'Mo': {92: (91.906811, 0.6105828854898718), 94: (93.9050883, 0.38156262918561384), 95: (94.9058421, 0.6572964034725093), 96: (95.9046795, 0.6895411326994626), 97: (96.9060215, 0.39520463001240186), 98: (97.9054082, 1.0), 100: (99.907477, 0.3997519636213311)},
        'Mt': {276: (276, 1.0)},
        'N': {14: (14.0030740048, 1.0), 15: (15.0001088982, 0.0036532980047372433)},
        'Na': {23: (22.9897692809, 1.0)},
        'Nb': {93: (92.9063781, 1.0)},
        'Nd': {142: (141.9077233, 1.0), 143: (142.9098143, 0.44852941176470584), 144: (143.9100873, 0.8749999999999999), 145: (144.9125736, 0.3051470588235294), 146: (145.9131169, 0.6323529411764705), 148: (147.916893, 0.20955882352941177), 150: (149.920891, 0.20588235294117646)},
        'Ne': {20: (19.9924401754, 1.0), 21: (20.99384668, 0.0029840848806366046), 22: (21.991385114, 0.10223253757736515)},
        'Ni': {58: (57.9353429, 1.0), 60: (59.9307864, 0.3851982096717095), 61: (60.931056, 0.01674429946134445), 62: (61.9283451, 0.05338815369089956), 64: (63.927966, 0.013596388789736314)},
        'No': {259: (259, 1.0)},
        'Np': {237: (237, 1.0)},
        'O': {16: (15.99491461956, 1.0), 17: (16.9991317, 0.00038092564932786676), 18: (17.999161, 0.002054993634531913)},
        'Os': {184: (183.9524891, 0.0004904364884747426), 186: (185.9538382, 0.038989700833742036), 187: (186.9557505, 0.048062775870524765), 188: (187.9558382, 0.3246689553702795), 189: (188.9581475, 0.3960274644433546), 190: (189.958447, 0.6439431093673369), 192: (191.9614807, 1.0)},
        'P': {31: (30.97376163, 1.0)},
        'Pa': {231: (231.035884, 1.0)},
        'Pb': {204: (203.9730436, 0.026717557251908396), 206: (205.9744653, 0.45992366412213737), 207: (206.9758969, 0.4217557251908397), 208: (207.9766521, 1.0)},
        'Pd': {102: (101.905609, 0.03732162458836444), 104: (103.904036, 0.40761068422978414), 105: (104.905085, 0.8170508598609587), 106: (105.903486, 1.0), 108: (107.903892, 0.9681668496158069), 110: (109.905153, 0.42883278448591294)},
        'Pm': {145: (145, 1.0)},
        'Po': {209: (209, 1.0)},
        'Pr': {141: (140.9076528, 1.0)},
        'Pt': {190: (189.959932, 0.0004138094112083234), 192: (191.961038, 0.0231142113974935), 194: (193.9626803, 0.9744324899503429), 195: (194.9647911, 1.0), 196: (195.9649515, 0.7460983684086071), 198: (197.967893, 0.21172262946323006)},
        'Pu': {244: (244, 1.0)},
        'Ra': {226: (226, 1.0)},
        'Rb': {85: (84.911789738, 1.0), 87: (86.909180527, 0.38561729250381044)},
        'Re': {185: (184.952955, 0.597444089456869), 187: (186.9557531, 1.0)},
        'Rf': {265: (265, 1.0)},
        'Rg': {280: (280, 1.0)},
        'Rh': {103: (102.905504, 1.0)},
        'Rn': {222: (222, 1.0)},
        'Ru': {96: (95.907598, 0.17559429477020602), 98: (97.905287, 0.05927099841521395), 99: (98.9059393, 0.4044374009508716), 100: (99.9042195, 0.3993660855784469), 101: (100.9055821, 0.5407290015847861), 102: (101.9043493, 1.0), 104: (103.905433, 0.5901743264659272)},
        'S': {32: (31.972071, 1.0), 33: (32.97145876, 0.007895567954521529), 34: (33.9678669, 0.044741551742288665), 36: (35.96708076, 0.00010527423939362039)},
        'Sb': {121: (120.9038157, 1.0), 123: (122.904214, 0.7479461632581715)},
        'Sc': {45: (44.9559119, 1.0)},
        'Se': {74: (73.9224764, 0.017939931465430357), 76: (75.9192136, 0.18887321104616006), 77: (76.919914, 0.15379963716992542), 78: (77.9173091, 0.4791372707115501), 80: (79.9165213, 1.0), 82: (81.9166994, 0.17597258617214273)},
        'Sg': {271: (271, 1.0)},
        'Si': {28: (27.9769265325, 1.0), 29: (28.9764947, 0.0508007763789944), 30: (29.97377017, 0.03352742808193184)},
        'Sm': {144: (143.911999, 0.11476635514018692), 147: (146.9148979, 0.5603738317757009), 148: (147.9148227, 0.42018691588785045), 149: (148.9171847, 0.5166355140186916), 150: (149.9172755, 0.2758878504672897), 152: (151.9197324, 1.0), 154: (153.9222093, 0.8504672897196262)},
        'Sn': {112: (111.904818, 0.02977286678944138), 114: (113.902779, 0.020257826887661142), 115: (114.903342, 0.010435850214855739), 116: (115.901741, 0.4462860650705955), 117: (116.902952, 0.23572744014732966), 118: (117.901603, 0.743400859422959), 119: (118.903308, 0.26365868631062006), 120: (119.9021947, 1.0), 122: (121.903439, 0.1421117249846532), 124: (123.9052739, 0.17771639042357276)},
        'Sr': {84: (83.913425, 0.006781302978929523), 86: (85.9092602, 0.11939937030758052), 87: (86.9088771, 0.08476628723661904), 88: (87.9056121, 1.0)},
        'Ta': {180: (179.9474648, 0.00012001440172820739), 181: (180.9479958, 1.0)},
        'Tb': {159: (158.9253468, 1.0)},
        'Tc': {98: (98, 1.0)},
        'Te': {120: (119.90402, 0.002640845070422535), 122: (121.9030439, 0.07482394366197183), 123: (122.90427, 0.026115023474178403), 124: (123.9028179, 0.13908450704225353), 125: (124.9044307, 0.2074530516431925), 126: (125.9033117, 0.5528169014084507), 128: (127.9044631, 0.9313380281690141), 130: (129.9062244, 1.0)},
        'Th': {232: (232.0380553, 1.0)},
        'Ti': {46: (45.9526316, 0.11190992946283235), 47: (46.9517631, 0.10092240911557243), 48: (47.9479463, 1.0), 49: (48.94787, 0.07338578404774825), 50: (49.9447912, 0.07026587086272382)},
        'Tl': {203: (202.9723442, 0.4188422247446084), 205: (204.9744275, 1.0)},
        'Tm': {169: (168.9342133, 1.0)},
        'U': {234: (234.0409521, 5.43947974398182e-05), 235: (235.0439299, 0.007256668902897228), 238: (238.0507882, 1.0)},
        'Uuh': {293: (293, 1.0)},
        'Uuo': {294: (294, 1.0)},
        'Uup': {288: (288, 1.0)},
        'Uuq': {289: (289, 1.0)},
        'Uus': {292: (292, 1.0)},
        'Uut': {284: (284, 1.0)},
        'V': {50: (49.9471585, 0.002506265664160401), 51: (50.9439595, 1.0)},
        'W': {180: (179.946704, 0.0039164490861618795), 182: (181.9482042, 0.8648825065274152), 183: (182.950223, 0.46703655352480417), 184: (183.9509312, 1.0), 186: (185.9543641, 0.9278720626631853)},
        'Xe': {124: (123.905893, 0.0035379023806515393), 126: (125.904274, 0.003307492771827594), 128: (127.9035313, 0.07098845722185473), 129: (128.9047794, 0.9811212772124898), 130: (129.903508, 0.15128992218101278), 131: (130.9050824, 0.7890562868376654), 132: (131.9041535, 1.0), 134: (133.9053945, 0.3878202507748452), 136: (135.907219, 0.32916242390908484)},
        'Y': {89: (88.9058483, 1.0)},
        'Yb': {168: (167.933897, 0.004084197298146402), 170: (169.9347618, 0.09550738297203895), 171: (170.9363258, 0.44863336475023563), 172: (171.9363815, 0.6858309770656612), 173: (172.9382108, 0.5067546339930883), 174: (173.9388621, 1.0), 176: (175.9425717, 0.4008796732642161)},
        'Zn': {64: (63.9291422, 1.0), 66: (65.9260334, 0.5795765310350542), 67: (66.9271273, 0.08498384022540814), 68: (67.9248442, 0.39413275876357007), 70: (69.9253193, 0.013072843291621778)},
        'Zr': {90: (89.9047044, 1.0), 91: (90.9056458, 0.21807580174927113), 92: (91.9050408, 0.33333333333333337), 94: (93.9063152, 0.33780369290573375), 96: (95.9082734, 0.054421768707483)}
        }

def get_isotopes(composition, threshold=0.001, tolerance=0.01, mass_dict=_rel_mass):
    """Calculate isotope distribution from composition object using
    the binomial method.

    input
    -----
    composition: Composition object
    threshold: minimum abundance for isotope
    tolerance: tolerance for combining adjacent peaks in daltons
    mass_dict: dictionary containing relative masses

    returns
    -------
    masses: list of isotope masses
    abun: list of isotope relative abundances
    """

    masses = []
    abun = []

    for element in composition:
        ele_name, ele_num = re.match(r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$', element).groups()
        count = composition[element]

        if ele_num:
            ele_num = int(ele_num)
            if not ele_num in mass_dict[ele_name]:
                raise AssertionError('Can not produce isotope pattern, specified isotope {} not abundant'.format(element))
                continue
            isotopes = {0: mass_dict[ele_name][ele_num]}
        else:
            isotopes = mass_dict[element]

        for i in range(count):

            if len(masses) == 0:
                masses = [isotopes[iso][0] for iso in isotopes]
                abun = [isotopes[iso][1] for iso in isotopes]
                continue

            _masses = []
            _abun = []

            for peak in range(len(masses)):

                if abun[peak] < threshold:
                    continue

                for iso in isotopes:
                    _masses.append(masses[peak] + isotopes[iso][0])
                    _abun.append(abun[peak] * isotopes[iso][1])

            _masses, _abun = (list(buf) for buf in zip(*sorted(zip(_masses, _abun))))
            masses = [_masses[0]]
            abun = [_abun[0]]

            for peak in range(1, len(_masses)):
                if _masses[peak] <= (masses[-1] + tolerance):
                    _abundance = abun[-1] + _abun[peak]
                    masses[-1] = (masses[-1]*abun[-1] + _masses[peak]*_abun[peak]) / _abundance
                    abun[-1] = _abundance
                else:
                    masses.append(_masses[peak])
                    abun.append(_abun[peak])
            _max = max(abun)
            abun = [i/_max for i in abun]

    masses, abun = (list(buf) for buf in zip(*[(masses[i], ab) for i, ab in enumerate(abun) if ab >= threshold]))

    return np.array(masses), np.array(abun)

def get_colours(colour, arr, threshold):
    """Convert array of data to colours ranging from colour to white.
    Colour must be RGB 0-1 tuple. Alt_colour is supplied for data that
    may stray into negative numbers"""

    multiplier = 0.25+(0.75/(1-threshold)*(arr - threshold))
    colours = np.zeros((len(arr), 3))
    colours[arr > threshold] = np.array(colour)*multiplier[arr > threshold, np.newaxis] + 1 - multiplier[arr > threshold, np.newaxis]
    colours[arr < threshold] = 0.9

    return colours

def get_figure_rows(peptide_ids):
    """Fix overlaps in peptides and
    return numpy array of rows for plotting
    """

    plotted = []
    rows_ = np.zeros(len(peptide_ids))
    row = 0
    last = 0
    while len(plotted) < len(peptide_ids):
        for i, pep in enumerate(peptide_ids):
            if plotted:
                if i not in plotted and pep[0] > last:
                    plotted.append(i)
                    rows_[i] = row
                    last = peptide_ids[i][1]
            else:
                plotted.append(i)
        row += 1
        last = 0

    return rows_

#


class PDB():
    """
    Class for parsing PDB files that follow the accepted format
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    Extracts ATOM and HETATM records.
    HETATM with HOH resname are removed.

    Dictionary format:

    Structure (list)>
        Model (odict)>
            Chain (odict)>
                Residue (list)>
                    Atom (dict)
    """

    def __init__(self, filename=None):
        if not os.path.splitext(filename)[1] == '.pdb':
            raise ValueError('Incorrect file extension, must be .pdb')
        else:
            self.filename = filename

        self.AA_3to1    = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E',
                           'GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                           'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
                           'TYR':'Y','VAL':'V'}

        self.ATOM_STRING = "{}{:5d} {:4}{:.1}{:.3} {:.1}{:>4d}{:.1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:.2}  \n"

        self.structure = [None]
        self.num_models = 0
        self.num_chains = 0
        self.chains = {}
        self._parse()

# FIXME: Add ability to handle PDBs with gaps
    def _parse(self):
        """Parse the PDB file as a series of record entries into a dictionary"""

        with open(self.filename, 'r') as entries:
            model = 0
            open_model = False
            chain = None
            resnum = None
            for entry in entries:
                record_type = entry[0:6]

                if record_type == 'ATOM  ' or record_type == 'HETATM' and not entry[17:20] == 'HOH':
                    if not open_model:
                        model += 1
                        open_model = True
                        self.structure.append(OrderedDict())
                    if not entry[21] == chain:
                        chain = entry[21]
                        if chain == ' ':
                            chain = 'A'
                        if not chain in self.structure[model]:
                            self.structure[model][chain] = OrderedDict()
                    if not int(entry[22:26]) == resnum:
                        resnum = int(entry[22:26])
                        self.structure[model][chain][resnum] = []
                    self.structure[model][chain][resnum].append({'type': record_type,
                                                                 'serial': int(entry[6:11]),
                                                                 'name': entry[12:16],
                                                                 'altLoc': entry[16],
                                                                 'resname': entry[17:20],
                                                                 'icode': entry[26],
                                                                 'x': float(entry[30:38]),
                                                                 'y': float(entry[38:46]),
                                                                 'z': float(entry[46:54]),
                                                                 'occupancy': float(entry[54:60]),
                                                                 'bfactor': -1.00,
                                                                 'element': entry[76:78]
                                                                 })
                elif record_type == 'MODEL':
                    open_model = True
                    model = int(entry[10:14])
                elif record_type == 'ENDMDL':
                    open_model = False
                    chain = None
                    resnum = None

            self.num_models = model
            self.chains = {chain:self._get_sequence(self.structure[1][chain]) for chain in self.structure[1]}
            self.num_chains = len(self.chains)

    def _get_sequence(self, chain):
        """Parse single character amino acid code from chain"""
        _sequence = ''

        for residue in chain:
            resname = chain[residue][0]['resname']
            if chain[residue][0]['type'] == 'ATOM  ':
                _sequence += self.AA_3to1[resname]

        return _sequence


# FIXME: Add ability to check for alignment with sequence and colour when gaps present
    def bfactor_by_residue(self, sequence, bfactors):
        """
        Modify b-factor of atoms in self.structure by residue.
        Takes input sequence to compare against self.chains and a
        list of corresponding b-factors
        """
        matches = []
        for chain in self.chains:
            if self.chains[chain] == sequence:
                matches.append(chain)

        if matches == []:
            raise ValueError('No chains match input sequence, please check that sequence matches PDB.')
            return

        for match in matches:
            for model in self.structure[1:]:
                for residue in model[match].items():
                    if residue[1][0]['type'] == 'ATOM  ': #Only count ATOM entries
                        bfactor = bfactors[residue[0]-1] #Account for one-indexing
                        for atom in residue[1]:
                            atom['bfactor'] = bfactor

    def _to_string(self):
        """Convert self.structure into string with PDB formatting"""

        _PDB_string = ''
        # Add REMARK for bookkeeping purposes
        _PDB_string += 'REMARK     Output from pepFoot.\n'
        _PDB_string += 'REMARK     {}\n'.format(os.path.split(self.filename)[1])
        _PDB_string += 'REMARK     {}\n'.format(datetime.today().strftime('%d %b %Y  %I:%M%p'))
        for i, model in enumerate(self.structure[1:]):
            _PDB_string += 'MODEL     {}\n'.format(i+1)
            for chain in model.items():
                for residue in chain[1].items():
                    for atom in residue[1]:
                        _line = self.ATOM_STRING.format(atom['type'], atom['serial'], atom['name'],
                                                    atom['altLoc'], atom['resname'],chain[0],
                                                    residue[0], atom['icode'], atom['x'],
                                                    atom['y'], atom['z'], atom['occupancy'],
                                                    atom['bfactor'], atom['element'])
                        _PDB_string += _line

            _PDB_string += 'ENDMDL\n'
        _PDB_string += 'END\n'

        return _PDB_string

    def write(self, filename):
        """Write self.structure to file"""

        with open(filename, 'w') as out_file:
            lines = self._to_string()
            out_file.write(lines)
