#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED_au565 {

static constexpr std::array<sunindextype, 779> dwdx_colptrs_SPARCED_au565_ = {
    0, 1, 5, 10, 18, 22, 26, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 81, 83, 85, 88, 90, 91, 92, 93, 104, 115, 117, 129, 133, 136, 141, 144, 148, 152, 161, 166, 169, 172, 175, 178, 181, 185, 189, 198, 201, 208, 210, 215, 217, 220, 223, 227, 230, 238, 241, 245, 248, 254, 256, 258, 264, 265, 266, 272, 274, 277, 280, 282, 284, 286, 289, 292, 294, 296, 300, 303, 307, 309, 311, 314, 317, 321, 323, 326, 328, 331, 335, 338, 340, 343, 344, 346, 349, 352, 359, 361, 364, 367, 369, 373, 378, 380, 384, 386, 390, 392, 393, 395, 398, 400, 402, 405, 407, 409, 412, 415, 417, 420, 422, 424, 428, 431, 434, 436, 438, 439, 443, 445, 447, 449, 451, 453, 457, 460, 462, 465, 467, 470, 472, 475, 484, 499, 503, 505, 509, 512, 515, 526, 528, 538, 540, 549, 560, 562, 566, 570, 572, 576, 578, 580, 583, 586, 590, 594, 596, 599, 602, 605, 608, 611, 614, 617, 620, 630, 640, 650, 654, 657, 661, 666, 670, 676, 682, 688, 693, 699, 702, 706, 711, 716, 720, 726, 732, 738, 744, 748, 754, 759, 764, 769, 775, 781, 784, 790, 795, 800, 805, 810, 814, 818, 821, 824, 827, 834, 841, 848, 855, 862, 869, 876, 883, 890, 897, 904, 911, 918, 925, 932, 939, 946, 953, 960, 967, 974, 981, 988, 995, 1002, 1006, 1010, 1014, 1018, 1024, 1030, 1036, 1042, 1044, 1046, 1048, 1050, 1052, 1054, 1056, 1058, 1060, 1062, 1064, 1066, 1068, 1070, 1072, 1074, 1076, 1078, 1080, 1082, 1084, 1086, 1088, 1090, 1092, 1094, 1096, 1098, 1100, 1103, 1106, 1109, 1112, 1115, 1118, 1121, 1124, 1127, 1130, 1133, 1136, 1139, 1142, 1145, 1148, 1151, 1154, 1157, 1160, 1163, 1166, 1169, 1172, 1175, 1178, 1181, 1184, 1187, 1190, 1193, 1196, 1199, 1202, 1205, 1208, 1211, 1214, 1217, 1220, 1223, 1226, 1229, 1232, 1235, 1238, 1241, 1244, 1247, 1250, 1253, 1256, 1259, 1262, 1265, 1268, 1271, 1274, 1277, 1280, 1283, 1286, 1289, 1292, 1295, 1298, 1301, 1304, 1307, 1310, 1313, 1316, 1319, 1322, 1325, 1328, 1331, 1334, 1337, 1340, 1343, 1346, 1349, 1352, 1355, 1358, 1361, 1364, 1367, 1370, 1373, 1376, 1379, 1382, 1385, 1388, 1391, 1394, 1397, 1400, 1403, 1406, 1409, 1412, 1415, 1418, 1421, 1424, 1427, 1430, 1433, 1436, 1439, 1442, 1445, 1448, 1451, 1454, 1457, 1460, 1463, 1466, 1469, 1472, 1475, 1478, 1481, 1484, 1487, 1490, 1493, 1496, 1499, 1502, 1505, 1508, 1511, 1514, 1517, 1520, 1523, 1526, 1529, 1532, 1535, 1538, 1541, 1544, 1547, 1550, 1553, 1556, 1559, 1562, 1565, 1568, 1571, 1574, 1577, 1580, 1583, 1586, 1589, 1592, 1595, 1598, 1601, 1604, 1607, 1610, 1613, 1616, 1619, 1622, 1625, 1628, 1631, 1634, 1637, 1640, 1643, 1646, 1649, 1652, 1655, 1658, 1661, 1664, 1667, 1670, 1673, 1676, 1679, 1682, 1685, 1688, 1691, 1694, 1697, 1700, 1703, 1706, 1709, 1712, 1715, 1718, 1721, 1724, 1727, 1730, 1733, 1736, 1739, 1742, 1745, 1748, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1778, 1781, 1784, 1787, 1790, 1793, 1796, 1799, 1802, 1805, 1808, 1811, 1814, 1817, 1820, 1823, 1826, 1829, 1832, 1835, 1838, 1841, 1844, 1847, 1850, 1853, 1856, 1859, 1862, 1865, 1868, 1871, 1874, 1877, 1880, 1883, 1886, 1889, 1892, 1895, 1898, 1901, 1904, 1907, 1910, 1913, 1916, 1919, 1922, 1925, 1928, 1931, 1934, 1937, 1940, 1943, 1946, 1949, 1952, 1955, 1958, 1961, 1964, 1967, 1970, 1973, 1976, 1979, 1982, 1985, 1988, 1991, 1994, 1997, 2000, 2003, 2006, 2009, 2012, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042, 2045, 2048, 2051, 2054, 2057, 2060, 2063, 2066, 2069, 2072, 2075, 2078, 2081, 2084, 2087, 2090, 2093, 2096, 2099, 2102, 2105, 2108, 2111, 2114, 2117, 2120, 2123, 2126, 2129, 2132, 2135, 2138, 2141, 2144, 2147, 2150, 2153, 2156, 2165, 2195, 2225, 2227, 2288, 2290, 2320, 2322, 2324, 2356, 2358, 2388, 2391, 2396, 2426, 2428, 2431, 2433, 2436, 2441, 2502, 2505, 2507, 2509, 2513, 2517, 2519, 2523, 2527, 2532, 2537, 2539, 2542, 2547, 2549, 2552, 2556, 2559, 2562, 2565, 2568, 2570, 2572, 2573, 2575, 2578, 2581, 2582, 2642, 2648, 2650, 2653, 2656, 2658, 2666, 2668, 2671, 2674, 2676, 2679, 2682, 2685, 2688, 2690, 2692, 2696, 2698, 2702, 2704, 2707, 2709, 2712, 2714, 2716, 2720, 2734, 2737, 2739, 2740, 2742, 2744, 2746, 2749, 2752, 2754, 2756, 2759, 2762, 2765, 2768, 2771, 2774, 2777, 2780, 2783, 2786, 2789, 2793, 2796, 2799, 2802, 2805, 2808, 2811, 2814, 2817, 2820, 2821, 2824, 2827, 2830, 2833, 2836, 2839, 2842, 2846, 2850, 2853, 2856, 2978, 2980, 2984, 2987, 2990, 2992, 2993, 2995, 2996, 2998, 3000, 3002, 3006, 3007, 3008, 3009, 3010
};

void dwdx_colptrs_SPARCED_au565(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_SPARCED_au565_));
}
} // namespace model_SPARCED_au565
} // namespace amici
