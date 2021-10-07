#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCED_au565_k.h"

namespace amici {
namespace model_SPARCED_au565 {

void x0_SPARCED_au565(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1897.3920000000001;
    x0[1] = 250.45461826117199;
    x0[2] = 7.9004432331299599;
    x0[3] = 162.51811464184499;
    x0[4] = 2.2139389155153801;
    x0[7] = 7.9004432331288204;
    x0[8] = 7.9004432331268504;
    x0[9] = 7.9004432331267296;
    x0[10] = 7.9004432331153698;
    x0[11] = 7.9004432331149204;
    x0[12] = 7.9004432331138696;
    x0[13] = 7.9004432331137302;
    x0[14] = 7.9004432331109999;
    x0[15] = 7.9004432331076;
    x0[16] = 7.9004432331031103;
    x0[17] = 7.9004432331272598;
    x0[18] = 7.9004432331231396;
    x0[19] = 7.9004432331175698;
    x0[20] = 7.9004432331142302;
    x0[21] = 7.9004432331096304;
    x0[22] = 7.9004432331000203;
    x0[23] = 7.9004432330964498;
    x0[24] = 7.9004432330942098;
    x0[25] = 7.9004432330904599;
    x0[26] = 7.90044323308903;
    x0[27] = 1.52773183218054;
    x0[28] = 47.949320187152701;
    x0[29] = 5.3101898808394203;
    x0[30] = 0.00012923990043548601;
    x0[31] = 1.85360238606969e-5;
    x0[34] = 4.5680949831399902;
    x0[35] = 238.81109984977701;
    x0[37] = 16.6697717051225;
    x0[38] = 1.9734535624212499;
    x0[39] = 7968.3249990129898;
    x0[40] = 2.9880985707322401;
    x0[41] = 0.031654087573086402;
    x0[42] = 1128.97873391962;
    x0[43] = 0.66341584324008296;
    x0[44] = 1.4569068029852701;
    x0[45] = 17.458972112442702;
    x0[46] = 3.5643860819979198;
    x0[47] = 0.25539278480611499;
    x0[48] = 0.89006573310397696;
    x0[49] = 2.74857991478707;
    x0[50] = 11.408668449590101;
    x0[51] = 1.1222882917045101;
    x0[52] = 141.29335044365601;
    x0[53] = 49.248916919292903;
    x0[54] = 16.019953924816502;
    x0[55] = 55.298216779885003;
    x0[56] = 4.0298702578317496;
    x0[57] = 14.044443064724099;
    x0[58] = 5.2088806264855503;
    x0[59] = 1.2761204897180001;
    x0[60] = 20.415784268573098;
    x0[61] = 0.594325853633206;
    x0[62] = 108.811400906412;
    x0[63] = 11.5836766345124;
    x0[64] = 4.2047887824151902;
    x0[65] = 41.501882267316397;
    x0[66] = 3.44351291813045;
    x0[67] = 3.90106612835176;
    x0[68] = 0.68426294396204101;
    x0[69] = 71.167470935605607;
    x0[70] = 0.73013601962133301;
    x0[71] = 9.4030028408023494;
    x0[72] = 26.518895499505;
    x0[73] = 16.740411350390598;
    x0[74] = 0.083818596756356001;
    x0[75] = 0.14977203180328799;
    x0[76] = 8996.0694688458807;
    x0[77] = 1.6867505849964799;
    x0[78] = 19.3797114902484;
    x0[79] = 2.0300996279564099;
    x0[80] = 0.21306676558587001;
    x0[81] = 0.24227189683255301;
    x0[82] = 0.015912978193499999;
    x0[84] = 5.0356028690955696;
    x0[87] = 0.52560696236975002;
    x0[89] = 15.7579351303589;
    x0[91] = 0.00060695062990999697;
    x0[92] = 0.000172796864583068;
    x0[94] = 86.430072592415399;
    x0[95] = 0.0082842417969008;
    x0[96] = 9.6101627306187403e-5;
    x0[97] = 4.8438637065825203;
    x0[99] = 2.7634792452749302e-5;
    x0[101] = 13.7454716515402;
    x0[102] = 8.2689118186216802e-5;
    x0[103] = 321.99326141996102;
    x0[105] = 0.38698462900981301;
    x0[106] = 63.633675602568601;
    x0[107] = 0.000110739584481628;
    x0[108] = 2.8922966219346299e-6;
    x0[109] = 36.7964766359013;
    x0[110] = 0.000326273677951481;
    x0[111] = 117.460830285151;
    x0[113] = 0.00013667903230321101;
    x0[114] = 0.0018475001448977301;
    x0[115] = 5.2549780483625002;
    x0[116] = 0.0297637987959283;
    x0[117] = 0.000107924210869074;
    x0[118] = 0.0017386923978764;
    x0[120] = 5.9332370014618998e-6;
    x0[121] = 2259.9696867745201;
    x0[122] = 2.6292297456222299e-6;
    x0[123] = 0.00026292282423237198;
    x0[124] = 17388.8342854031;
    x0[125] = 0.026285015317822999;
    x0[126] = 0.26278073910577199;
    x0[127] = 1827.63461886191;
    x0[128] = 0.0027626581037228401;
    x0[129] = 0.018944780566453601;
    x0[130] = 0.018392221681234999;
    x0[131] = 0.082254837081730398;
    x0[132] = 2.3894681208266302e-6;
    x0[133] = 0.39427098167140501;
    x0[137] = 0.00071868708817458599;
    x0[139] = 0.21258726019113799;
    x0[140] = 6.0891180052038596;
    x0[162] = 53.668380059968896;
    x0[163] = 26.6738175027238;
    x0[170] = 14.1373876930511;
    x0[171] = 0.00042236335782909099;
    x0[172] = 2.0347283726983099;
    x0[173] = 0.029297108863339999;
    x0[174] = 0.029786378667274599;
    x0[175] = 2.8802513859916798;
    x0[184] = 0.199864422314662;
    x0[185] = 0.0041400924463773402;
    x0[186] = 1.31105706142958;
    x0[187] = 0.30873598713625799;
    x0[225] = 2.6677858258338998;
    x0[642] = 1.59625294392691;
    x0[643] = 6.9843494952573896e-5;
    x0[644] = 5.29474048049416;
    x0[645] = 80.400155089168095;
    x0[646] = 0.093319209221500596;
    x0[647] = 0.46336130863057301;
    x0[648] = 21.322577404081098;
    x0[649] = 0.032436170047875197;
    x0[650] = 13.4332909399619;
    x0[651] = 4.1572832993659601;
    x0[653] = 3.53547422212801;
    x0[659] = 0.0015031101410492599;
    x0[661] = 2.4750061956951401;
    x0[662] = 262.42133273477202;
    x0[663] = 1.9430129974389201;
    x0[665] = 35.165139279178803;
    x0[666] = 7.0810126973175702;
    x0[667] = 1.74781847969143;
    x0[668] = 0.40092822451774901;
    x0[669] = 0.172591713640279;
    x0[670] = 55.140641814642002;
    x0[671] = 16.3998242257436;
    x0[672] = 48.767985332669099;
    x0[674] = 2.0048002470142001;
    x0[675] = 36.086542048847001;
    x0[676] = 158.33623444139999;
    x0[677] = 7.8707265257343799;
    x0[678] = 19.974573209629501;
    x0[683] = 19.875662682114601;
    x0[685] = 0.00283840258368856;
    x0[686] = 33.8630460910476;
    x0[687] = 16.944082802292598;
    x0[688] = 181.679600977175;
    x0[690] = 9998.31573086957;
    x0[691] = 1.6081692611652001;
    x0[692] = 2.3138978729505899;
    x0[693] = 0.39790044354362902;
    x0[694] = 26.050572781416601;
    x0[695] = 0.020998664423034798;
    x0[696] = 0.20999872727133401;
    x0[697] = 10.9314973201974;
    x0[698] = 0.175784646857296;
    x0[699] = 0.023009707790052202;
    x0[700] = 0.079074654102476993;
    x0[701] = 3.7192690387263001;
    x0[702] = 2.7424036765999902;
    x0[703] = 0.0021350199371874602;
    x0[704] = 41.385459762536499;
    x0[705] = 0.78961480959807095;
    x0[706] = 0.025203004338202498;
    x0[707] = 0.90340290212743202;
    x0[708] = 0.34165881442171198;
    x0[709] = 0.22756286233334;
    x0[710] = 21.0748745172165;
    x0[713] = 342.36608861183799;
    x0[715] = 227.05161152573601;
    x0[716] = 49.344995597795801;
    x0[717] = 54.680331674173203;
    x0[718] = 0.22412483827677701;
    x0[719] = 0.205420021476072;
    x0[720] = 66.117786658081002;
    x0[723] = 26.057403172015398;
    x0[725] = 1.5150247255992499;
    x0[727] = 0.00019413105418770801;
    x0[728] = 0.0046387940403034098;
    x0[729] = 0.35198975995533999;
    x0[731] = 0.0080149440940076797;
    x0[733] = 0.15861337666442099;
    x0[734] = 0.047174487121356502;
    x0[735] = 5.5364153303104198;
    x0[736] = 1.2032259461447199;
    x0[737] = 7.8707734583224704;
    x0[743] = 0.0033828535843692499;
    x0[744] = 0.00023314918760366699;
    x0[745] = 0.00021033958926375401;
    x0[746] = 0.0079007705689829906;
    x0[747] = 0.000948564569187201;
    x0[748] = 0.024698524128223698;
    x0[752] = 0.00021393492992434501;
    x0[757] = 0.0987677634781126;
    x0[758] = 0.00097549161922655598;
    x0[759] = 0.00048899260833589701;
    x0[760] = 0.00014543524631712099;
    x0[761] = 102.908410444059;
    x0[762] = 15.543821698904599;
    x0[763] = 0.30547349043127098;
    x0[764] = 0.84219616408309805;
    x0[765] = 0.25048437232536402;
    x0[766] = 6.7209592293805596;
    x0[771] = 62.343545279350103;
    x0[772] = 64.155683083277296;
    x0[773] = 0.0022137362333277998;
    x0[774] = 0.0161286496999597;
    x0[775] = 0.00063249606666508596;
    x0[776] = 0.00347872836665797;
    x0[777] = 0.0031624803333254298;
    x0[778] = 0.00537621656665323;
    x0[779] = 0.00537621656665323;
    x0[780] = 0.00537621656665323;
    x0[781] = 0.00537621656665323;
    x0[782] = 0.020556122166615301;
    x0[783] = 0.00094874409999762797;
    x0[784] = 0.00126499213333017;
    x0[785] = 0.00537621656665323;
    x0[786] = 0.00537621656665323;
    x0[787] = 0.00537621656665323;
    x0[788] = 0.00537621656665323;
    x0[789] = 0.00537621656665323;
    x0[790] = 0.00537621656665323;
    x0[791] = 0.00537621656665323;
    x0[792] = 0.00537621656665323;
    x0[793] = 0.00537621656665323;
    x0[794] = 0.00537621656665323;
    x0[795] = 0.00537621656665323;
    x0[796] = 0.00537621656665323;
    x0[797] = 0.00537621656665323;
    x0[798] = 0.0047437204999881402;
    x0[799] = 0.011068681166639;
    x0[800] = 0.00126499213333017;
    x0[801] = 0.0113849291999715;
    x0[802] = 0.00094874409999762797;
    x0[804] = 0.00031624803333254298;
    x0[805] = 0.0037949763999905101;
    x0[806] = 0.0056924645999857698;
    x0[807] = 0.0031624803333254298;
    x0[808] = 0.00094874409999762797;
    x0[809] = 0.0072737047666484797;
    x0[810] = 0.00537621656665323;
    x0[811] = 0.0028462322999928801;
    x0[812] = 0.0075899527999810298;
    x0[813] = 0.0091711929666437406;
    x0[814] = 0.0079062008333135696;
    x0[815] = 0.0044274724666555996;
    x0[816] = 0.00031624803333254298;
    x0[817] = 0.00094874409999762797;
    x0[818] = 0.0028462322999928801;
    x0[819] = 0.00126499213333017;
    x0[820] = 0.10341310689974099;
    x0[821] = 0.0031624803333254298;
    x0[822] = 0.00063249606666508596;
    x0[823] = 0.00063249606666508596;
    x0[824] = 0.0056924645999857698;
    x0[825] = 0.00031624803333254298;
    x0[826] = 0.00347872836665797;
    x0[827] = 0.00031624803333254298;
    x0[830] = 0.00126499213333017;
    x0[831] = 0.00031624803333254298;
    x0[832] = 0.00031624803333254298;
    x0[835] = 0.0044274724666555996;
    x0[837] = 0.00063249606666508596;
    x0[838] = 0.00094874409999762797;
    x0[840] = 0.0022137362333277998;
    x0[841] = 0.00094874409999762797;
    x0[842] = 0.0060087126333183104;
    x0[843] = 0.00094874409999762797;
    x0[845] = 0.00063249606666508596;
    x0[846] = 0.00126499213333017;
    x0[848] = 0.00063249606666508596;
    x0[849] = 0.00126499213333017;
    x0[850] = 0.0028462322999928801;
    x0[851] = 0.00252998426666034;
    x0[852] = 0.00094874409999762797;
    x0[853] = 0.0047437204999881402;
    x0[854] = 0.00031624803333254298;
    x0[855] = 0.00031624803333254298;
    x0[856] = 0.01865863396662;
    x0[857] = 0.0028462322999928801;
    x0[858] = 0.00063249606666508596;
    x0[859] = 0.0018974881999952601;
    x0[860] = 0.00252998426666034;
    x0[861] = 0.00031624803333254298;
    x0[862] = 0.0018974881999952601;
    x0[863] = 0.0015812401666627099;
    x0[864] = 0.00094874409999762797;
    x0[865] = 0.00031624803333254298;
    x0[866] = 0.00094874409999762797;
    x0[867] = 0.0022137362333277998;
    x0[868] = 0.00031624803333254298;
    x0[869] = 0.0015812401666627099;
    x0[870] = 0.00031624803333254298;
    x0[871] = 0.0018974881999952601;
    x0[872] = 0.00126499213333017;
    x0[873] = 0.0192911300332851;
    x0[874] = 0.00252998426666034;
    x0[875] = 0.00063249606666508596;
    x0[876] = 0.00094874409999762797;
    x0[877] = 0.00126499213333017;
    x0[878] = 0.00126499213333017;
    x0[879] = 0.0022137362333277998;
    x0[880] = 0.0015812401666627099;
    x0[881] = 0.0031624803333254298;
    x0[882] = 0.00063249606666508596;
    x0[883] = 0.0047437204999881402;
    x0[884] = 0.00031624803333254298;
    x0[886] = 0.00063249606666508596;
    x0[887] = 0.028462322999928801;
    x0[888] = 0.0037949763999905101;
    x0[889] = 0.00094874409999762797;
    x0[890] = 0.00126499213333017;
    x0[891] = 0.0113849291999715;
    x0[892] = 0.00031624803333254298;
    x0[893] = 0.0015812401666627099;
    x0[894] = 0.00031624803333254298;
    x0[895] = 0.0094874409999762804;
    x0[896] = 0.0015812401666627099;
    x0[897] = 0.00031624803333254298;
    x0[898] = 0.0044274724666555996;
    x0[899] = 0.0031624803333254298;
    x0[900] = 0.00094874409999762797;
    x0[903] = 0.0246673465999383;
    x0[904] = 0.00126499213333017;
    x0[905] = 0.00126499213333017;
    x0[908] = 0.00094874409999762797;
    x0[909] = 0.00347872836665797;
    x0[910] = 0.0015812401666627099;
    x0[911] = 0.00031624803333254298;
    x0[912] = 0.0015812401666627099;
}

} // namespace model_SPARCED_au565
} // namespace amici