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
    x0[1] = 257.97745903797602;
    x0[2] = 7.5879057271332604;
    x0[3] = 169.172532940362;
    x0[4] = 2.2164381536211701;
    x0[7] = 7.5879057271449497;
    x0[8] = 7.5879057271651504;
    x0[9] = 7.5879057271663797;
    x0[10] = 7.58790572728262;
    x0[11] = 7.5879057272872803;
    x0[12] = 7.58790572729795;
    x0[13] = 7.5879057272994199;
    x0[14] = 7.5879057273274002;
    x0[15] = 7.5879057273622301;
    x0[16] = 7.5879057274081303;
    x0[17] = 7.5879057271609502;
    x0[18] = 7.58790572720316;
    x0[19] = 7.5879057272601296;
    x0[20] = 7.58790572729436;
    x0[21] = 7.5879057273413899;
    x0[22] = 7.5879057274398196;
    x0[23] = 7.5879057274763397;
    x0[24] = 7.5879057274992503;
    x0[25] = 7.5879057275375796;
    x0[26] = 7.5879057275522497;
    x0[27] = 0.49684511458747299;
    x0[28] = 66.574203604846304;
    x0[29] = 14.059279609982401;
    x0[30] = 0.00039740007246086799;
    x0[31] = 1.22434425848785e-5;
    x0[34] = 52.7206194223889;
    x0[35] = 151.82689432238499;
    x0[37] = 7.4874616258725997;
    x0[38] = 1.6950145270872601;
    x0[39] = 7970.2084661261197;
    x0[40] = 3.2976014985085902;
    x0[41] = 0.0341160128634291;
    x0[42] = 1128.72922195657;
    x0[43] = 0.68429384555805295;
    x0[44] = 1.5022539550981999;
    x0[45] = 18.820595214769099;
    x0[46] = 3.8424172322215502;
    x0[47] = 0.25487426563927901;
    x0[48] = 0.88981667418199095;
    x0[49] = 2.7476092434044799;
    x0[50] = 11.4087020666635;
    x0[51] = 1.1219051240578;
    x0[52] = 141.29061296089;
    x0[53] = 49.246160962195802;
    x0[54] = 16.0212696111391;
    x0[55] = 55.298294313577699;
    x0[56] = 4.0227435472208599;
    x0[57] = 14.0441965582449;
    x0[58] = 5.2098706285527498;
    x0[59] = 1.2763770430989001;
    x0[60] = 20.425626432638801;
    x0[61] = 0.59447346152757896;
    x0[62] = 108.811152278902;
    x0[63] = 11.5850959529162;
    x0[64] = 4.20471661099337;
    x0[65] = 41.507400735924499;
    x0[66] = 3.4436651934934202;
    x0[67] = 3.90060674216903;
    x0[68] = 0.68419238749168898;
    x0[69] = 71.167510625785894;
    x0[70] = 0.73005533784832;
    x0[71] = 9.4028985967315393;
    x0[72] = 26.518949797775299;
    x0[73] = 16.740396907084101;
    x0[74] = 0.083899711906506802;
    x0[75] = 0.14977203180328799;
    x0[76] = 8991.9741146814904;
    x0[77] = 1.86017388633794;
    x0[78] = 2.8981458594063998;
    x0[79] = 0.32803917955686801;
    x0[80] = 0.031870640811638297;
    x0[81] = 0.0363227670454703;
    x0[82] = 0.0023724698853255801;
    x0[84] = 5.0403742099450204;
    x0[87] = 0.77074498490576104;
    x0[89] = 15.7556232368406;
    x0[91] = 0.00060678965907372403;
    x0[92] = 0.0028629039676284799;
    x0[93] = 5.4802967992892598e-6;
    x0[94] = 86.427004863884903;
    x0[95] = 0.0082817638512529598;
    x0[96] = 9.5934475964226006e-5;
    x0[97] = 4.8445258513871403;
    x0[99] = 2.7572380420539399e-5;
    x0[101] = 13.765333502817899;
    x0[102] = 8.2664571289245694e-5;
    x0[103] = 322.321717901189;
    x0[105] = 0.38645263047267697;
    x0[106] = 63.699681615341099;
    x0[107] = 0.000110752327496312;
    x0[108] = 2.8889977779826299e-6;
    x0[109] = 36.846428184867897;
    x0[110] = 0.00032634395662796699;
    x0[111] = 117.514183158755;
    x0[113] = 0.000136480184794505;
    x0[114] = 0.0018446943254569301;
    x0[115] = 5.2621117443592702;
    x0[116] = 0.029758939513901201;
    x0[117] = 0.000107596641014086;
    x0[118] = 0.00173576827644512;
    x0[120] = 5.9052799303181403e-6;
    x0[121] = 2259.9698588209199;
    x0[122] = 2.6132935933620202e-6;
    x0[123] = 0.00026132920391362099;
    x0[124] = 17392.510413294502;
    x0[125] = 0.026131220498285701;
    x0[126] = 0.26124319681820501;
    x0[127] = 1827.7031282657199;
    x0[128] = 0.0027460160905985001;
    x0[129] = 0.018826032059744102;
    x0[130] = 0.018284607926628799;
    x0[131] = 0.082774625310500102;
    x0[132] = 2.39049850929638e-6;
    x0[133] = 0.39444097260002497;
    x0[137] = 0.00071370996112675802;
    x0[139] = 0.21142008400526499;
    x0[140] = 6.0873856291964801;
    x0[162] = 28.809707256709601;
    x0[163] = 53.107954175722597;
    x0[164] = 0.00022315588573685799;
    x0[165] = 0.00041143492587758398;
    x0[166] = 5.1248091512557902e-5;
    x0[170] = 14.1575100646776;
    x0[171] = 0.00021605423937469601;
    x0[172] = 2.0377035913480102;
    x0[173] = 0.029318647051415899;
    x0[174] = 0.0121525427824524;
    x0[175] = 0.82998665893918699;
    x0[176] = 6.4289605358221703e-6;
    x0[177] = 1.4764207718115201e-6;
    x0[184] = 0.200433779035015;
    x0[185] = 0.0041522087425835;
    x0[186] = 1.3129854514461901;
    x0[187] = 0.051390892433901297;
    x0[225] = 5.3115999378286203;
    x0[226] = 4.1143352691067902e-5;
    x0[642] = 1.5473681324796;
    x0[643] = 7.1832960015396496e-5;
    x0[644] = 6.1226991104817197;
    x0[645] = 78.972534882181805;
    x0[646] = 0.106881755854983;
    x0[647] = 1.968364450952;
    x0[648] = 21.354016131336898;
    x0[649] = 0.029372783431205401;
    x0[650] = 14.869659356257101;
    x0[651] = 4.16719095413114;
    x0[653] = 3.54006341390614;
    x0[659] = 0.0042004386567622096;
    x0[661] = 15.8546205217425;
    x0[662] = 1653.8852097858801;
    x0[663] = 1.90224763153794;
    x0[665] = 21.8626835240389;
    x0[666] = 1.1869498337937301;
    x0[667] = 1.85847306712872;
    x0[668] = 0.025450648021143801;
    x0[669] = 0.073380669410695404;
    x0[670] = 491.72819885208003;
    x0[671] = 115.023713085164;
    x0[672] = 269.01838710161002;
    x0[674] = 7.4357565373497501;
    x0[675] = 133.84412776517101;
    x0[676] = 118.54960199062;
    x0[677] = 21.8569042302098;
    x0[678] = 55.469128792993899;
    x0[683] = 19.902897358434299;
    x0[685] = 0.64672107843670501;
    x0[686] = 0.171526037878763;
    x0[687] = 0.085826864742217596;
    x0[688] = 210.135591414076;
    x0[690] = 9996.2884031562699;
    x0[691] = 3.5492737931415799;
    x0[692] = 0.44348169077844701;
    x0[693] = 0.78688001258987506;
    x0[694] = 24.7111591560178;
    x0[695] = 0.090069145019649396;
    x0[696] = 0.90074323995501604;
    x0[697] = 10.7423337835877;
    x0[698] = 0.38124804216574998;
    x0[699] = 0.098574070574033695;
    x0[700] = 0.079068744225070403;
    x0[701] = 3.7244860637583299;
    x0[702] = 2.7460267001808201;
    x0[703] = 0.0093332081997871495;
    x0[704] = 183.85876335854499;
    x0[705] = 15.0465302765398;
    x0[706] = 0.038797592394246197;
    x0[707] = 0.50736991301765599;
    x0[708] = 0.74258680630715102;
    x0[709] = 0.19674205378457801;
    x0[710] = 21.100495301256501;
    x0[713] = 342.39787217471599;
    x0[715] = 41.258151573623998;
    x0[716] = 33.177199154586397;
    x0[717] = 202.80804465529101;
    x0[718] = 0.178958476355144;
    x0[719] = 0.70350803123920702;
    x0[720] = 66.204435900083297;
    x0[723] = 26.0958111973762;
    x0[725] = 0.36029721737467302;
    x0[727] = 0.00043244954605721302;
    x0[728] = 0.019705690191402499;
    x0[729] = 0.21883715425654299;
    x0[731] = 0.050265661825492303;
    x0[733] = 0.60138810910575002;
    x0[734] = 0.140675058855874;
    x0[735] = 5.5495911780160103;
    x0[736] = 4.4626306491461998;
    x0[737] = 21.8570345497739;
    x0[743] = 0.0014309433957557;
    x0[744] = 0.00099998774673656097;
    x0[745] = 0.00090228966363420501;
    x0[746] = 0.15055332496794799;
    x0[747] = 0.0022850435270359401;
    x0[748] = 0.0514479737547099;
    x0[752] = 0.00073270361454721099;
    x0[757] = 0.039526596541679199;
    x0[758] = 0.0001562253907147;
    x0[759] = 0.00069836693462002402;
    x0[760] = 0.00016336008002677799;
    x0[761] = 112.157231783874;
    x0[762] = 4.0287866285269596;
    x0[763] = 0.345288121403656;
    x0[764] = 8.4893544673468906;
    x0[765] = 1.9858065386730299;
    x0[766] = 5.3665271818731197;
    x0[771] = 59.457614311451998;
    x0[772] = 66.684900658607702;
    x0[778] = 0.0042238153645795598;
    x0[779] = 0.00340662916346592;
    x0[780] = 0.000409858094676183;
    x0[781] = 0.0014990180314545601;
    x0[782] = 0.0026058878859041299;
    x0[783] = 0.00537621656665323;
    x0[784] = 0.00537621656665323;
    x0[785] = 0.00537621656665323;
    x0[786] = 0.00537621656665323;
    x0[787] = 0.0220867973242166;
    x0[788] = 3.7949823581127899e-6;
    x0[789] = 0.011998469222233301;
    x0[790] = 0.00537621656665323;
    x0[791] = 0.00537621656665323;
    x0[792] = 0.00537621656665323;
    x0[793] = 0.00537621656665323;
    x0[794] = 0.00537621656665323;
    x0[795] = 0.00537621656665323;
    x0[796] = 0.00537621656665323;
    x0[797] = 0.00537621656665323;
    x0[798] = 0.00537621656665323;
    x0[799] = 0.00537621656665323;
    x0[800] = 0.00537621656665323;
    x0[801] = 0.00537621656665323;
    x0[802] = 0.00537621656665323;
    x0[803] = 0.00057177734195566201;
    x0[804] = 0.0094204112069553696;
    x0[805] = 0.0047032481358211297;
    x0[806] = 0.0149876503263069;
    x0[807] = 0.00027703371214223498;
    x0[808] = 0.000208724029696205;
    x0[809] = 0.0013611336724431301;
    x0[810] = 0.00150660799617079;
    x0[811] = 0.0046311434710169797;
    x0[812] = 0.0030815256747875998;
    x0[813] = 0.00018468914142815701;
    x0[814] = 0.0050789513892742898;
    x0[815] = 0.0027311223037218501;
    x0[816] = 0.0037127577403536899;
    x0[817] = 0.0010815699720621501;
    x0[818] = 0.0018962261849370401;
    x0[819] = 0.0244902861510213;
    x0[820] = 0.0047158980770148403;
    x0[821] = 7.5899647162256103e-6;
    x0[822] = 0.0054888094839504803;
    x0[823] = 0.014709351620045301;
    x0[824] = 0.0082110768288367608;
    x0[825] = 0.041061709114780701;
    x0[826] = 0.0050131716950670103;
    x0[827] = 0.00037443825933379499;
    x0[828] = 0.00050473265362900295;
    x0[829] = 0.0099517087370911504;
    x0[830] = 0.00067930184210219201;
    x0[831] = 0.00040859310055681202;
    x0[832] = 0.00178490670243239;
    x0[833] = 0.00246926852101207;
    x0[834] = 2.0239905909934999e-5;
    x0[835] = 0.0040302712643158099;
    x0[836] = 0.39695515465859899;
    x0[837] = 0.0133558079123184;
    x0[838] = 1.2649941193709201e-6;
    x0[840] = 0.0038253422169777099;
    x0[841] = 2.5299882387418698e-6;
    x0[842] = 1.2649941193709201e-6;
    x0[843] = 3.28898471036443e-5;
    x0[844] = 0.00097151548367687996;
    x0[845] = 9.6139553072191204e-5;
    x0[846] = 0.0010954849073752299;
    x0[847] = 0.0054331497426981703;
    x0[848] = 0.0033939792222722098;
    x0[849] = 0.00015053430020514199;
    x0[850] = 0.0011789745192537101;
    x0[851] = 0.00239210387973044;
    x0[852] = 2.65648765067897e-5;
    x0[853] = 0.000189749117905641;
    x0[854] = 0.00048828773007718004;
    x0[855] = 0.00209989023815576;
    x0[856] = 0.0049448620126209897;
    x0[857] = 0.0031928451572922202;
    x0[858] = 0.00220741473830229;
    x0[859] = 5.0599764774837396e-6;
    x0[860] = 6.9574676565401503e-5;
    x0[861] = 0.0026425727153658899;
    x0[862] = 0.0017178620141057299;
    x0[863] = 0.0138086758070532;
    x0[864] = 0.00161160250807858;
    x0[865] = 0.0063009357085866396;
    x0[866] = 0.0016976221081958;
    x0[867] = 0.0027880470390935501;
    x0[868] = 0.0070574021919704603;
    x0[869] = 0.000254263817993558;
    x0[870] = 0.0050599764774837399;
    x0[871] = 0.00153949784327443;
    x0[872] = 0.00135227871360753;
    x0[873] = 0.0048196275948032604;
    x0[874] = 0.00066791689502785398;
    x0[875] = 0.00076532144221941605;
    x0[876] = 0.00163690239046599;
    x0[877] = 0.0121705084224678;
    x0[878] = 0.0092192771419753808;
    x0[879] = 0.0013965535077855101;
    x0[880] = 0.026752095636456599;
    x0[881] = 0.00701692238015057;
    x0[882] = 0.0044755491943343504;
    x0[883] = 0.00078809133636809396;
    x0[884] = 0.00346734888119573;
    x0[885] = 0.0044704892178568803;
    x0[886] = 0.00117770952513434;
    x0[887] = 0.0062136511143500396;
    x0[888] = 0.0010031403366611499;
    x0[889] = 1.2649941193709201e-6;
    x0[890] = 1.01199529549675e-5;
    x0[891] = 0.0194429596147313;
    x0[892] = 0.017183680117534801;
    x0[893] = 0.0059922771434601299;
    x0[894] = 0.00816047706406192;
    x0[895] = 0.00074508153630948096;
    x0[896] = 0.018642218337169499;
    x0[897] = 0.0019480909438312399;
    x0[898] = 0.0015799776550943;
    x0[899] = 0.0088524288473578305;
    x0[900] = 0.0030878506453844599;
    x0[901] = 0.0035356585636417501;
    x0[902] = 0.0035015037224187299;
    x0[903] = 0.0027374472743187102;
    x0[904] = 2.90948647455316e-5;
    x0[905] = 0.0022618094854352398;
    x0[906] = 0.000165714229637593;
    x0[907] = 5.0599764774837396e-6;
    x0[908] = 0.0199287173565698;
    x0[909] = 0.000220108976770543;
    x0[910] = 6.3249705968546902e-6;
    x0[911] = 5.0599764774837396e-6;
    x0[912] = 4.9334770655466402e-5;
    x0[913] = 0.00063882203028232097;
    x0[914] = 0.0046336734592557398;
    x0[915] = 0.00049334770655466405;
    x0[916] = 0.000803271265800546;
    x0[917] = 0.000336488435752669;
}

} // namespace model_SPARCED_au565
} // namespace amici
