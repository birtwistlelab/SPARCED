#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCED_u87i_x_rdata.h"

namespace amici {
namespace model_SPARCED_u87i {

void total_cl_SPARCED_u87i(realtype *total_cl, const realtype *x_rdata){
    total_cl[0] = m_INS;
    total_cl[1] = m_INSR;
    total_cl[2] = m_MGMT;
    total_cl[3] = m_BRCA2;
    total_cl[4] = m_MSH6;
    total_cl[5] = m_IGF1R;
    total_cl[6] = m_IGF2;
    total_cl[7] = m_IGF1;
    total_cl[8] = m_IRS2;
    total_cl[9] = m_IRS1;
    total_cl[10] = m_EIF4E;
    total_cl[11] = m_FGF2;
    total_cl[12] = m_FGF1;
    total_cl[13] = m_FGFR2;
    total_cl[14] = m_FGFR1;
    total_cl[15] = m_MDM4;
    total_cl[16] = m_CDKN2A;
    total_cl[17] = m_SOS1;
    total_cl[18] = m_EIF4EBP1;
    total_cl[19] = m_RPS6KB2;
    total_cl[20] = m_RPS6KB1;
    total_cl[21] = m_RPTOR;
    total_cl[22] = m_RHEB;
    total_cl[23] = m_FOXO3;
    total_cl[24] = m_MAPK3;
    total_cl[25] = m_MAPK1;
    total_cl[26] = m_PEBP1;
    total_cl[27] = m_PRKCD;
    total_cl[28] = m_PRKCG;
    total_cl[29] = m_PRKCB;
    total_cl[30] = m_PRKCA;
    total_cl[31] = m_TSC2;
    total_cl[32] = m_TSC1;
    total_cl[33] = m_GSK3B;
    total_cl[34] = m_MTOR;
    total_cl[35] = m_RICTOR;
    total_cl[36] = m_PDPK1;
    total_cl[37] = m_AKT2;
    total_cl[38] = m_AKT1;
    total_cl[39] = m_PTEN;
    total_cl[40] = m_CTNNB1;
    total_cl[41] = m_MYC;
    total_cl[42] = m_JUN;
    total_cl[43] = m_FOS;
    total_cl[44] = m_DUSP1;
    total_cl[45] = m_RPS6KA4;
    total_cl[46] = m_RPS6KA3;
    total_cl[47] = m_RPS6KA2;
    total_cl[48] = m_RPS6KA1;
    total_cl[49] = m_DUSP6;
    total_cl[50] = m_MAP2K2;
    total_cl[51] = m_MAP2K1;
    total_cl[52] = m_BRAF;
    total_cl[53] = m_RAF1;
    total_cl[54] = m_NF1;
    total_cl[55] = m_HRAS;
    total_cl[56] = m_KRAS;
    total_cl[57] = m_NRAS;
    total_cl[58] = m_RASGRP3;
    total_cl[59] = m_RASGRP1;
    total_cl[60] = m_PIK3C2A;
    total_cl[61] = m_PIK3R4;
    total_cl[62] = m_PIK3R3;
    total_cl[63] = m_PIK3R2;
    total_cl[64] = m_PIK3R1;
    total_cl[65] = m_PIK3CD;
    total_cl[66] = m_PIK3CG;
    total_cl[67] = m_PIK3CB;
    total_cl[68] = m_PIK3CA;
    total_cl[69] = m_PLCG2;
    total_cl[70] = m_PLCG1;
    total_cl[71] = m_GRB2;
    total_cl[72] = m_CBL;
    total_cl[73] = m_SPRY2;
    total_cl[74] = m_PDGFB;
    total_cl[75] = m_PDGFRB;
    total_cl[76] = m_PDGFRA;
    total_cl[77] = m_HGF;
    total_cl[78] = m_MET;
    total_cl[79] = m_EGFRvIII;
    total_cl[80] = m_ERBB4;
    total_cl[81] = m_ERBB3;
    total_cl[82] = m_ERBB2;
    total_cl[83] = m_EGFR;
    total_cl[84] = m_NRG1;
    total_cl[85] = m_EGF;
    total_cl[86] = m_BCL2L11;
    total_cl[87] = m_PMAIP1;
    total_cl[88] = m_BBC3;
    total_cl[89] = m_BAD;
    total_cl[90] = m_CASP9;
    total_cl[91] = m_APAF1;
    total_cl[92] = m_DIABLO;
    total_cl[93] = m_CYCS;
    total_cl[94] = m_BAX;
    total_cl[95] = m_MCL1;
    total_cl[96] = m_BCL2L1;
    total_cl[97] = m_BCL2;
    total_cl[98] = m_BID;
    total_cl[99] = m_PARP1;
    total_cl[100] = m_XIAP;
    total_cl[101] = m_CASP6;
    total_cl[102] = m_CASP7;
    total_cl[103] = m_CASP3;
    total_cl[104] = m_BFAR;
    total_cl[105] = m_CASP10;
    total_cl[106] = m_CASP8;
    total_cl[107] = m_CFLAR;
    total_cl[108] = m_TNFRSF10B;
    total_cl[109] = m_TNFRSF10A;
    total_cl[110] = m_TNFSF10;
    total_cl[111] = m_CDK6;
    total_cl[112] = m_CDK4;
    total_cl[113] = m_CDK2;
    total_cl[114] = m_CDK1;
    total_cl[115] = m_CDKN1A;
    total_cl[116] = m_CHEK1;
    total_cl[117] = m_WEE1;
    total_cl[118] = m_CDC20;
    total_cl[119] = m_CCNB1;
    total_cl[120] = m_CDH1;
    total_cl[121] = m_CDKN1B;
    total_cl[122] = m_CCNA2;
    total_cl[123] = m_CDC25C;
    total_cl[124] = m_CDC25B;
    total_cl[125] = m_CDC25A;
    total_cl[126] = m_SKP2;
    total_cl[127] = m_CCNE2;
    total_cl[128] = m_CCNE1;
    total_cl[129] = m_CCND3;
    total_cl[130] = m_CCND2;
    total_cl[131] = m_CCND1;
    total_cl[132] = m_E2F3;
    total_cl[133] = m_E2F2;
    total_cl[134] = m_E2F1;
    total_cl[135] = m_RB1;
    total_cl[136] = m_ATR;
    total_cl[137] = m_ATM;
    total_cl[138] = m_PPM1D;
    total_cl[139] = m_MDM2;
    total_cl[140] = m_TP53;
}

} // namespace model_SPARCED_u87i
} // namespace amici
