#include "sundials/sundials_types.h"

void dJydy_rowvals_BigModel_FULL_v4(sunindextype *rowvals, int index){
    switch(index) {
        case 0:
                rowvals[0] = 0;
            break;
        case 1:
                rowvals[0] = 0;
            break;
        case 2:
                rowvals[0] = 0;
            break;
    }
}