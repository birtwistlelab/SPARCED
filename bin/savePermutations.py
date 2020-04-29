#!/usr/bin/env python3

import sys

params = ['']
sweepString = ""

if len(sys.argv) < 2:
    sweepString = ""
else:
    sweepString = sys.argv[1]

if sweepString == None or sweepString == "":
    print("Warning: no non-default values specified")
else:
    products = sweepString.split(",")
    for product in products:
        if ':' not in product:
            print("MalformedConfigError: Missing \':\' in params.speciesVals")
            print(product)
            exit(1)

        if '-' not in product:
            for i in range(0,len(params)):
                params[i]+=str(product+',')
        else:
            productID, value = tuple(product.split(":"))
            lower_bound, rest = tuple(value.split("-"))
            upper_bound, increment = tuple(rest.replace(")","").split("("))
            valueRange = []
            if round(((float(upper_bound) - float(lower_bound)) / float(increment)),3) % 1.0 != 0.0:
                print("MalformedConfigError: invalid increment for bounds")
                print(product)
                exit(1)
            else:
                current = float(lower_bound)
                while round(current,7) <= float(upper_bound):
                    valueRange.append(productID + ":" + str(round(current,6)))
                    current+=round(float(increment),7)
                originalLen = len(params)
                params = params * len(valueRange)
                for idx, value in enumerate(valueRange):
                    for i in range(0,originalLen):
                        params[(idx*originalLen)+i]+=str(value+",")

for idx, entry in enumerate(params):
    with open(str("sweep"+str(idx)+".txt"),"w") as f:
        if entry != "":
            f.write(entry[:-1])
        else:
            f.write("")
