#!/usr/bin/env python3

import sys

params = ['']
speciesString = ""
ratelawString = ""

if len(sys.argv) == 3:
    speciesString = sys.argv[1]
    ratelawString = sys.argv[2]
else:
    print("Error in parsing species and ratelaw configurations")

if speciesString == "None" :
    print("Warning: no non-default species values specified")
else:
    products = speciesString.split(",")
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



if ratelawString == "None" :
    print("Warning: no non-default ratelaw values specified")
else:
    for i in range(0,len(params)):
        params[i] = str(params[i])[::-1].replace(",","",1)[::-1] + "\n"
    products = ratelawString.split(",")
    for product in products:
        if product.count(':')%2!=0:
            print("MalformedConfigError: Missing \':\' in params.ratelawVals")
            print(product)
            exit(1)
        if '-' not in product:
            for i in range(0,len(params)):
                params[i]+=str(product+',')
        else:
            productID, value = tuple(product.rsplit(":",1))
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



print(params)

for idx, entry in enumerate(params):
    with open(str("sweep"+str(idx)+".txt"),"w") as f:
        if entry != "":
            f.write(entry[:-1])
        else:
            f.write("")
