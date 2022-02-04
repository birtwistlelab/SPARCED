#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 14:51:46 2022

@author: arnab
"""
from mpi4py.futures import MPICommExecutor

def f(x,y):
    print(f"{x} + {y} = {x + y}")
    
    results = {}
    
    results['xval'] = x
    results['yval'] = y
    results['sum'] = x + y
    
    return results

xvals = [1,2,3]
yvals = [10,20,30]


#%%
# print("Using submit!")
# futures = []
# with MPICommExecutor() as executor:
#     if executor is not None:
#         for x,y in zip(xvals, yvals):
#             fut = executor.submit(f,x,y)
#             futures.append(fut)
#         for fut in futures:
#             fut.result()

#%%
print("Using map!")
with MPICommExecutor() as executor:
    if executor is not None:
        results = list(executor.map(f,xvals,yvals))
        # results = executor.map(f,xvals,yvals)