# coding: utf-8
#!/usr/bin/env python
"""Extractoutput module.
Utility functions for selecting data in output files
"""
import csv
import os, os.path
from ast import literal_eval 
from pprint import pprint
#from decimal import Decimal
#gtvl = groundtruthvalue = [5,6,9]
#gtv = groundtruthvalue = 7
#tvl = testvalue = [5,6,8]
#tv = testvalue = 7
def scalarcmp(l,r):
    return l == r
#scalarcmp(gtv,tv)
def iterablecmp(l,r):
    return [le==ri for le,ri in zip(l,r)]
#iterablecmp(gtvl,tvl)
class Skipspace(csv.excel):
    skipinitialspace = True    

def csv_like(filepath, selectors, extractors, headers=None, verbose=1): #single value for now
   # (,pass_lambda_as_selector=False) #not yet implemented
    """selectors is a dictionary with the format {column_header:row_value} if headers==True, otherwise the format {zero_indexed_column_number:row_value}
extractors is a list with values corresponding to the keys portion of the selectors
"""
    if os.path.exists(filepath):
        filepath = os.path.abspath(filepath)
    else:
        print ("%s not found. (Working from %s)" % (filepath, os.getcwd()))
        return None
   
    #with open(filepath, 'rb') as csvfile:
    with open(filepath, 'rU') as csvfile:
        linereader = csv.reader(csvfile, 
                                delimiter=' ', 
                                quotechar='|', 
                                dialect=Skipspace)

        if not isinstance(extractors, list):
                extractors = [extractors]

        if headers:
            header_number_dict = {s:i for s,i in zip(headers, range(len(headers)))}
            number_header_dict = {v:k for k,v in header_number_dict.items()}
            selectors = new_selectors = {header_number_dict[k]:selectors[k] for k in selectors.keys()}
            
            
            extractors_new2old = {headers.index(ev):ev for ev in extractors}
            extractors = new_extractors = [headers.index(ev) for ev in extractors]


        line = 0
        li = True
        while li:
            # for proper comparisons we'll need to eval the contents
            # add a step with literal_eval if this proves needed
            # (e.g. Do we want to do extractors=[dim2: ">= 250, <= 600"])
            # passing functions would probably be a good solution
            try:
                li = next(linereader)
                #li = linereader.next()
                
            except StopIteration as SI:
                print ("stopped iteration trying to get to line ", line+1, "\n", SI)
                break
            line += 1

            try:
                if all( [literal_eval(li[k]) == v for k,v in selectors.items()] ):
                    # added float conversion. Take care to change if comparing strings or other types
                    return [ float(li[x]) for x in extractors ]
                    #return [ Decimal(li[x]) for x in extractors ]
                    #return [ li[x] for x in extractors ]
                        
            except SyntaxError as se:
                if verbose >=1:
                    print ("SyntaxError occured while processing line %s of %s. Presumably literal_eval of a non-number\n" % (line,filepath) , se)

    if verbose >=0:
        print ("Nothing extracted")
    return None
if __name__ == "__main__":
    pass
