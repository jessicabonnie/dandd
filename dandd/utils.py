#! /usr/bin/env python3
from __future__ import annotations
import sys
import hashlib
import csv
from typing import List, Dict
from math import factorial
from itertools import permutations
from random import sample, shuffle

def insert_pre_ext(filename, string):
    toks = filename.split('.')
    return '.'.join(toks[:-1] + [string] + [toks[-1]])

def write_listdict_to_csv(outfile: str, listdict:List[Dict], suffix:str="", last_col: str = None):
    '''
    Take a list of dictionaries and write them to a csv files with the keys as headers. If a particular column needs to be last, provide the name.
    '''
    writer = open(outfile+suffix, "w") if outfile is not None and outfile != '-' else sys.stdout
    fieldnames=set()
    for x in listdict:
        fieldnames.update(x.keys())
    fieldnames=list(fieldnames)
    # The fastas field should be at the end since it sometimes has commas
    if "fastas" in fieldnames:
        last_col="fastas"
    if "files" in fieldnames:
        last_col="files"
    if last_col:
        i=fieldnames.index(last_col)
        fieldnames=fieldnames[:i]+fieldnames[i+1:]+ [fieldnames[i]]
    dict_writer = csv.DictWriter(writer, fieldnames=fieldnames)
    dict_writer.writeheader()
    dict_writer.writerows(listdict)
    writer.close()


def blake2b(fname):
    '''Create a blake2b hexsum from a file'''
    hash_blake2b = hashlib.blake2b()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_blake2b.update(chunk)
    return hash_blake2b.hexdigest()


def canon_command(canon:bool, tool='dashing'):
    '''Determine what should be added to sketching command when not canonicalizing
    '''
    outstr=''
    if not canon:
        if tool == 'dashing':
            outstr='--no-canon'
        if tool == 'kmc':
            outstr='-b'
    return outstr

def permute(length, norder, preexist=set(), exhaust=False, verbose=False)-> Set[Tuple[int]]:
    '''Create a set of ordering tuples. One approach will be used if the maximum number of possible permutations is low or desired. Another will be used if not.'''
    fact=factorial(length)
    newset=preexist.copy()
    norder = min(norder, fact)
    if verbose:
        print(f"{norder} permutations will be produced.")
    # if the number of permutations is low enough, generate all of them
    if fact < 5041 or norder==fact or exhaust==True:
        if verbose:
            print(f"Factorial of {norder} less than 5041 (7! + 1). Randomized list of all permutations will be produced and then subset.")
        # create a randomized list of all permutations
        newpermute=list(permutations(range(length)))
        shuffle(newpermute)
    else:
        newpermute=list()
        while len(newpermute) < norder:
            # tmpset=set([tuple(sample(list(range(length)),length)) for _ in range(norder)])
            newpermute.extend([tuple(sample(list(range(length)),length)) for _ in range(norder)])
            newpermute=list(set(newpermute))
    index = 0
    while len(newset) < norder:
        newset.update([newpermute[index]])
        index+=1
    return newset