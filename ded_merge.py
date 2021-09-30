#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge distributed analysis sets from a FileHandler.

Usage:
    ded_merge.py d1 d2 d3 ...

"""
import os
import sys
from dedalus.tools import post

n=len(sys.argv) - 1
k=0
while (k<n):
    k=k+1
    name =  sys.argv[k]
    if os.path.isdir(name):
        print('post.merge_process_files "{}"'.format(name))
        #post.merge_process_files(name, cleanup=True)
        #post.merge_analysis(name, cleanup=True)
        post.merge_process_files_single_set(name, cleanup=True)
#    else:
#        print('"{}" is not a directory'.format(name))
    
