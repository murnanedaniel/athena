#! /usr/bin/env python
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# Author : Benjamin Trocme (CNRS/IN2P3 - LPSC Grenoble)- 2017 - 2022
# Python 3 migration by Miaoran Lu (University of Iowa)- 2022
#
# For each year/tag/system, an output directory YearStats-[system]/[year]/[tag] is required.
# This directory also contains some subdirectories
# This script creates all needed directories for a new set of system/year/tag.
# If the set already exists, the script simply exits
#
# Documentation: https://twiki.cern.ch/twiki/bin/viewauth/Atlas/DataQualityDemo
#############################################################################################

from __future__ import print_function
import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-y','--year',dest='parser_year',default = "2017",help='Year [Default: 2017]. May also include special conditions such as 5TeV, 5TeV... Check the runList.py file',action='store')
parser.add_argument('-t','--tag',dest='parser_tag',default = "Tier0_2017",help='Defect tag [Default: "Tier0_2017"]',action='store')
parser.add_argument('-s','--system',dest='parser_system',default="LAr",help='System: LAr, CaloCP [Default: "LAr"]',action='store')

args = parser.parse_args()
parser.print_help()

direct = "YearStats-%s"%args.parser_system
if not os.path.exists(direct):
  print("%s system directory does not exists. Creating it"%direct)
  os.system("mkdir %s"%direct)

direct = "YearStats-%s/%s"%(args.parser_system,args.parser_year)
if not os.path.exists(direct):
  print("%s year directory does not exists. Creating it"%direct)
  os.system("mkdir %s"%direct)

direct = "YearStats-%s/%s/%s"%(args.parser_system,args.parser_year,args.parser_tag)
if not os.path.exists(direct):
  print("%s tag directory does not exists. Creating it"%direct)
  os.system("mkdir %s"%direct)
  os.system("mkdir %s/Run"%direct)
  os.system("mkdir %s/Weekly"%direct)
