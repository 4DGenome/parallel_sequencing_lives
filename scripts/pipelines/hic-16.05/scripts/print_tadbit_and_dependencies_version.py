#!/usr/bin/python


#==================================================================================================
# Created on: 2016-04-05
# Usage: ./print_tadbit_and_dependencies_version.py
# Author: Javier Quilez (GitHub: jaquol)
# Goal: print the version of TADbit and its dependencies
#==================================================================================================

# import packages
import pytadbit
import re

# print TADbit and dependencies versions
print re.sub(r"\n+", ";", pytadbit.get_dependencies_version()).replace(" ", "")