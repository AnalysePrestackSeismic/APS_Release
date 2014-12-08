#!/bin/sh
## ------------------ Disclaimer  ------------------
# 
# BG Group plc or any of its respective subsidiaries, affiliates and 
# associated companies (or by any of their respective officers, employees 
# or agents) makes no representation or warranty, express or implied, in 
# respect to the quality, accuracy or usefulness of this repository. The code
# is this repository is supplied with the explicit understanding and 
# agreement of recipient that any action taken or expenditure made by 
# recipient based on its examination, evaluation, interpretation or use is 
# at its own risk and responsibility.
# 
# No representation or warranty, express or implied, is or will be made in 
# relation to the accuracy or completeness of the information in this 
# repository and no responsibility or liability is or will be accepted by 
# BG Group plc or any of its respective subsidiaries, affiliates and 
# associated companies (or by any of their respective officers, employees 
# or agents) in relation to it.
## ------------------ License  ------------------ 
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
## github
# https://github.com/AnalysePrestackSeismic/
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
umask 002
rm -rf /localcache/mcr_cache_umask_friendly/$1

exit
