#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
umask 002
rm -rf /localcache/mcr_cach*;

exit