#!/bin/sh#
cat $1 | sed '/^[a-z]/ d' | sed '/^\	/ d' | sed '/^\-/ d' |
	sed '/^Black\-Scholes/ d' |
	sed '/^Trials/ d' |
	sed '/^Confidence/ d' |
	sed '/^Average/ d' |
	sed '/^Arverage/ d' |
	sed '/^Standard Deviation    / d' |
	sed 's/(.*)//g' |
	if [ $# -eq 1 ]; then
			  sed 's/Total[a-zA-Z ]*/Total/g' | 
			  sed 's/PRNG[a-zA-Z ]*/PRNG/g' |
			  sed 's/BS[a-zA-Z ]*/BS/g'  |
			  sed 's/Black[a-zA-Z ]*/BS/g' |
			  sed 's/Standard[a-zA-Z ]*/Stddev/g' 
	else 
			  # print total only
			  if [ $2 -eq 0 ] ; then
						 sed '/^PRNG/ d' |
						 sed '/^BS/ d' |
						 sed '/^Black/ d' |
						 sed '/^Standard/ d' |
						 sed 's/Total[a-zA-Z ]*/Total/g'
			  fi
			  # print PRNG only
			  if [ $2 -eq 1 ] ; then
						 sed '/^Total/ d' |
						 sed '/^BS/ d' |
						 sed '/^Black/ d' |
						 sed '/^Standard/ d' |
						 sed 's/PRNG[a-zA-Z ]*/SPAWN/g'
			  fi
			  #print BS only
			  if [ $2 -eq 2 ] ; then
						 sed '/^Total/ d' |
						 sed '/^PRNG/ d' |
						 sed '/^Standard/ d' |
						 sed 's/BS[a-zA-Z ]*/BS/g' |
						 sed 's/Black[a-zA-Z ]*/BS/g'
			  fi
			  #print std dev
			  if [ $2 -eq 3 ] ; then
						 sed '/^Total/ d' |
						 sed '/^PRNG/ d' |
						 sed '/^Black/ d' |
						 sed 's/Standard[a-zA-Z ]*/Stddev/g'
			  fi
	fi
