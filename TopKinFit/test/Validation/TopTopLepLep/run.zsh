#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../:../../IPHCNtuple/NtupleProducer/:../../IPHCNtuple/NtupleProducer/obj

#./test list.txt output ttH
#./test list.txt output ttbar | egrep -v V | egrep -v MINIMUM | egrep -v FUNCTION | egrep -v "====" | egrep -v '\\n'
#./test list.txt output ttbar 100000

./gen list.txt output
