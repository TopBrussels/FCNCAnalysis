#if there is one arg compile only this .cc
#export PATH=$PATH:$ROOTSYS/bin:/user/kderoove/Software/MadAnalysis5/v1.1.12beta/bin/:/user/kderoove/Software/LHAPDF/local/bin/:/user/kderoove/Software/LHAPDF/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib:/user/kderoove/Software/LHAPDF/local/lib

if [[ -n $1 ]]
then
    ccfile=$1

    ofile=`echo $ccfile |sed 's/\.cc$//g'`
    echo "compiling : " $ccfile ", executible name: " $ofile
    g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent80 -l TopTreeAna80 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
#    g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent80 -l TopTreeAna80 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile `lhapdf-config --cflags --ldflags`
#    cp ~/lib/libTopTreeAnaContent74.so /localgrid/qpython/lib/
#    cp ~/lib/libTopTreeAna74.so /localgrid/qpython/lib/
    

# if there is no arg compile all .cc
else
    for ccfile in ./*.cc
    do
	ofile=`echo $ccfile |sed 's/\.cc$//g'`
	echo "compiling : " $ccfile ", executible name: " $ofile
	g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent80 -l TopTreeAna80 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile

    done
fi
