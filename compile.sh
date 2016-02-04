for ccfile in ./*.cc
do
    ofile=`echo $ccfile |sed 's/\.cc$//g'`
    echo "compiling : " $ccfile ", executible name: " $ofile
    g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent74 -l TopTreeAna74 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
done
