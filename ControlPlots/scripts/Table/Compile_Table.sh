#echo -e "\033[40;32m----------------------------------------------------------"
#echo -n "Have you defined Signal.txt, Background.txt and Variables.txt (y/n)? "
#read -n 1 decision
#echo
#
#if [ "$decision" != "y" ]; then
#	echo -e "\033[40;36mPlease update your signal, background and viarbles txt files"
#	echo "Aborting compilation"
#	echo -e "\033[40;32m----------------------------------------------------------\033[37m"
#	exit 1
#fi
#
#
if [ -d Output_Table ]; then
	:
else
	mkdir Output_Table
fi
echo -e "-- Compiling macro to analyze variables tables"
if [ -s Table.o ]; then
	rm Table.o
fi
g++ -c `root-config --cflags` Table.cpp
g++ -o Table `root-config --libs` Table.o
if [ -s Table.o ]; then

echo -e "COMPILATION FINISHED"
echo "Flags for running ./Table:"
./Table --help
echo -e "----------------------------------------------------------"
else
	exit 1
fi
