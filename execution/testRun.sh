#!/bin/bash
echo " "
echo " "
echo "----------Beginning Shell Script-----------"
origPath=$(pwd)
echo Initial Directory: $origPath

if [ ! -d /tmp/jtc ];
  then
	mkdir /tmp/jtc
  fi

dirID=testRuns_
dirdate=$(date +%Y-%m-%d_%H-%M-%S-%N)
dirname=$dirID$dirdate 
echo "Execution Directory Name: $dirname"

mkdir /tmp/jtc/$dirname

# Copy Desired Files
copyList=(../qmPerc_LINUX)
echo Files Copied: ${copyList[*]}
cp -rf "${copyList[@]}" /tmp/jtc/$dirname  

echo "-----------Preparation Complete------------"
echo "Executing..."
cd /tmp/jtc/$dirname
{ time ./qmPerc_LINUX 100 1 1 0 1 ; } &> log_$dirID$dirdate.txt
cd $origPath

echo "------------Execution Complete-------------"
echo "Copying..."

# Copy output directory to work directory
mkdir ~/code/QuantumPercolation/work/$dirname
cp -rf /tmp/jtc/$dirname/* ~/code/QuantumPercolation/work/$dirname
rm -rf /tmp/jtc/$dirname
echo "-------------Copying Complete--------------"

echo "-------------------Done--------------------"
echo " "
echo " "
