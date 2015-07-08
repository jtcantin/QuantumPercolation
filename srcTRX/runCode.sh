cp * /tmp/jtc/QPerc/run4/
cd /tmp/jtc/QPerc/run4

echo $PWD

export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

echo $OMP_NUM_THREADS
echo $MKL_NUM_THREADS

echo "----------------"
echo "N30 q60 w13 Dis0"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 0

echo "----------------"
echo "N30 q60 w13 Dis1"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 1

echo "----------------"
echo "N30 q60 w13 Dis2"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 2

echo "----------------"
echo "N30 q60 w13 Dis3"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 3

echo "----------------"
echo "N30 q60 w13 Dis4"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 4

echo "----------------"
echo "N30 q60 w13 Dis5"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 5

echo "----------------"
echo "N30 q60 w13 Dis6"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 6

echo "----------------"
echo "N30 q60 w13 Dis7"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 7

echo "----------------"
echo "N30 q60 w13 Dis8"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 8

echo "----------------"
echo "N30 q60 w13 Dis9"
/usr/bin/time ./QmPerc_LINUX_TRX dis_in_N30_q60.txt 9

mv log4.txt log4_old.txt
cd /home/jtcantin/code/QuantumPercolation/srcTRX/runDir/run4
cp log4.txt /tmp/jtc/QPerc/run4/
cp /tmp/jtc/QPerc/run4/* /home/jtcantin/code/QuantumPercolation/srcTRX/data/files_n030q60w013_2015-07-06_18-38-19_TB/

