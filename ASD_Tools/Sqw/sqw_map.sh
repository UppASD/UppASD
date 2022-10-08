#!/bin/bash
mypath=.
iname=`ls -l sqw.*.out* | tail -1 | awk '{ print $NF }'`
if [[ "$iname" == ../sqw.*.out.gz ]]; then
    zcat $iname > sqw.dat
else
    cat $iname > sqw.dat
fi
dim1=`tail -1 sqw.dat | awk '{ print $1 }' ` 
dim2=`tail -1 sqw.dat | awk '{ print $5 }' ` 
dim2half=`tail -1 sqw.dat | awk '{ print $5/2 }' ` 
echo $dim1 $dim2 > tmp1
time $mypath/decodeS.x < tmp1 
awk '{ print $3 }' sqw.norm.dat > sqw.arr
$mypath/col2arr.awk sqw.arr $dim2half | $mypath/flip.awk > sqw.mat
echo $dim2
dt=`grep timestep inpsd.dat | tail -1 | awk '{ print $2}' ` 
ss=`grep sc_step inpsd.dat | tail -1 | awk '{ print $2}'` 
ns=`grep sc_nstep inpsd.dat | tail -1 | awk '{ print $2}'` 
# energy in meV below
yscale=`echo $dt $ss $ns| awk '{ print 4.13566733e-12/$1/$2/$3}'`
sed "s/YSCALE/$yscale/" $mypath/base.gnu > $mypath/plot.gnu
gnuplot $mypath/plot.gnu 
rm -f sqw.dat tmp1 sqw.arr
