#!/bin/bash

cat << "EOF"
Created by: Haiyang (Skyler) Chang (｡•̀ᴗ-)✧
Method: Branch
Date: 04/01/18
EOF
####checking the input file##########

echo -e "Please input sequences file in Fasta Format"
if [ $# -ne 1 ]
then
echo -e "Error input: Just one file with many fasta format sequences!"
exit 1
fi

###################  find the best 2 anchor sets##################
function Findanchorsets(){
awk '/.fitness/ {close(x); x="anchor_group"++i;}{print > x;}' best.pcld
for i in anchor_group*
do awk 'NR==FNR{a[$0]=1;next}!a[$0]' $i CLU.txt > CLU$i.txt
done
for i in anchor_group*
do tail -n +2 $i > $i.tmp && mv $i.tmp $i.txt
done
mkdir result
for i in {1..30}
do perl index_anchor.pl anchor_group$i.txt CLUanchor_group$i.txt > result/Anchor_Set$i.txt
done
id=($(Rscript Check_Ratio.R))
}


################Find seq's coors############
function Findseqcoors(){
id=${id[0]}
perl an_seq_amb.pl anchor_group$id CLUanchor_group$id.txt > Clusters.txt
sed 1d anchor_group$id | sed 's/^ //g' >> anchors$id.txt
rm CLUanchor_group*.txt
rm anchor_group*
}

############ cut Clusters by their coors########
############ place in dir include Clusters.txt##########
function Cutclusters(){
awk '/>/ {close(x); x="cluster"++i;}{print > x;}' Clusters.txt
for i in cluster*
do
sed 1d $i | sed 's/^ //g' >> $i.txt
awk '{print;} FNR>1 { print ">seq";}' $i | sed -e '$ d' | sed 's/^ //g' >> $i.fasta
rm $i
done
for t in cluster*.*; do mkdir -p "${t%.*}"; mv "$t" "${t%.*}"; done
}

########place in dir include all cluster#########
function Copyscripts(){
for k in cluster*/
do
mkdir -p $k/anchordis${ancdis}/Conwayreportanchordis${ancdis}
cp ~/Desktop/shell_script/a.out ~/Desktop/shell_script/SelectDNAII.cpp ~/Desktop/shell_script/stat.h ~/Desktop/shell_script/stat.cpp $k/anchordis${ancdis}
cp ~/Desktop/shell_script/an_seq_amb.pl $k/anchordis${ancdis}
cp ~/Desktop/shell_script/index_anchor.pl ~/Desktop/shell_script/Check_Ratio.R ~/Desktop/shell_script/check_HD.R $k/anchordis${ancdis}
done
}

##########test if more than 100 seqs########
function Clustersize(){
num=$(wc -l cluster*.txt | awk '{print $1}')
echo "num is $num"
if [ $num -gt "100" ]
then
./a.out cluster*.fasta 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
mv *.txt CLU.txt
Findanchorsets
Findseqcoors
rm cluster*.fasta
Cutclusters
ancdis=`expr $ancdis - 3`
Copyscripts
ancdis=`expr $ancdis + 3`
else
mv cluster*.txt end_cluster.txt
Rscript check_HD.R end_cluster.txt
echo "******the end ********"
fi
}

######Find anchors###########
echo -e "*******Prepare the Conway Algorithm"
mkdir -p Output2
g++ -lm -O3 SelectDNAII.cpp stat.cpp
seqcount=$(grep ">" $1 | wc -l)
seqlength=$(awk 'NR==2 {print length($0)}' $1)
ancdis=`expr $seqlength - 3`
echo "Total $seqcount sequences and length is $seqlength "
echo -e "********Running Conway Crossover Operator********"
mkdir -p Output2/anchordis${ancdis}/Conwayreportanchordis${ancdis}
cp ./a.out SelectDNAII.cpp stat.h stat.cpp $1 Output2/anchordis${ancdis}
cd Output2/anchordis${ancdis}
./a.out $1 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
bdac=$(grep '1' best.pcld | wc -l)
while [ "$bdac" -gt "28" ]
do
echo $bdac
ancdis=`expr $ancdis - 3`
echo $ancdis
pwd
cd ../
mkdir -p anchordis${ancdis}/Conwayreportanchordis${ancdis}
pwd
cp ../a.out ../SelectDNAII.cpp ../stat.h ../stat.cpp ../$1 anchordis${ancdis}
cp ../index_anchor.pl ../Check_Ratio.R anchordis${ancdis}
cp ../an_seq_amb.pl anchordis${ancdis}
cd anchordis${ancdis}
grep -v '^>' $1 > CLU.txt
./a.out $1 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
bdac=$(grep '1' best.pcld | wc -l)
done
echo "Start anchor distance is $ancdis"
Findanchorsets
Findseqcoors
Cutclusters
ancdis=`expr $ancdis - 3`
Copyscripts
pwd ##### in Output/anchordis39/######
echo "ancdis is $ancdis"      #36


for k in cluster*/
do
mv $k/cluster*.* $k/anchordis${ancdis}/
cd $k
cd anchordis${ancdis}
num=$(wc -l cluster*.txt | awk '{print $1}')
echo "num is $num"
if [ $num -gt "100" ]
then
./a.out cluster*.fasta 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
mv *.txt CLU.txt
Findanchorsets
Findseqcoors
rm cluster*.fasta
Cutclusters
echo "After cut in "
pwd            #/Output2/anchordis39/cluster3/anchordis36
ancdis=`expr $ancdis - 3`
Copyscripts
dir="$PWD"
for i in $dir/cluster*/
do
cd $i
mv cluster*.* anchordis*/   ###Output2/anchordis39/cluster1/anchordis36/cluster1
cd anchordis*/
Clustersize
echo " after loop in"
pwd
done
ancdis=`expr $ancdis + 3`
cd ../../
else
mv cluster*.txt end_cluster.txt
Rscript check_HD.R end_cluster.txt
echo "******the end ********"
fi
cd ../../
done

echo "in the end in "
pwd








