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
###################  find the best anchor set##################
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
rm -rf result
}
################Find seq's coors############
function Matchancseqs(){
id=${id[0]}
perl an_seq_amb.pl anchor_group$id CLUanchor_group$id.txt > Anc_clu.txt
sed 's/^ //g' Anc_clu.txt > Anc_clus.txt
rm CLUanchor_group*.txt
rm anchor_group*
}
############ cut Clusters by their coors########
############ place in dir include Clusters.txt##########
function Cutclusters(){
awk '/^>/{$0=$0",'$clustername'."(++i)}1' Anc_clus.txt > Anc_label_clus$clustername.txt
grep '^>' Anc_label_clus$clustername.txt | sed "s/$/,ancdis$ancdis/g" > Anchors$clustername.txt
awk '/>/ {close(x); x="cluster'$clustername'."++i;}{print > x;}' Anc_label_clus$clustername.txt
}

########place in dir include all cluster#########
function Copyscripts(){
ancdis=`expr $ancdis - 3`
mkdir -p ../anchordis${ancdis}/Conwayreportanchordis${ancdis}
cp ~/Desktop/shell_script/a.out ~/Desktop/shell_script/SelectDNAII.cpp ~/Desktop/shell_script/stat.h ~/Desktop/shell_script/stat.cpp ../anchordis${ancdis}
cp ~/Desktop/shell_script/an_seq_amb.pl ../anchordis${ancdis}
cp ~/Desktop/shell_script/index_anchor.pl ~/Desktop/shell_script/Check_Ratio.R ~/Desktop/shell_script/check_HD.R ../anchordis${ancdis}
}

######Find  start anchor distance ###########
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
Matchancseqs
awk '/^>/{$0=$0","(++i)}1' Anc_clus.txt > Anc_label_clus.txt
grep '^>' Anc_label_clus.txt | sed "s/$/,ancdis$ancdis/g" > Anchorsmothers.txt
awk '/^>/ {close(x); x="cluster"++i;}{print > x;}' Anc_label_clus.txt
for i in cluster*
do
sed 1d $i | sed 's/^ //g' >> $i.txt
rm $i
done
pwd
#########################Clustering the original input large-scale sequences#################
while :
do
Copyscripts
for k in cluster*.txt
do
num=$(wc -l $k | awk '{print $1}')
echo "num at first is $num"
if [ $num -gt "150" ]
then
cp $k ../anchordis${ancdis}
fi
done
cd ../anchordis${ancdis}
for j in cluster*.txt
do
awk 'BEGIN { print ">seq";}{print;} { print ">seq";}' $j | sed -e '$ d' | sed 's/^ //g' >> run.fasta
./a.out run.fasta 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
cp $j CLU.txt
Findanchorsets
Matchancseqs
echo "$j"
clustername=$(echo "${j%.*}" | sed 's/[^0-9.]*//g')
rm $j
echo "name is $clustername"
Cutclusters
done
for i in cluster*
do
sed 1d $i | sed 's/^ //g' >> $i.txt
rm $i
done
n=$(find . -name "cluster*"  -type f -exec bash -c '[[ $(wc -l < "$1") -gt 150 ]] && echo "$1"' _ '{}' \; | wc -l)
echo "n is $n"
if [ $n -eq 0 ]
then
echo "**********Clustering Work Done*************"
break
fi
done
#############Collect the anchor information######################
pwd
cd ../
find . -type f -name 'Anchors*.txt' -exec cat {} + >> outputanchors.txt
awk -F',' '{print $3,$2,$1}' outputanchors.txt | sort -gk1,1r -gk2,2 > anchorsfile.txt
find . -type f -not -name 'cluster*.txt' -not -name 'anchor*file.txt' -print0 | xargs -0 rm --
