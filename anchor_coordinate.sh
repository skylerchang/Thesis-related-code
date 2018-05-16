#!/bin/bash

cat << "EOF"
Created by: Haiyang (Skyler) Chang (｡•̀ᴗ-)✧
Method: COORDINATES
Date: 04/01/18
EOF
####checking the input file##########

echo -e "Please input sequences file in Fasta Format"
if [ $# -ne 1 ]
then
echo -e "Error input: Just one file with many fasta format sequences!"
exit 1
fi
###################  find the best 2 anchor sets ##################
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
id=($(Rscript check_ratio1.R))
rm -rf result
}
################Find seq's coors############
function Findseqcoors(){
id1=${id[0]}
id2=${id[1]}
awk 'FNR!=NR && FNR==1 {next} 1' anchor_group{$id1,$id2} > an_set.txt;
awk 'NR==FNR{a[$0]=1;next}!a[$0]' an_set.txt CLU.txt > CLU_wo_an_set.txt;
for i in {anchor_group$id1,anchor_group$id2}
do perl seq_coordinate.pl $i CLU_wo_an_set.txt > Coor$i.txt
done
paste Cooranchor_group* > coorpairs.txt
awk '{print $1"-"$3,$2"-"$4 }' coorpairs.txt > coors.txt
perl seq_cluster.pl coors.txt CLU_wo_an_set.txt > Anccp_clu.txt
sed 's/^ //g' Anccp_clu.txt > Anccp_clus.txt;
rm CLUanchor_group*.txt
rm anchor_group*
rm Cooranchor_group*.txt
rm coorpairs.txt
rm coors.txt
rm Anccp_clu.txt
rm CLU_wo_an_set.txt
}

############ cut Clusters by their coors########
############ place in dir include Clusters.txt##########
function Cutclusters(){
awk '/^>/ { gsub($1,">'$clustername'."$1)} {print $0}' Anccp_clus.txt > Anc_label_clus$clustername.txt
grep '^>' Anc_label_clus$clustername.txt | sed "s/$/ ancdis$ancdis/g" | sed 's/>//g' > Anchors$clustername.txt
awk '$1 ~ /^>/ { close(file); file = "cluster"$1".txt"} {print > file;}' Anc_label_clus$clustername.txt
}



#awk '/^>/{$0=$0",$anchorpair"}1' Anccp_clus.txt > Anc_label_clus.txt
#grep '^>' Anc_label_clus.txt | sed "s/$/,ancdis$ancdis/g" > Anchorsmothers.txt
#awk '$1 ~ /^>/ { file = "cluster"$1".txt"} {print > file}' Anc_label_clus.txt



########place in dir include all cluster#########
function Copyscripts(){
ancdis=`expr $ancdis - 1`
mkdir -p ../anchordis${ancdis}/Conwayreportanchordis${ancdis}
cp ~/Desktop/shell_script/a.out ~/Desktop/shell_script/SelectDNAII.cpp ~/Desktop/shell_script/stat.h ~/Desktop/shell_script/stat.cpp ../anchordis${ancdis}
cp ~/Desktop/shell_script/seq_coordinate.pl ~/Desktop/shell_script/seq_cluster.pl ../anchordis${ancdis}
cp ~/Desktop/shell_script/index_anchor.pl ~/Desktop/shell_script/check_ratio1.R  ../anchordis${ancdis}
}


######Find anchors###########
echo -e "*******Prepare the Conway Algorithm"
mkdir -p Output_coor
g++ -lm -O3 SelectDNAII.cpp stat.cpp
seqcount=$(grep ">" $1 | wc -l)
seqlength=$(awk 'NR==2 {print length($0)}' $1)
ancdis=`expr $seqlength - 1`
echo "Total $seqcount sequences and length is $seqlength "
echo -e "********Running Conway Crossover Operator********"
mkdir -p Output_coor/anchordis${ancdis}/Conwayreportanchordis${ancdis}
cp ./a.out SelectDNAII.cpp stat.h stat.cpp $1 Output_coor/anchordis${ancdis}
cd Output_coor/anchordis${ancdis}
./a.out $1 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
bdac=$(grep '1' best.pcld | wc -l)
while [ "$bdac" -gt "28" ]
do
echo $bdac
ancdis=`expr $ancdis - 1`
echo $ancdis
pwd
cd ../
mkdir -p anchordis${ancdis}/Conwayreportanchordis${ancdis}
pwd
cp ../a.out ../SelectDNAII.cpp ../stat.h ../stat.cpp ../$1 anchordis${ancdis}
cp ../index_anchor.pl ../check_ratio1.R anchordis${ancdis}
cp ../seq_coordinate.pl ../seq_cluster.pl anchordis${ancdis}
cd anchordis${ancdis}
grep -v '^>' $1 > CLU.txt
./a.out $1 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
bdac=$(grep '1' best.pcld | wc -l)
done
echo "Start anchor distance is $ancdis"
Findanchorsets
Findseqcoors
grep '^>' Anccp_clus.txt | sed 's/>//g'| sed "s/$/ ancdis$ancdis/g" > Anchorsmothers.txt
awk '$1 ~ /^>/ { file = "cluster"$1".txt"} {print > file}' Anccp_clus.txt
for i in cluster*
do
k=$(echo $i | sed 's/>//g')
sed 1d $i | sed 's/^ //g' >> $k
rm $i
done
pwd

########################Clustering the original input large-scale sequences#################

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
echo "ancdis is $ancdis"
./a.out run.fasta 100 $ancdis 20 > Conwayreportanchordis${ancdis}/Conway_report.txt
rm run*
bdac=$(grep '1' best.pcld | wc -l)

if [ "$bdac" -lt "29" ]
then
echo "bdac is $bdac"
cp $j CLU.txt
Findanchorsets
Findseqcoors
echo "$j"
clustername=$(echo "${j%.*}" | sed 's/[a-z]*//g')
rm $j
echo "name is $clustername"
Cutclusters
else
mv $j XXX$j
fi
done



for i in cluster*
do
k=$(echo $i | sed 's/>//g')
sed 1d $i | sed 's/^ //g' >> $k
rm $i
done

files=(XXX*.txt)
if [ -e "${files[0]}" ]
then
for i in XXX*.txt
do
k=$(echo $i | sed 's/XXX//g')
#echo "k"
mv $i $k
done
fi

n=$(find . -name "cluster*"  -type f -exec bash -c '[[ $(wc -l < "$1") -gt 150 ]] && echo "$1"' _ '{}' \; | wc -l)
echo "n is $n"
if [ $n -eq 0 ]
then
echo "*******yes in the end ******"
break
fi
done
#############Collect the anchor information######################
pwd
cd ../
find . -type f -name 'Anchors*.txt' -exec cat {} + >> outputanchors.txt
awk -F' ' '{print $3,$1,$2}' outputanchors.txt | sort -gk1,1r -gk2,1 > anchorsfile.txt
find . -type f -not -name 'cluster*.txt' -not -name 'anchor*file.txt' -print0 | xargs -0 rm --














