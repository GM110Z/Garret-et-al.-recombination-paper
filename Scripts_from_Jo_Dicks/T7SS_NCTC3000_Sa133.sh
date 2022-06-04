#!/bin/bash
echo Identifying putative T7SS regions within 133 NCTC Staphylococcus aureus genome assemblies

strain=( 1803 2669 3750 3761 4136 4137 4163 5655 5656 5657 5658 5660 5663 5664 6131 6133 6134 6135 6136 6137 6507 6571 6966 7121 7361 7414 7415 7428 7445 7446 7447 7485 7712 7718 7791 7856 7878 7972 7988 8004 8178 8317 8325 8399 8507 8530 8531 8532 8723 8724 8725 8726 8765 8766 9309 9369 9393 9546 9547 9551 9552 9553 9555 9556 9611 9612 9613 9614 9752 9944 10035 10344 10345 10399 10442 10443 10649 10652 10654 10655 10656 10657 10702 10703 10724 10788 10804 10833 10988 11150 11561 11939 11940 11962 11963 11965 12035 12232 12233 12880 12981 13131 13132 13133 13134 13135 13136 13137 13138 13139 13140 13141 13142 13143 13196 13277 13297 13298 13299 13300 13373 13394 13395 13434 13435 13552 13616 13626 13758 13811 13812 13813 13841 )

for index in $(seq 0 132);

do 


#locate putative essC, esaD, esaG and focA sequences in the current assembly

nhmmer --tblout NCTC${strain[index]}_essC.tblout essC.hmm NCTC${strain[index]}.fasta
nhmmer --tblout NCTC${strain[index]}_esaD.tblout esaD.hmm NCTC${strain[index]}.fasta
nhmmer --tblout NCTC${strain[index]}_esaG.tblout esaG.hmm NCTC${strain[index]}.fasta
nhmmer --tblout NCTC${strain[index]}_focA.tblout focA.hmm NCTC${strain[index]}.fasta

echo "Identification of putative T7SS regions in 133 NCTC Staphylococcus aureus strains" >> t7ssout.txt

essChere=0
findessC=$(grep -v "^#" NCTC${strain[index]}_essC.tblout | head -n 1 | awk -v awkvar=$essChere '{if ($13 < 1.0e-100){awkvar++;print awkvar}}')
if [ "$findessC" == "" ] 
then
   echo essC not found on strain NCTC ${strain[index]} >> t7ssout.txt
elif [ $findessC -eq 1 ]
then 
   contigessC=$(grep -v "^#" NCTC${strain[index]}_essC.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $1}}')
   startessC=$(grep -v "^#" NCTC${strain[index]}_essC.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $9}}')
   endessC=$(grep -v "^#" NCTC${strain[index]}_essC.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $10}}')
   echo essC gene found on strain NCTC ${strain[index]}: contig $contigessC between co-ordinates $startessC and $endessC >> t7ssout.txt
fi

esaDhere=0
findesaD=$(grep -v "^#" NCTC${strain[index]}_esaD.tblout | head -n 1 | awk -v awkvar=$esaDhere '{if ($13 < 1.0e-100){awkvar++;print awkvar}}')
if [ "$findesaD" == "" ] 
then
   echo esaD not found on strain NCTC ${strain[index]} >> t7ssout.txt
elif [ $findesaD -eq 1 ]
then 
   contigesaD=$(grep -v "^#" NCTC${strain[index]}_esaD.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $1}}')
   startesaD=$(grep -v "^#" NCTC${strain[index]}_esaD.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $9}}')
   endesaD=$(grep -v "^#" NCTC${strain[index]}_esaD.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $10}}')
   echo esaD gene found on strain NCTC ${strain[index]}: contig $contigesaD between co-ordinates $startesaD and $endesaD >> t7ssout.txt
fi

focAhere=0
findfocA=$(grep -v "^#" NCTC${strain[index]}_focA.tblout | head -n 1 | awk -v awkvar=$focAhere '{if ($13 < 1.0e-100){awkvar++;print awkvar}}')
if [ "$findfocA" == "" ] 
then
   echo focA not found on strain NCTC ${strain[index]} >> t7ssout.txt
elif [ $findfocA -eq 1 ]
then 
   contigfocA=$(grep -v "^#" NCTC${strain[index]}_focA.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $1}}')
   startfocA=$(grep -v "^#" NCTC${strain[index]}_focA.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $9}}')
   endfocA=$(grep -v "^#" NCTC${strain[index]}_focA.tblout | head -n 1 | awk '{if ($13 < 1.0e-100){print $10}}')
   echo focA gene found on strain NCTC ${strain[index]}: contig $contigfocA between co-ordinates $startfocA and $endfocA >> t7ssout.txt
fi

esaGcount=0
coords=( )
grep -v "^#" NCTC${strain[index]}_esaG.tblout | awk '{if ($13 < 1.001){print $0}}' > temp_esaG.tblout
while read line
do
   contigesaG=$(echo "$line" | awk '{print $1}')
   startesaG=$(echo "$line" | awk '{print $9}')
   endesaG=$(echo "$line" | awk '{print $10}')
   scoreesaG=$(echo "$line" | awk '{print $13}')
   if [ "$contigesaG" == "$contigessC" ] 
   then
      esaGcount=$((esaGcount+1))
      coords+=($startesaG)
      coords+=($endesaG)
   fi
done < temp_esaG.tblout
echo "$esaGcount esaG-like sequence(s) found for strain NCTC ${strain[index]} on the recorded essC contig" >> t7ssout.txt
echo "" >> t7ssout.txt

done
