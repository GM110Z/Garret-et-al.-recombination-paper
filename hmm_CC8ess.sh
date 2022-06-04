#!/bin/bash
echo Analysing 31 Staphylococcus aureus T7SS loci using a DNA hmm profile built with hmmbuild

assembly=(NCTC10703	NCTC13133	NCTC13394	NCTC8325
NCTC10442	NCTC10724	NCTC13136	NCTC13395	NCTC8532
NCTC10443	NCTC10833	NCTC13138	NCTC13758	NCTC9309
NCTC10652	NCTC10988	NCTC13139	NCTC13812	NCTC9369
NCTC10654	NCTC12232	NCTC13140	NCTC4136
NCTC10657	NCTC12233	NCTC13141rev	NCTC8004
NCTC10702	NCTC13132	NCTC13196	NCTC8178)

strain=(NCTC10703	NCTC13133	NCTC13394	NCTC8325
NCTC10442	NCTC10724	NCTC13136	NCTC13395	NCTC8532
NCTC10443	NCTC10833	NCTC13138	NCTC13758	NCTC9309
NCTC10652	NCTC10988	NCTC13139	NCTC13812	NCTC9369
NCTC10654	NCTC12232	NCTC13140	NCTC4136
NCTC10657	NCTC12233	NCTC13141rev	NCTC8004
NCTC10702	NCTC13132	NCTC13196	NCTC8178)


#locate esaG gene sequences in the current strain
hmmbuild CC8esaG.hmm esaG_aln.fasta

for index in $(seq 0 30);

do 

#index the current strain for manipulation with HMMER's accessory tool Easel
esl-sfetch --index ${assembly[index]}.gb

#locate esaG gene sequences in the current strain
nhmmer --tblout ${strain[index]}_esaG.tblout CC8esaG.hmm ${assembly[index]}.gb

#extract all esaG-like sequences with an E-value < 1.0e-100
grep -v "^#" ${strain[index]}_esaG.tblout | awk '{if ($13 < 1.0e-100){print $1"/"$9"-"$10, $9, $10, $1}}' | esl-sfetch -Cf ${assembly[index]}.gb - > ${strain[index]}_esaG.fa

#add in turn each esaG-like sequence to a fasta file
cat ${strain[index]}_esaG.fa >> CC8esaG_all.fa

done
