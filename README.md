# extra credit
extra credit repository for bio312
### Lab 3 Blast Search
```
git clone https://github.com/Bio312/lab03-$MYGIT: this clones the lab 3 repository

create the BLAST database

cd ~/lab03-$MYGIT:this takes you to the lab 3 folder
gunzip proteomes/*.gz: This uncompresses the proteome
cat  proteomes/*.faa > allprotein.fas: this command puts all the protein sequences into a single file
makeblastdb -in allprotein.fas -dbtype prot: makes BLAST database

mkdir ~/lab03-$MYGIT/globins
cd ~/lab03-$MYGIT/globins
pwd: navigate to the new directory produced

Download the query protein and produce a blast search

ncbi-acc-download -F fasta -m protein "NP_000549.1"
blastp -db ../allprotein.fas -query NP_000549.1.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typical.out
grep -c H.sapiens globins.blastp.detail.out:can count total hits in the file
awk '{if ($6< 1e-30)print $1 }' globins.blastp.detail.out > globins.blastp.detail.filtered.out: filter out outputs
```
### Lab 4 Gene family sequence alignment
```
cd lab04-$MYGIT
pwd
mkdir /mygene
mkdir /home/bio312-user/lab04-s3-smotaganahalli/mygene
obtains the sequences that are in the BLAST file
seqkit grep --pattern-file /home/bio312-user/lab03-s3-smotaganahalli/mygene/peptidase.blastp.detail.filtered.out /home/bio312-user/lab03-s3-smotaganahalli/allprotein.fas | seqkit grep -v -p "carpio" > /home/bio312-user/lab04-s3-smotaganahalli/mygene/peptidase.homologs.fas

perform global multiple alignment
muscle -align ~/lab04-$MYGIT/peptidase/peptidase.homologs.fas -output ~/lab04-$MYGIT/peptidase/peptidase.homologs.al.fas
ls /home/bio312-user/lab04-s3-smotaganahalli/mygene/

save and view
alv -kli  ~/lab04-$MYGIT/mygene/peptidase.homologs.al.fas | less -RS
alv -kli --majority ~/lab03-$MYGIT/peptidase/peptidase.homologs.al.fas | less -RS

print to a large PDF
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/peptidase.homologs.al.fas

calculate length of alignment
alignbuddy  -al  ~/lab04-$MYGIT/peptidase.homologs.al.fas

length after removing any gaps
alignbuddy -trm all  ~/lab04-$MYGIT/peptidase.homologs.al.fas | alignbuddy  -al

length after removing invariants
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/peptidase.homologs.al.fas | alignbuddy  -al

average percent identity using tcoffee
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/peptidase.homologs.al.fas -output sim

average percent identity using alignbuddy
alignbuddy -pi ~/lab04-$MYGIT/peptidase.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
```
