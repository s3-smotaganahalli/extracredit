# extracredit
extra credit repository for bio312
### Lab 3 BLAST SEARCH
```
git clone https://github.com/Bio312/lab03-$MYGIT: this clones the lab 3 repository

create the BLAST database
```
cd ~/lab03-$MYGIT: this takes you to the lab 3 folder
gunzip proteomes/*.gz:this uncompresses the proteome
cat  proteomes/*.faa > allprotein.fas: this command puts all the protein sequences into a single file
makeblastdb -in allprotein.fas -dbtype prot: makes BLAST database

mkdir ~/lab03-$MYGIT/globins
cd ~/lab03-$MYGIT/globins
pwd: navigate to the new directory produced

Download the query protein and produce a blast search
```
ncbi-acc-download -F fasta -m protein "NP_000549.1"
blastp -db ../allprotein.fas -query NP_000549.1.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typical.out
grep -c H.sapiens globins.blastp.detail.out:can count total hits in the file
awk '{if ($6< 1e-30)print $1 }' globins.blastp.detail.out > globins.blastp.detail.filtered.out: filter out outputs
```
