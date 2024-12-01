# extra credit
extra credit repository for BIO312
### Lab 3 Blast Search
The goal of this lab is to use the Basic Local Alignment Search Tool to find homologous sequences throughout the proteomes of eleven gnathostome species, such as orthologs and paralogs of your own gene family. Using Unix-based tools, we learnt how to create BLAST databases, run BLAST searches, filter results according to e-value criteria, and examine the distribution of paralogs in each species.
```
git clone https://github.com/Bio312/lab03-$MYGIT: this clones the lab 3 repository

create the BLAST database

cd ~/lab03-$MYGIT
#this takes you to the lab 3 folder
gunzip proteomes/*.gz
#This uncompresses the proteome
cat  proteomes/*.faa > allprotein.fas
#this command puts all the protein sequences into a single file
makeblastdb -in allprotein.fas -dbtype prot
#makes BLAST database

mkdir ~/lab03-$MYGIT/peptidase
cd ~/lab03-$MYGIT/peptidase
pwd
#this allows us to navigate to the new directory produced

Download the query protein and produce a blast search

ncbi-acc-download -F fasta -m protein "NP_000549.1"
blastp -db ../allprotein.fas -query NP_000549.1.fa -outfmt 0 -max_hsps 1 -out peptidase.blastp.typical.out
grep -c H.sapiens peptidase.blastp.detail.out
#can count total hits in the file
awk '{if ($6< 1e-30)print $1 }' peptidase.blastp.detail.out > peptidase.blastp.detail.filtered.out
#filter out outputs
```
### Lab 4 Gene family sequence alignment
Working with FASTA files and aligning sequences with programs like MUSCLE, T-Coffee, and AlignBuddy are the main topics of this lab. In order to compute metrics such as alignment width, invariant locations, and percent identity, we interpreted the generated alignments. 
```
cd lab04-$MYGIT
pwd
mkdir /mygene
mkdir /home/bio312-user/lab04-s3-smotaganahalli/mygene
#now the directory is set up 

# the followinf code obtains the sequences that are in the BLAST file
seqkit grep --pattern-file /home/bio312-user/lab03-s3-smotaganahalli/mygene/peptidase.blastp.detail.filtered.out /home/bio312-user/lab03-s3-smotaganahalli/allprotein.fas | seqkit grep -v -p "carpio" > /home/bio312-user/lab04-s3-smotaganahalli/mygene/peptidase.homologs.fas

# this code perform global multiple alignment
muscle -align ~/lab04-$MYGIT/peptidase/peptidase.homologs.fas -output ~/lab04-$MYGIT/peptidase/peptidase.homologs.al.fas
ls /home/bio312-user/lab04-s3-smotaganahalli/mygene/

# use this to save and view
alv -kli  ~/lab04-$MYGIT/mygene/peptidase.homologs.al.fas | less -RS
alv -kli --majority ~/lab03-$MYGIT/peptidase/peptidase.homologs.al.fas | less -RS

# if you want a different view like print to a large PDF
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/peptidase.homologs.al.fas

# this will calculate length of alignment
alignbuddy  -al  ~/lab04-$MYGIT/peptidase.homologs.al.fas

# this will calculate the length after removing any gaps
alignbuddy -trm all  ~/lab04-$MYGIT/peptidase.homologs.al.fas | alignbuddy  -al

# this is the length after removing invariants
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/peptidase.homologs.al.fas | alignbuddy  -al

# average percent identity using tcoffee
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/peptidase.homologs.al.fas -output sim

# average percent identity using alignbuddy
alignbuddy -pi ~/lab04-$MYGIT/peptidase.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
```
### Lab 5 Gene Family Phylogeny using IQ-TREE
In this lab, we created a phylogenetic tree for our gene family and used IQ-TREE to handle the alignment data. This program determines the optimal substitution model, estimates branch lengths, and computes bootstrap support. The generated tree was rooted using midpoint or outgroup techniques for meaningful interpretation, and it was displayed as a cladogram or phylogram.
```
git clone https://github.com/Bio312/lab05-$MYGIT
mkdir ~/lab05-$MYGIT/peptidase
cd ~/lab05-$MYGIT/peptidase
#the directory is made and ready to use

#removes any sequences containing duplicate label tags
sed 's/ /_/g'  ~/lab04-$MYGIT/peptidase/peptidase.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas

#estimates ultrafast bootstrap support  levels
iqtree -s ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas -bb 1000 -nt 2

#updates newick display
nw_display ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile

#draws the tree unrooted
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile.pdf 0.4 15

#reroots the tree
gotree reroot midpoint -i ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile

#look at the rooted tree and make it a graphic
nw_order -c n ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile  | nw_display -
nw_order -c n ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile.svg -

#coverts svg image to a pdf
convert  ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile.svg

#Switch the view to a cladogram
nw_nw_order -c n ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > 
convert ~/lab05-$MYGIT/globins/globins.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/globins/globins.homologsf.al.midCl.treefile.pdforder -c n ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.midCl.treefile.svg -
nw_display ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.midCl.treefile.svg

#Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -25608.830
Optimal log-likelihood: -25608.823
Proportion of invariable sites: 0.050
Site proportion and rates:  (0.197,0.160) (0.310,0.543) (0.352,1.249) (0.091,3.965)
Parameters optimization took 1 rounds (0.519 sec)
BEST SCORE FOUND : -25608.823
Creating bootstrap support values...
Split supports printed to NEXUS file /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.splits.nex
Total tree length: 18.897

Total number of iterations: 106
CPU time used for tree search: 310.401 sec (0h:5m:10s)
Wall-clock time used for tree search: 157.547 sec (0h:2m:37s)
Total CPU time used: 2021.872 sec (0h:33m:41s)
Total wall-clock time used: 1023.575 sec (0h:17m:3s)

#Computing bootstrap consensus tree...
Reading input file /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.splits.nex...
37 taxa and 159 splits.
Consensus tree written to /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.contree
Reading input trees file /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.contree
Log-likelihood of consensus tree: -25611.075

Analysis results written to: 
  IQ-TREE report:                /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.iqtree
  Maximum-likelihood tree:       /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.treefile
  Likelihood distances:          /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.splits.nex
  Consensus tree:                /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.contree
  Screen log file:               /home/bio312-user/lab05-s3-smotaganahalli/peptidase/peptidase.homologsf.al.fas.log

Date and Time: Thu Sep 26 19:32:26 2024
```
### Lab 6 Reconciling a Gene and Species Tree
To distinguish between orthologs and paralogs, we used Notung to estimate duplication and loss events in order to reconcile gene and species trees in this lab. In order to run reconciliations and analyze results like duplication and loss counts using graphical tools, we learnt how to set up a Python environment as well. 
```
git clone https://github.com/Bio312/lab06-$MYGIT
cd lab06-$MYGIT
#set up

# this is the code for the Notung package
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar --help
 
#prune command to isolate ingroups
gotree prune -i ~/lab06-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile -o ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta

#name and save tree
mkdir ~/lab06-$MYGIT/peptidase
cp ~/lab05-$MYGIT/peptidase/peptidase.homologs.al.mid.treefile

#Notung Reconciliation
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGITpeptidase/peptidase.homologsf.pruned.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/peptidase/

#internal node search
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -

#generates RecPhyloXML object 
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.ntg --include.species

#view the gene tree reconciliation and convert to PDF
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.svg
convert  -density 150 ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.svg ~/lab06-$MYGIT/peptidase/peptidase.homologsf.pruned.treefile.rec.pdf
```
### Lab 8 Protein Domain Prediction
In order to find Pfam domains in protein sequences, this lab focuses on protein domain prediction using RPS-BLAST. To complete primary analyses for our research projects, we examined functional patterns across lineages, produced domain annotations plotted alongside phylogenies, and examined peptidase proteins.

```
git clone https://github.com/Bio312/lab08-$MYGIT
mkdir ~/lab08-$MYGIT/peptidase && cd ~/lab08-$MYGIT/peptidase
#set up the right directories

#makes a copy of the unaligned sequence and removes stop codons
sed 's/*//' ~/lab04-$MYGIT/mygene/peptidase.homologs.al.fas > ~/lab08-$MYGIT/peptidase/peptidase.homologs.al.fas

#RSP Blast command
rpsblast -query ~/lab08-$MYGIT/peptidase/peptidase.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001

#copy final tree and run R-script
cp ~/lab05-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile ~/lab08-$MYGIT/peptidase
Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out ~/lab08-$MYGIT/peptidase/peptidase.tree.rps.pdfRscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/peptidase/peptidase.homologsf.al.fas.treefile ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out ~/lab08-$MYGIT/peptidase/peptidase.tree.rps.pdf

#delinates annotations in blast file
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out | tail -n +2 | less -S

#shortcut so you don't have to count annotations and calculations by hand
cut -f 1 ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out | sort | uniq -c
cut -f 6 ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out | sort | uniq -c
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out |  sort  -k2nr

#pulls out just the evalues and sorts
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out 
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/peptidase/peptidase.rps-blast.out | sort
```
