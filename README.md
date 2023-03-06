# Models_Epidemiology

[Epidemiological Models and visualization using shiny](https://marcosch.shinyapps.io/epidemiologicmodels_shiny/?fbclid=IwAR3TaFEqaH8Q0ouyNG_2STH_Mkq1fJwl1tXM4P72xQg-NCRlSaZku007wuk)

fastafile='Aedes-aegypti-LVP_AGWG_PEPTIDES_AaegL5.2.fa'
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $fastafile > linear.fasta
#Do grep
#Listgenes.txt contains the list of the union (9231 genes)  
IDs='Listgenes.txt'
while read IDS ; do grep -m1 "\b$IDS\b" linear.fasta ; done < $IDs  |  sed 's/\t/\n/g' > AedesFastaFilt.fa
