# Metascan

Metascan is a Metagenomic Scanning and Annotation tool with an emphasis on metabolic genes.

## **Metabolic scanning and annotation of Metagenomes**

The heart of Metascan is a set of metabolic core genes, that are used to paint a picture of the metabolic capacity of the sample.
Furthermore, it utilizes the Kegg pathways for a complete metabolic overview of each sample.

Samples can be analyzed as either a binned or unbinned metagenome.

Metascan consists of a perl script, a few auxillary (text-)files and a set of HMM profiles, created by clustering TrEmbl proteins, based on Kegg K-numbers.
Recently three scripts have been added. One to download the databases and two to retrieve fastas from the analysed datasets. 

When Metascan was developped, it was based on Prokka, as Prokka already had most of the functionality that was needed and I had never written a script before.
Because of this, Metascan would run if Prokka worked on the sytem.

However, due a cascade of changes starting with a change of tbl2asn to tabel2asn, it is alot more difficult to use Metascan if Prokka can run in an environment or system.
Therefor, Metascan now has its own Conda environment. This update was also used to include a viral sequence algorithm which is VERY! beta (see below).

## **INSTALLATION:**

To install Metascan, first clone Metascan to your system:

```bash
# Clone the git repository to your system;
git clone https://github.com/gcremers/metascan.git
# and cd into this directoy
cd metascan/
# Create an enironment and activate the environment
conda env create -f metascan.yaml -n metascan
conda activate metascan

# NB, the creation throws a number of error messages on my system. They do not seem to harm the installation.

# Unfortunately some hardcoding is required. I hope to find a solution for this at some point.
# Open the metascan script with (for instance) gedit.
gedit ./metascan
# and in line 50, change the path to the directory where you want the databases located in the line
my $databasedir="/path/to/databases/metascan";

# create that directory
mkdir /path/to/databases/metascan

# Put metascan in the bin folder of the environment (change the second part accordingly).
# There is one extra script now for downloading the databases
cp ./metascan /path/to/Anaconda/envs/metascan/bin/
cp ./metascan_db_download /path/to/Anaconda/envs/metascan/bin/
```

Metascan is now installed

Now we have to build the databases

```bash
# Still working from the github directory we created

# Download the databases
metascan . --download
# Index the databases
metascan . --setupdb
# Check databases with 
metascan . --listdb

#if you want to save some space on your system and if you already have a checkm database, you can (re)set the path with
checkm data setRoot <DB_PATH>
see checkm data -h
```

Metascan should now be ready to use.

##### Some notes:


* The script to download databases can be used on its own, to install the databases seperately.

* Databases and auxillary files can also be found -> [here](https://www.microbiology.science.ru.nl/gcremers/)_ <- or ->[here](https://zenodo.org/record/6365663)<-

* The file used for coverage should be parsed to the following format (checkm_lineage optional):

`binname (tab) coverage (tab) checkm_lineage`
 
Where binname is the name of the bin in the directory.

* Although Metascan has to option to use mapping (BWA) and Checkm, I'd generally advice against it, since nowadays, much more specialized tools and pipelines are available that will probably work much faster as well.

* Due to Github issues on my part, the initial Github repository is now named metascan-old.


## **So, what can it do?**


- First and foremost, Metascan is intended to get an overview of the metabolic process within a given sample (metagenome). This could be an unbinned core assembly, but it can also be a binned metagenome. In which case you should not forget to include the unbinned leftover contigs as a bin in your analysis as these are still part of the metagenome.

- Besides analysing the sample as a whole, Metascan can also be used to completely annotate genomes (with an emphasis on metabolic processes).

- Thirdly, Metascan can be used to retrieve specific genes from a metagenome. Users can also submit their own (not necessarily metabolic genes) for Metascan to search and retrieve in fasta format.
There are two scripts added to that will allow you to retrieve fastas from the fasta files, based on either the header or patterns in the bases/amino-acids.

To iterate over multiple directories, go the base directory where your subdirs are and run: `for d in ./*/ ; do (cd "$d" && perl fastaextract_desc_id.pl *.hmm.faa searchterm ../xxx); done`
This will create a file called xxx.fasta, which contains all genes that have 'searchterm' in the header.

If you want proteins containing certain motifs, for instance CxxCH (cytochrome C), you could run (iteratively if desired) `perl fastaextract_desc_id.pl *.hmm.faa c..ch ../cytc`.

- Fourth is the option to search for viral contigs/areas within the metagenome (although not yet fully validated).

## **How does Metascan work?**

Metascan works on a basis of ~180 different metabolic genes, that are key-genes or important genes in metabolic processes. McrA for instance, is needed for the conversion of methyl-CoM into methane. Without mcrA, there is no methanogenesis. Therefore, the amount of mcrA in a sample can be seen as a measure of the potential for Methanogenesis.
The complete key-gene set is divided into 8 main subsets:
- Methane cycle
- Nitrogen cycle
- Sulphur cycle
- Oxygen respiration
- C1 compounds
- Carbon fixation
- Hydrogen
- Miscellaneous metal cycling


Each (most) of these subsets are composed of different processes and these processes are often formed by modules.
Fore more information about the processes and modules I refer to the KEGG website, as the data and setup used in Metascan comes directly from KEGG.
This makes it also easy to load your data into the KEGG website for further analysis and a visual reference.
https://www.genome.jp/kegg/mapper.html

## **Can I add my own genes?**

Yes, you can. Metascan has an option `--hmms` that will able you to use your own HMM profile. You can simultaneously choose to use it with all the datasets or none at all. The last option makes Metascan convenient to look for specific genes of interest in (large meta-) genomes. You will however need a HMM profile to do so. If you already have one, or you found one online, then Bob’s your uncle.
If however, you have a new gene, or you want to update an old one, you’ll need to make one yourself. For this, the first thing you need to do is to gather all the relevant protein fastas of your gene of interest. Once you have those, you need to align them. The resulting alignment can then be used with hmmbuild to create a hmm profile. After indexing the HMM profile with hmmpress your profile is ready to go.
Please be aware that multiple HMMs need to be in one file. You can use the command `cat` for this if you have multiple profiles.

When you are finished, add 1 line to the start of the HMM profile:

`#CYCLE (tab) name`

In order for Metascan to be able to output the generated data into the overview files

## **Viral contigs algorithm:**

_How does it work?_

Viral HMM profiles were downloaded from https://vogdb.org/ and can now be used as a separate subset in the analysis. During the analysis, Metascan marks all the proteins that have a hit with this database. After the analysis is done, Metascan then does a ‘rolling average’ for all the proteins/loci on the contigs over a certain window.



Theoretical DNA stretch, containing possible viral proteins 

|**I**  |A 	|B 	|C 	|D 	|E 	|F 	|G 	|H 	|I 	|J 	|K 	|L 	|M 	|N 	|O 	|P 	|Q 	|R 	|S 	|T 	|U 	|V 	|W     |
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:----:|
|**II** |1 	|0 	|0 	|0 	|0 	|0 	|1 	|0 	|0 	|1 	|1 	|1 	|1 	|0 	|1 	|0 	|1 	|0 	|0 	|0 	|0 	|0 	|0     |
|**III**|0.2 	|0.2 	|0.2 	|0 	|0.2 	|0.2 	|0.2 	|0.4 	|0.6 	|0.6 	|0.8 	|0.8 	|0.8 	|0.6 	|0.6 	|0.4 	|0.4 	|0.2 	|0.2 	|0 	|0 	|0 	|0     |

Line **I**: Loci,
Line **II**: hit yes(1) or no(0),
Line **III**: rolling average over a window of 5. 



This means that for each protein, the average is calculated over the amount of proteins of that window. So with a window of 5, the average is calculated over 2 proteins in front, the protein itself and 2 proteins thereafter. Missing numbers at the ends are counted as 0.
So:
for J, the calculation is (0+0+1+1+1)/5=0.6
for L, the calculation is (1+1+1+1+0)/5=0.8
for A, the calculation is (0+0+1+0+0)/5=0.2

This results in a line of numbers with an average for each locus. From the sample above it is clear that with the area G to Q, there are a lot of viral genes and thus it is an area of interest. Metascan then collect the areas that have a continuous line of positive averages (E to S). When the line reaches 0, the area is ended.

Metascan then outputs all the loci for all the stretches that are found, but only the stretches that reach a certain threshold are thoroughly reported, protein by protein and with the inclusion of fasta files. These stretches and fasta files can then be manually investigated to your own leisure.

Because some viral genes are, as of yet, exclusively found in viruses, while others are also part bacterial genomes, the VOG database comes with an index value for each profile. This number indicates how exclusive the gene is for viruses (3 being only found in viruses, 0 not found in viruses). This number is also taken into account in the calculation used by Metascan. 
