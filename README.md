# Metascan


> _**Databases and auxillary files can be found -> [here](https://www.microbiology.science.ru.nl/gcremers/)_ <- or ->[here](https://zenodo.org/record/6365663)<-**


> _**The file used for coverage should be parsed to the following format (Checkm Optional):_
> 
> `binname (tab) coverage (tab) checkm_lineage`
> 
> _Where binname is the name of the bin in the directory._**






**Metabolic scanning and annotation of Metagenomes**

Metascan is an Metagenomic Scanning and Annotation tool, with an emphasis on metabolic genes.
The heart of Metascan is a set of metabolic core genes, that are used to paint a picture of the metabolic capacity of the sample.
Furthermore, it utilizes the Kegg pathways for a complete metabolic overview of each sample.

Samples can be analyzed as eiter binned or unbinned metagenome.

Metascan consists of a perl script, a few auxillary (text-)files and a set of HMM profiles, created by clustering TrEmbl proteins, based on Kegg K-numbers.

**INSTALLATION pointers:**

Since it is a Prokka adaptation, it should be able run on any system that can already run Prokka, just by downloading the script and the databases. The only modification that needs to be done is to direct the script to the right location of the database.

The conda Prokka environment can also be used.

However, when using the Prokka conda environment, there can be some issues. 

SignalP has to be requested before it can be downloaded from the website and can be found here: https://services.healthtech.dtu.dk/service.php?SignalP-5.0. It contains a bin file that can be put in PATH

The maximum version for hmmpress is 3.1b2 and can be found  here http://eddylab.org/software/hmmer3/3.1b2/
Higer versions are too stricked to index the Metascan databases.
Version hmmer-3.1b2-linux-intel-x86_64.tar.gz contains precompiled bin files.

Once the databases are indexed, the version doesn really matter anymore. That is why there is no maxiumum version set to this in the script.

When the cmpress and BLASTP don't have the right path, the easiest way to fix this is to run: prokka --listdb
This will show a  line like this:
[09:44:41] Looking for databases in: /usr/local/bioinfo/prokka/db

in lines 51, you should enter that location in the $prokkaloc placeholder.

Krona is not needed to run Metascan right now, but it will leave you with some additional files if you don't.

**So, what can it do?**


- First and foremost, Metascan is intended to get an overview of the metabolic process within a given sample (metagenome). This could be an unbinned core assembly, but it can also be a binned metagenome. In which case you should not forget to include the unbinned leftover contigs as a bin in your analysis as these are still part of the metagenome.

- Besides analysing the sample as a whole, Metascan can also be used to completely annotate genomes (with an emphasis on metabolic processes).

- Thirdly, Metascan can be used to retrieve specific genes from a metagenome. Users can also submit their own (not necessarily metabolic genes) for Metascan to search and retrieve in fasta format.

- Fourth (coming up in metascan2) is the option to search for viral contigs/areas within the metagenome (although not yet fully validated).

**How does Metascan work?**

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

**Can I add my own genes?**

Yes, you can. Metascan has an option `--hmms` that will able you to use your own HMM profile. You can simultaneously choose to use it with all the datasets or none at all. The last option makes Metascan convenient to look for specific genes of interest in (large meta-) genomes. You will however need a HMM profile to do so. If you already have one, or you found one online, then Bob’s your uncle.
If however, you have a new gene, or you want to update an old one, you’ll need to make one yourself. For this, the first thing you need to do is to gather all the relevant protein fastas of your gene of interest. Once you have those, you need to align them. The resulting alignment can then be used with hmmbuild to create a hmm profile. After indexing the HMM profile with hmmpress your profile is ready to go.
Please be aware that multiple HMMs need to be in one file. You can use the command `cat` for this if you have multiple profiles.

When you are finished, add 1 line to the start of the HMM profile:

`#CYCLE (tab) name`

In order for Metascan to be able to output the generated data into the overview files

**Dependencies:**

Prokka:
- parallel    MINVER  => "20130422"
- prodigal    MINVER  => "2.6",     MAXVER  => "2.69"
- hmmsearch   MINVER  => "3.1"
- hmmpress    MINVER  => "3.1", #3.2 causes problems in building databases because of similar DESC fields
- tbl2asn     MINVER  => "24.3"
- rnammer     MINVER  => "1.2"
- barrnap     MINVER  => "0.4"
- blastn      MINVER  => "2.1", #This is actually 2.2, but since were at 2.10 now, the comp thinks were at 2.1 again
- blastp      MINVER  => "2.1"
- signalp     MINVER  => "3.0" max 5.0
- aragorn     MINVER  => "1.2"
- minced      MINVER  => "1.6"
- cmscan      MINVER  => "1.1"
- cmpress     MINVER  => "1.1"


Optional:
- bwa         MINVER  => "0.7"
- checkm      MINVER  => "1"
- samtools    MINVER  => "1.1" #This is actually 1.6, but since were at 1.13 now


To be removed from Metascan:
- makeblastdb MINVER  => "2.2"
- metabat2    MINVER  => "2.12"

