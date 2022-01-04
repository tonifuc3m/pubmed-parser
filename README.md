# biocreative-background-set-pipeline <img src="images/biocreative-logo.png" alt="Logo" align="right"  width="150" height="150">

<p align="left">
    Library to download PubMed abstracts with metadata. Originally created to obtain the DrugProt (BioCreative VII) background set.
    <br />
    <a href="https://github.com/tonifuc3m/biocreative-background-set-pipeline"><strong>Explore the docs »</strong></a>
</p>

## Requirements

+ Python3
+ [biopython](https://biopython.org/wiki/Download)

## Usage

It has 2 modes:

 + *get_pmids* mode. This mode is intented to be used when we have a set of PubMed queries and we want to extract the PMIDs that match them. It returns a list (or several lists) of PMIDs.
 + *fetch* mode. This mode receives a list of PMIDs and downloads the PubMed titles and abstracts together with their metadata. It returns a tab-separated file (or several ones) with PMID, title, abstract, PMC id, MeSH terms and language. It also stores the complete object downloaded from PubMed into a JSON file (or several files). The titles and abstracts have UFT-8 encoding, with NFKC Unicode normalization and all whitespaces are normalized (meaning, there are no tabs, new lines, etc).

To modify the execution mode, line 241 must be changed.


#### Usage in the *get_pmids* mode

1. Go to line 241 and write 
```
    mode = 'get_pmids'
```

2. Execute python code
```
python get_background.py -i /path/ -o toy-data/queries-example --logfile ~/outfolder/log.log
```

Script Arguments
+ ```-i```: directory where the pubmed queries file is. It is also the directory where the output will be created
+ ```-o```: name of file with pubmed queries
+ ```--logfile```: path to logfile


#### Usage in the *fetch* mode

1. Go to line 241 and write 
```
    mode = 'fetch'
```

2. Execute python code
```
python get_background.py --input toy-data/pmids-example --output ~/outfolder --logfile ~/outfolder/log.log
```

Script Arguments
+ ```--input```: text file with the list of PMIDs you want to download, one per line
+ ```--output```: folder where we will store the output
+ ```--logfile```: path to logfile 


<p align="center">
    <a href="https://github.com/tonifuc3m/biocreative-background-set-pipeline/issues">Report Bug</a>
    ·
    <a href="https://github.com/tonifuc3m/biocreative-background-set-pipeline/issues">Request Feature</a>
</p>
