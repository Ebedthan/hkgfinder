# butler <img src="img/logo.png" align="right" width="120"/>

[![License: GPL v3](https://img.shields.io/badge/License-MIT-blue.svg)](https://www.gnu.org/licenses/MIT)


## Find housekeeping genes in prokaryotic (meta)genomes


### Introduction



### Installation

You will have first to install [Prodigal](https://github.com/hyattpd/Prodigal) and [HMMER 3](https://hmmer.org) to be able to run butler.


#### Install from source

```bash
# Download butler development version
git clone https://github.com/Ebedthan/butler.git butler

# Navigate to directory
cd butler

# Install with poetry: see https://python-poetry.org
poetry install --no-dev

# Enter the Python virtual environment with
poetry shell

# Test butler is correctly installed
butler -h
```

If you do not want to go into the virtual environment just do:

```bash
poetry run butler -h
```

## Test

* Type `butler -h` and it should output something like:

```
usage: butler [options] [<FILE>]

optional arguments:
  -o [FILE]      output result to FILE [stdout]
  -a             activate anonymous mode [false]
  -x             deactivate protein-coding gene prediction [false]
  --faa FILE     output matched proteins sequences to FILE
  --fna FILE     output matched DNA sequences to FILE
  -t INT         number of threads [1]
  -q             decrease program verbosity
  -v, --version  show program's version number and exit
  -h, --help     show this help message and exit
```


## Invoking butler

```
butler --faa housekeeping.faa --fna housekeeping.fna file.fa.gz
```

  
## Bugs

Submit problems or requests to the [Issue Tracker](https://github.com/Ebedthan/butler/issues).


## Dependencies

### Mandatory

* [**Prodigal**](https://github.com/sib-swiss/pftools3)  
  Used for protein-coding gene prediction.    
  *Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119*

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*


## Licence

[MIT](https://github.com/Ebedthan/butler/blob/main/LICENSE).


## Author

* [Anicet Ebou](https://orcid.org/0000-0003-4005-177X)

