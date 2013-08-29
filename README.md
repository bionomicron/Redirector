## About
The Redirector framework is a set of methods based on a steady state model of bacterial metabolism (Flux Balance Analysis FBA) for finding metabolic engineering target to optimize microbial factories, or a desired flux phenotype.  The core concept of the Redirector approach is to model the impact of alterations to the organisms genetic / enzyme level components on the flux of metabolites through its metabolic network by reconstructing the metabolic objective.

## Important
This software package has only just been opened to general usage so please be aware that it has not been thoroughly tested by public users.  I will be working to make it less buggy, better commented and easier to use.  Constructive feedback and suggestions are appreciated but development time is limited, thanks for your patience. 

## Guides
*See the Redirctor GIT Wiki
'http://github.com/bionomicron/Redirector/wiki'
*Also try auto documentation from in line python documentation

## Useful Links
* PuLP-or: http://code.google.com/p/pulp-or/
* Installing PuLP-or: http://code.google.com/p/pulp-or/wiki/InstallingPulpatHome

### Key Redirector objects
* ModelFactory: Factory object for parsing and generating LinearModel objects from flat files
* LinearModel: Main data file for flux balance analysis and objective construction models
* OptimizationControlRedirector: Analysis object which performs progressive target discovery

## Usage

### Redirector Config File
* Default config file is 'Redirector.config'
* Format: standard python config setup
* Changing '''Redirector Model''' configuration in this file will change the setting for all other configurations
* Setting up different Config Tags will allow for easy development and tracking of all framework features related to a particular optimization.

### Command line Redirector
* Command line help information
```bash
  $ python Redirector.py -h
```

* Most used command line inputs information
``` bash
  $ python Redirector.py -c <config file> -n <configuration tag>\
			---model_config <model config file> -m <model name>\
			-b <biological objective> -s <production objective>\
			 --iter <number of iterations> --sn <search neighborhood>\
			--section <subsection of loaded model to use>\
			--control <flat,binary,sense>\
			--simocontrol <number of simulatanious metabolic alterations to a target>\
			--report
```

* Pre-made configuration 
``` bash
  $ python Redirector.py --n <configuration tag>
```

* Example configuration using very simple model (included in repository)
'''
$ python Redirector.py -n "Test Simple Model"
'''

* Example command line of iAF1260 model test configuration
``` bash
$ python Redirector.py -n "Test iAF1260"
```
