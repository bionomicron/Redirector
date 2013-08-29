
## Redirector Usage

### Redirector Config File
* Default config file is Redirector.config
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
  $ python Redirector.py -m <model name> -b <biological objective>\
			-s <production objective> -c <model config file>\ 
			-n 'Default' --iter <number of iterations> --sn <search neighborhood>\
			--section <subsection of loaded model to use>\
			--control <flat,binary,sense>\
			--simocontrol <number of simulatanious metabolic alterations to a target>\
			--report
```

* Pre-made configuration 
``` bash
  $ python Redirector.py --n <configuration tag>
```

* Example command line of iAF1260 model test configuration
``` bash
$ python Redirector.py -n "Test iAF1260"
```
