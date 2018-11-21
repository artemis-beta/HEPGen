# HEPGen

A particle physics event simulation library in Python (in Alpha!). This version has compatibility with the Herne configuration system.

## Herne
This version of HEPGen allows event generation and setup to be performed via Herne scripts of options. A Herne App `Generator' has been created which can be configured by creating a script in the following form:

```
from hepgen.Configurables import Generator

Generator().OutFile = 'my_file.dtf'
Generator().OutTree = 'name_of_data_tree'
Generator().DecID   = 'ID of Decay in decay_files'
Generator().nEvts   = number_of_events
Generator().Energy  = boost_energy_mev
```
