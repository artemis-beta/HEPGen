from hepgen.Configurables import Generator

Generator().OutFile = 'KPMMMC.dtf'
Generator().OutTree = 'KPiMuMuTree'
Generator().DecID = 'SK020001'
Generator().nEvts = 2000
Generator().Energy = 100
