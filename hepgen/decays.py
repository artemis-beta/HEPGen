from hepgen import particle
import os
_mod_loc = os.path.abspath(__file__).replace(__file__, '')

def get_particle(symbol):
    _dict = {'*' : 'star', '+' : 'plus',
             '-' : 'minus'}
    _symb = symbol
    for d in _dict:
        _symb = _symb.replace(d, _dict[d])
    return getattr(particle.pdg, _symb)

class Decay(object):
    def __init__(self, _id):
        pass

    def __str__(self):
         return self.Descriptor

class DecayList:
    def __init__(self):
        self._decays = {}
        self._parse_decay_files()

    def __call__(self, _id):
        print(list(self._decays.keys()))
        for d in self._decays:
            if self._decays[d].ID == _id:
                return self._decays[d]
        raise Exception("Could Not Find Decay")

    def __iter__(self):
        for d in self._decays:
            yield self._decays[d]

    def __str__(self):
        return '\n'.join([self._decays[i].Descriptor for i in self._decays])

    def _parse_decay_files(self):
        from glob import glob
        import os
        _files = glob(os.path.join(os.path.join(_mod_loc, 'decay_files'), '*.dcf'))

        assert len(_files) > 0
     
        for file_addr in _files:
           _tmp = {}
           with open(file_addr, 'r') as f:
               for l in f.readlines():
                   l = l.replace('\n', '')
                   try:
                       _tmp[l.split(':')[0]] = l.split(':')[1]
                   except IndexError:
                       continue

           _decay = Decay(_tmp['ID'])
           _decay.ID = _tmp['ID'].replace(' ','')
           _decay.BR = float(_tmp['Branching Ratio'])
           _decay.Descriptor = _tmp['Decay']
           _decay.Description = _tmp['Description']
           _decay.Mother = get_particle(_decay.Descriptor.split(' -> ')[0].replace(' ', ''))
           _decay.Daughters = [get_particle(i.replace(' ','')) for i in _decay.Descriptor.split(' -> ')[1].split(' ')]

           self._decays[_decay.ID] = _decay

decay_table = DecayList()
