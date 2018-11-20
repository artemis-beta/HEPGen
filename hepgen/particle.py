import pypdt

_all_particles = list(pypdt.PDT().ids())

def convert_name(name):
    _symbols = {')(' : '_', '(' : '_', ')' : '_',
                '*' : 'star', '+' : 'plus',
                '-' : 'minus'}
    for s in _symbols:
        name = name.replace(s, _symbols[s])
    return name
         
class particle(object):
     def __init__(self, _id):
          _part = pypdt.get(_id)
          self.symbol   = _part.name
          self.name     = convert_name(_part.name)
          self.charge   = _part.charge
          self.ctau     = _part.ctau 
          self.mass     = _part.mass*1000. #MeV
          self.lifetime = _part.get_lifetime()

class ParticleList(object):
    def __init__(self):
        self.__iter_list__ = []
        for p in _all_particles:
            _particle = particle(p)
            setattr(self, _particle.name, particle(p))
            self.__iter_list__.append(_particle.name)
            _particle2 = particle(p)
            _anti = _particle2.name.replace('plus', 'minus') if 'plus' in _particle.name else _particle.name.replace('minus', 'plus')
            _particle2.symbol = _particle2.symbol.replace('+', '-') if '+' in _particle2.symbol else _particle2.symbol.replace('-', '+')
            _particle2.charge = _particle2.charge*-1
            _particle2.name = _anti
            setattr(self, _anti, _particle2)
            self.__iter_list__.append(_anti)

    def __iter__(self):
        for i in self.__iter_list__:
            yield getattr(self, i)

    def get(self, symb):
        for p in self:
            print(p.symbol)
            if p.symbol == symb:
               return p
        raise AttributeError("Could not retrieve particle {}".format(symb))

pdg = ParticleList()
