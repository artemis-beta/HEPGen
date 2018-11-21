from setuptools import setup

setup(name                =  'HEPGen'                                      ,
      version             =  '0.1.3'                                       ,
      description         =  'Particle Physics Monte Carlo.'               ,
      url                 =  ''                                            ,
      author              =  'Kristian Zarebski'                           ,
      author_email        =  'krizar312@yahoo.co.uk'                       ,
      license             =  'MIT'                                         ,
      packages            =  ['hepgen', 'hepgen.Herne']                    ,
      zip_safe            =  False                                         ,
      include_package_data = True                                          ,
      install_requires    =  ['pypdt', 'matplotlib', 'scipy', 'pyyaml']
     )
