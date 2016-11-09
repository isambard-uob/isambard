import inspect as _inspect
import os as _os
import subprocess as _subprocess
import sys as _sys
import pyximport; pyximport.install()

_sys.path.append('')
_starting_dir = _os.getcwd()
_cmd_folder = _os.path.realpath(_os.path.abspath(_os.path.split(_inspect.getfile(_inspect.currentframe()))[0]))
_os.chdir(_cmd_folder)

try:
    if 'settings.json' not in _os.listdir(_cmd_folder):
        print('No configuration file (settings.json) found in {}.\nRunning configure.py...\n'.format(_cmd_folder))
        _subprocess.call(['python', 'configure.py'])
    from settings import global_settings
    import ampal
    import ampal.specifications as specifications
    from ampal import analyse_protein
    from ampal import interactions
    import buff
    import external_programs
    import tools
    import tools.geometry as geometry
    import optimisation
    with open('logo.txt', 'r') as inf:
        logo = ''.join(inf.readlines()[:51])
finally:
    _os.chdir(_starting_dir)

__version__ = "2016.2"
