"""Loads a json file containing settings that the ISAMBARD use.

A settings.json file must be made in this directory, template_settings.py can be copied and paths for your system
can be added. Json files look like python dictionaries, and in fact all json files are valid python dictionary code
(although this isn't necessarily true in reverse). The settings.json file is ignored by git to it avoid being
overwritten with another systems paths. When the isambard is imported a isambard.settings dictionary is generated, which
contains the user chosen settings, but it can be explicitly loaded by an internal module to get these values too.
"""

import json
import os
import subprocess

global_settings = None
package_dir = os.path.dirname(os.path.abspath(__file__))


def configure():
    subprocess.call(['python',
                     os.path.join(package_dir, 'configure.py'),
                     '-o'])
    load_global_settings()
    return


def load_global_settings():
    with open(os.path.join(package_dir, 'settings.json'), 'r') as settings_f:
        global global_settings
        settings_json = json.loads(settings_f.read())
        if global_settings is None:
            global_settings = settings_json
            global_settings[u'package_path'] = package_dir
        else:
            for k, v in settings_json.items():
                if type(v) == dict:
                    global_settings[k].update(v)
                else:
                    global_settings[k] = v
    global_settings['scwrl']['available'] = None
    global_settings['dssp']['available'] = None

if 'settings.json' not in os.listdir(package_dir):
    print('No configuration file (settings.json) found in {}.\nRunning configure.py...\n'.format(package_dir))
    configure()
else:
    load_global_settings()


__author__ = 'Christopher W. Wood'
