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
home_dir = os.path.expanduser('~')
settings_path = os.path.join(home_dir, '.isambard_settings')


def configure():
    """Runs configure.py in the ISAMBARD directory, creates settings file."""
    subprocess.call(['python',
                     os.path.join(package_dir, 'configure.py'),
                     '-o'])
    load_global_settings()
    return


def load_global_settings():
    """Loads settings file containing paths to dependencies and other optional configuration elements."""
    with open(settings_path, 'r') as settings_f:
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

if not os.path.isfile(settings_path):
    print("No configuration file ('.isambard_settings') found in '{}'.\nRunning configure.py...\n".format(home_dir))
    configure()
else:
    load_global_settings()


__author__ = 'Christopher W. Wood'
