"""Loads a json file containing settings that the ISAMBARD use.

A settings.json file must be made in this directory, template_settings.py can be copied and paths for your system
can be added. Json files look like python dictionaries, and in fact all json files are valid python dictionary code
(although this isn't necessarily true in reverse). The settings.json file is ignored by git to it avoid being
overwritten with another systems paths. When the isambard is imported a isambard.settings dictionary is generated, which
contains the user chosen settings, but it can be explicitly loaded by an internal module to get these values too.
"""

import json as _json
import os as _os

with open(_os.path.dirname(_os.path.abspath(__file__)) + '/settings.json', 'r') as _settings_f:
    global_settings = _json.loads(_settings_f.read())
global_settings[u'package_path'] = _os.path.dirname(_os.path.abspath(__file__), )

__author__ = 'Christopher W. Wood'
