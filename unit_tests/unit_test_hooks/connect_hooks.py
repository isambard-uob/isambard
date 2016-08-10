import shutil
import os

from isambard_dev.settings import global_settings


def main(args):
    hooks_dir = os.path.join(global_settings['package_path'], '.git', 'hooks')
    if not args.unhook:
        add_hooks(hooks_dir)
        print('Finished: Hooks created, unit tests will now be run before "git commit".')
    else:
        remove_hooks(hooks_dir)
        print('Finished: Hooks deleted.')
    return


def add_hooks(dest_path):
    hook_folder = os.path.join(global_settings['package_path'], 'unit_tests', 'unit_test_hooks')
    hook_scripts = [script for script in os.listdir(hook_folder) if script[-3:] != '.py']
    for script in hook_scripts:
        print('Copying "{}" to "{}"'.format(script, dest_path))
        shutil.copy(os.path.join(hook_folder, script), dest_path)
    return


def remove_hooks(dest_path):
    hook_folder = os.path.join(global_settings['package_path'], 'unit_tests', 'unit_test_hooks')
    hook_scripts = [script for script in os.listdir(hook_folder) if script[-3:] != '.py']
    for script in hook_scripts:
        print('Removing "{}" from "{}"'.format(script, dest_path))
        try:
            os.remove(os.path.join(dest_path, script))
        except FileNotFoundError as fnfe:
            print(fnfe)
    return

if __name__ == '__main__':
    import argparse

    description = "Connects hooks to ISAMBARD git repository to enforce that unit tests must run before a commit."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-u', '--unhook', help="Removes hooks from git repository.", action="store_true")
    arguments = parser.parse_args()

    main(arguments)
