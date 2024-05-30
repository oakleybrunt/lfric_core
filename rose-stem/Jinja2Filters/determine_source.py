##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Function to read through the dependencies.sh in an lfric_apps working copy and
determine the specified source for the argument, repo. Uses a subprocess ssh
command as working copy may be remote when this needs to be run by cylc.
Intended as a Jinja2 Custom Filter
'''

import os
import subprocess

def determine_source(wc_loc, repo):
    '''
    Main function for Jinja2 Custom Filter
    Get a source for a particular repo from the dependencies.sh file in the
    source for this rose-stem.
    '''
    # Determine the host and the path to the dependencies file, and populate the
    # subprocess command
    try:
        host, path = wc_loc.split(":")
        path = os.path.join(path, "dependencies.sh")
        cat_command = f"ssh {host} cat {path}"
    except ValueError:
        path = os.path.join(wc_loc, "dependencies.sh")
        cat_command = f"cat {path}"
    # Run the subprocess command and read the stdout
    result = subprocess.run(
        cat_command.split(), capture_output=True, text=True
    )
    dependencies_file = result.stdout.split("\n")
    source = ''
    rev = ''

    # Read through the dependencies file and populate revision and source
    # variables for requested repo
    for line in dependencies_file:
        line = line.strip()
        if line.startswith(f"export {repo}_rev"):
            rev = line.split("=")[1]
        if line.startswith(f"export {repo}_sources"):
            source = line.split("=")[1]
    # If source not set then default to trunk
    if source == '':
        # lfric_core doesn't match the url
        if repo == "lfric_core":
            source = "fcm:lfric.xm_tr"
        else:
            source = f"fcm:{repo}.xm_tr"
    # If a revision set then append to source
    # Defaults to the head of the source
    if rev != '':
        source = f"{source}@{rev}"
    return source
