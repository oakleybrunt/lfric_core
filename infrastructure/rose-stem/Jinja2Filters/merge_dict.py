#!/usr/bin/env python
##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
A Jinja2 filter to combine two dictionaries.
'''


def merge_dict(context, new_dict):
    '''
    Takes a dictionary and replaces key/value pairs with those from the new
    dictionary.

    @param[in] context Dictionary
    @param[in] new_dict Dictionary
    @return Dictionary combined from inputs.
    '''
    if new_dict is not None:
        if context is None:
            context = {}
        for key, value in new_dict.items():
            context[key] = value
    return context
