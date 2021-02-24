###############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
from enum import Enum


class StandardSynonyms(str, Enum):
    # validated standards
    CF = "CF"
    CMIP6 = "CMIP6"
    # non validated Synonyms
    AMIP = "AMIP"
    GRIB = "GRIB"
    STASH = "STASH"
