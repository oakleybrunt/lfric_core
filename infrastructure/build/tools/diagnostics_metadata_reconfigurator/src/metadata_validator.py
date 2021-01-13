##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
"""
Validation step for diagnostic fields

Checks that all fields have either a long or a standard name
"""

from logging import getLogger, WARNING, ERROR

from diagnostics_metadata_collection import Metadata

LOGGER = getLogger("reconfigurator.validator")


class InvalidMetadataError(ValueError):
    """Error showing field(s) have invalid metadata"""


def validate_metadata(metadata: Metadata, force: bool = False) -> None:
    """Validator that reports any diagnostic field(s) with invalid metadata"""
    invalid_fields = []
    LOGGER.info("Checking if fields have long or standard name")

    for field in metadata.get_fields():
        LOGGER.debug("Checking %s", field.unique_id)

        if not (field.long_name or field.standard_name):
            LOGGER.log(WARNING if force else ERROR,
                       "Field '%s' missing long and standard names",
                       field.unique_id)
            invalid_fields.append(field.unique_id)

    if invalid_fields and not force:
        raise InvalidMetadataError(f"Field(s) {invalid_fields} missing long "
                                   f"and standard names")
    else:
        LOGGER.info("Metadata valid")
