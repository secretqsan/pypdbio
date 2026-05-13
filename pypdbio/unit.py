# pylint: disable=missing-module-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
import os
import warnings

class UnitConfig:
    def __init__(self) -> None:
        unit = os.getenv("PYPDBIO_UNIT", "A")
        self.set_unit(unit)

    def set_unit(self, unit):
        self._unit = unit
        if unit == "nm":
            self.conversion_factor = 0.1 # A to nm
        elif unit == "angstrom" or unit == "A":
            self.conversion_factor = 1.0
        else:
            warnings.warn(f"Invalid unit: {unit}, using default unit: nm")
            self.conversion_factor = 1.0

unit_config = UnitConfig()

def set_unit(unit):
    unit_config.set_unit(unit)
