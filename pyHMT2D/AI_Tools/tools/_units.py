# -*- coding: utf-8 -*-
"""Unit resolution shared by the AI tools.

Both SRH-2D and HEC-RAS projects come in English or SI flavours, and pyHMT2D
reports lengths, areas, and discharges in whatever unit the project itself
uses — it never converts. Tools therefore have to label their output from the
open project instead of assuming metres, which is what they used to do.
"""

from __future__ import annotations

# 1 ft = 0.3048 m exactly, hence 1 ft² = 0.09290304 m².
_SQFT_TO_SQM = 0.3048 ** 2

# SRH-2D writes the OutputUnit keyword as "EN" or "SI".
_ENGLISH_SRH = {"EN", "ENGLISH"}
# RAS_2D_Data.units is set to 'Feet', 'Meter', or 'Unknown units'.
_ENGLISH_RAS = {"FEET", "FOOT", "ENGLISH"}


def is_english(session) -> bool:
    """Return True when the open project uses English units.

    Falls back to SI whenever the project is closed or does not state its
    units — that preserves the historical default, and only labels (plus the
    km² conversion) depend on it.
    """
    data = getattr(session, "data", None)
    if data is None:
        return False

    if getattr(session, "model_type", None) == "SRH-2D":
        try:
            _, unit = data.get_output_format_unit()
        except Exception:
            return False
        return str(unit).strip().upper() in _ENGLISH_SRH

    return str(getattr(data, "units", "")).strip().upper() in _ENGLISH_RAS


def length_unit(session) -> str:
    """'ft' or 'm' — the project's length unit."""
    return "ft" if is_english(session) else "m"


def area_unit(session) -> str:
    """'ft2' or 'm2' — the project's area unit."""
    return "ft2" if is_english(session) else "m2"


def discharge_unit(session) -> str:
    """'ft3/s' or 'm3/s' — the project's discharge unit."""
    return "ft3/s" if is_english(session) else "m3/s"


def area_to_km2(area: float, session) -> float:
    """Convert an area given in the project's own units to km²."""
    return area * (_SQFT_TO_SQM if is_english(session) else 1.0) / 1e6
