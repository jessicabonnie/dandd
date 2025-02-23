"""
DandD: Efficient measurement of sequence growth and similarity

A tool to estimate deltas for sequence sets and answer questions about relative contribution.
"""

from __future__ import annotations
# from .version import __version__
from .utils import *
from .species_specifics import SpeciesSpecifics
from .sketch_filepath import SketchFilePath
from .sketch_base import SketchObj
from .sketch_dashing import DashSketchObj
from .sketch_kmc import KMCSketchObj
from .delta_node import DeltaNode
from .delta_tree import DeltaTree, DeltaSpider, SubSpider

__all__ = [
    'DandD',
    'SketchObj',
    'DashSketchObj', 
    'KMCSketchObj',
    'SketchFilePath',
    'DeltaTree',
    'DeltaSpider',
    'SubSpider',
    'DeltaNode'
]
