"""
TCGS-SEQUENTION Protocol
========================

Three-Gate Falsification Protocol for Synchronous Parallel Emergence Testing.

Author: Henry Arellano-Peña
License: CC BY 4.0
"""

__version__ = "3.0.0"
__author__ = "Henry Arellano-Peña"
__email__ = "harellano@unal.edu.co"

from .protocol import ThreeGateProtocol, ProtocolConfig
from .ltee_analysis import run_full_protocol

__all__ = [
    "ThreeGateProtocol",
    "ProtocolConfig", 
    "run_full_protocol",
]
