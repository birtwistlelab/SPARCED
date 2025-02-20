#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filename: model_factory.py
Created: 2025-02-19
Author(s): Jonah R. Huggins
Description: Factory function that selects and instantiates the appropriate model handler.
"""
import os

def create_model_handler(filepath: os.PathLike, amici_path: os.PathLike) -> "ModelHandler":
    """
    Factory function that selects and instantiates the appropriate model handler.
    """
    try:
        import amici  # check if AMICI is available
        from model_handler import AMICIModelHandler
        return AMICIModelHandler(amici_path)
    except ImportError:
        try:
            import tellurium  # check if Tellurium is available
            from model_handler import TelluriumModelHandler
            return TelluriumModelHandler(filepath)
        except ImportError:
            raise ImportError("No supported model handler library is available.")