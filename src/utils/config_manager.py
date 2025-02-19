#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script name: config_manager.py
Created: 02/19/2025
Author: Jonah R. Huggins

Description: Handles loading, storing, and sending configuration file information.
"""
# -----------------------Package Import & Defined Arguements-------------------#
import os
import yaml

class ConfigLoader:
    """
    Handles loading, storing, and sending configuration file information.
    """
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.config = self._load_config()

    def _load_config(self):
        """Reads YAML file and returns a dictionary."""
        try:
            with open(self.filepath, encoding = 'utf-8', mode = 'r') as file:
                return yaml.safe_load(file) or {}  # Ensures no `NoneType`
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {self.filepath}")
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML file: {e}")

    def get(self, key, default=None):
        """Fetches a key from the config dictionary using dot notation."""
        keys = key.split(".")
        value = self.config
        for k in keys:
            value = value.get(k, default) if isinstance(value, dict) else default
        return value
    