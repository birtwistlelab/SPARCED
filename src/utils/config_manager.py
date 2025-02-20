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
import json

class ConfigManager:
    """
    Handles loading, storing, and sending configuration file information.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.config = {}
        self.config_loader = self.get_config_loader(filepath)
        self.config = self.config_loader.config

    def load_config(self):
        """Placeholder method to be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement this method.")

    def get(self, key, default=None):
        """Fetch a nested config value using dot notation."""
        keys = key.split(".")
        value = self.config
        for k in keys:
            value = value.get(k, default) if isinstance(value, dict) else default
        return value

class YAMLConfigLoader(ConfigManager):
    def __init__(self, filepath):
        super().__init__(filepath)  # Call parent constructor
        self.load_config()  # Load config upon initialization

    def load_config(self):
        """Loads configuration from a YAML file."""
        try:
            with open(self.filepath, encoding = 'utf-8', mode = 'r') as file:
                self.config = yaml.safe_load(file) or {}
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {self.filepath}")
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML file: {e}")
        
class JSONConfigLoader(ConfigManager):
    def __init__(self, filepath):
        super().__init__(filepath)  # Call parent constructor
        self.load_config()  # Load config upon initialization

    def load_config(self):
        """Loads configuration from a JSON file."""
        try:
            with open(self.filepath, encoding = 'utf-8', mode = 'r') as file:
                self.config = json.load(file) or {}
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {self.filepath}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Error parsing JSON file: {e}")

@staticmethod
def get_config_loader(filepath):
    """Factory method to get the appropriate config loader based on file extension."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".yaml":
        return YAMLConfigLoader(filepath)
    elif ext == ".json":
        return JSONConfigLoader(filepath)
    else:
        raise ValueError(f"Unsupported file extension: {ext}")
    
# -----------------------------End of config_manager.py-------------------------#
