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
import json
import yaml

class ConfigManager:
    """
    Handles loading, storing, and sending configuration file information.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.config = {}
        self.config_loader = get_config_loader(filepath)
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
    """
    Loads configuration from a YAML file.
    """
    def __init__(self, filepath):
        super().__init__(filepath)  # Call parent constructor
        self.load_config()  # Load config upon initialization

    def load_config(self):
        """Loads configuration from a YAML file."""
        try:
            with open(self.filepath, encoding = 'utf-8', mode = 'r') as file:
                self.config = yaml.safe_load(file) or {}
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"Configuration file not found: {self.filepath}") from exc
        except yaml.YAMLError as exc:
            raise ValueError(f"Error parsing YAML file: {exc}") from exc
        
class JSONConfigLoader(ConfigManager):
    """
    Loads configuration from a JSON file.
    """
    def __init__(self, filepath):
        super().__init__(filepath)  # Call parent constructor
        self.load_config()  # Load config upon initialization

    def load_config(self):
        """Loads configuration from a JSON file."""
        try:
            with open(self.filepath, encoding = 'utf-8', mode = 'r') as file:
                self.config = json.load(file) or {}
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"Configuration file not found: {self.filepath}") from exc
        except json.JSONDecodeError as exc:
            raise ValueError(f"Error parsing JSON file: {exc}") from exc

@staticmethod
def get_config_loader(filepath):
    """Factory method to get the appropriate config loader based on file extension."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".yaml":
        return YAMLConfigLoader(filepath)
    if ext == ".json":
        return JSONConfigLoader(filepath)

    raise ValueError(f"Unsupported file extension: {ext}")

# -----------------------------End of config_manager.py-------------------------#
