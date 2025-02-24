#!/bin/bash python3
"""
Filename: threshold_handler.py
Created: 2025-02-23
Author(s): Jonah R. Huggins

Description: File handles model-specific (user-defined) threshold rules for 
            how SPARCED should handle certain model components.
"""
import re

class ThresholdHandler:
    """Handles model-specific (user-defined) threshold rules for how SPARCED should handle certain model components."""
    def __init__(self, model, exceptions):
        self.model = model
        self.exceptions = exceptions


    def _handle_thresholds(self):
        """Evaluates user-defined threshold rules dynamically."""

        threshold_rules = self.exceptions.thresholds  # Load thresholds from config

        for rule_name, expression in threshold_rules.items():

            # Extract species/variables from expression (e.g., ['cPARP', 'PARP'])
            species_in_rule = re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', expression)

            # Get their actual values from the model
            species_values = {species: self.model.get_species_concentration(species) for species in species_in_rule}

            # Convert the expression into an evaluable string
            eval_expression = expression
            for species, value in species_values.items():
                eval_expression = eval_expression.replace(species, str(value))

            # Evaluate the threshold condition
            threshold_met = eval(eval_expression)  # SAFE: all components are numerical values

            print(f"Threshold '{rule_name}': {expression} → {eval_expression} = {threshold_met}")

            if threshold_met:
                self._trigger_event(rule_name)  # Call some function if threshold is met

    def _trigger_event(self, rule_name):
        """Triggers an event based on a threshold rule."""
        print(f"Threshold '{rule_name}' met! Triggering event...")

        # Do something based on the threshold rule
        pass