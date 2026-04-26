"""Custom Lipinski Filter with configurable thresholds and per-parameter enable/disable.
Reads filter values from ~/final_project/custom_filter_config.json.
Each parameter (MW, LogP, H-donors, H-acceptors, Rotatable Bonds) can be
individually enabled/disabled and configured with min/max ranges.
"""
import __future__
import json
import os

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.Lipinski as Lipinski
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Descriptors as Descriptors

rdkit.RDLogger.DisableLog("rdApp.*")

from autogrow.operators.filter.filter_classes.parent_filter_class import ParentFilter

# Default values (standard Lipinski Strict, all enabled)
_DEFAULTS = {
    "mw_enabled": True,
    "mw_min": 0.0,
    "mw_max": 500.0,
    "logp_enabled": True,
    "logp_min": -999.0,
    "logp_max": 5.0,
    "h_donors_enabled": True,
    "h_donors_min": 0,
    "h_donors_max": 5,
    "h_acceptors_enabled": True,
    "h_acceptors_min": 0,
    "h_acceptors_max": 10,
    "rot_bonds_enabled": True,
    "rot_bonds_min": 0,
    "rot_bonds_max": 10,
    "tpsa_enabled": False,
    "tpsa_min": 0.0,
    "tpsa_max": 140.0,
    "aromatic_rings_enabled": False,
    "aromatic_rings_min": 0,
    "aromatic_rings_max": 4,
}

# Also support old config format for backwards compatibility
_OLD_KEY_MAP = {
    "max_mw": ("mw_max", float),
    "max_logp": ("logp_max", float),
    "min_logp": ("logp_min", float),
    "max_h_donors": ("h_donors_max", int),
    "max_h_acceptors": ("h_acceptors_max", int),
    "max_rot_bonds": ("rot_bonds_max", int),
}

def _load_config():
    """Load custom filter config from JSON file."""
    config_path = os.path.expanduser("~/final_project/custom_filter_config.json")
    config = dict(_DEFAULTS)
    try:
        if os.path.exists(config_path):
            with open(config_path) as f:
                user_config = json.load(f)
            # Handle old format (max_mw, max_logp, etc.)
            for old_key, (new_key, cast) in _OLD_KEY_MAP.items():
                if old_key in user_config and new_key not in user_config:
                    user_config[new_key] = cast(user_config[old_key])
            config.update({k: v for k, v in user_config.items() if k in _DEFAULTS})
    except Exception as e:
        print(f"Warning: Could not load custom filter config: {e}")
    return config


class CustomLipinskiFilter(ParentFilter):
    """
    Custom Lipinski-style filter with user-configurable thresholds.
    Each parameter can be individually enabled/disabled and set with min/max ranges.
    Reads config from ~/final_project/custom_filter_config.json.
    """

    def run_filter(self, mol):
        """
        Filter molecule based on custom thresholds.
        Only checks parameters that are enabled in the config.

        Returns:
        :returns: bool: True if molecule passes all enabled criteria
        """
        config = _load_config()

        # Molecular Weight check
        if config.get("mw_enabled", True):
            exact_mwt = Descriptors.ExactMolWt(mol)
            if exact_mwt < config.get("mw_min", 0):
                return False
            if exact_mwt > config.get("mw_max", 500):
                return False

        # LogP check
        if config.get("logp_enabled", True):
            mol_log_p = Crippen.MolLogP(mol)
            if mol_log_p < config.get("logp_min", -999):
                return False
            if mol_log_p > config.get("logp_max", 5):
                return False

        # H-Bond Donors check
        if config.get("h_donors_enabled", True):
            num_h_donors = Lipinski.NumHDonors(mol)
            if num_h_donors < config.get("h_donors_min", 0):
                return False
            if num_h_donors > config.get("h_donors_max", 5):
                return False

        # H-Bond Acceptors check
        if config.get("h_acceptors_enabled", True):
            num_h_acceptors = Lipinski.NumHAcceptors(mol)
            if num_h_acceptors < config.get("h_acceptors_min", 0):
                return False
            if num_h_acceptors > config.get("h_acceptors_max", 10):
                return False

        # Rotatable Bonds check
        if config.get("rot_bonds_enabled", True):
            num_rot_bonds = Lipinski.NumRotatableBonds(mol)
            if num_rot_bonds < config.get("rot_bonds_min", 0):
                return False
            if num_rot_bonds > config.get("rot_bonds_max", 10):
                return False

        # TPSA check
        if config.get("tpsa_enabled", False):
            tpsa = Descriptors.TPSA(mol)
            if tpsa < config.get("tpsa_min", 0):
                return False
            if tpsa > config.get("tpsa_max", 140):
                return False

        # Aromatic Rings check
        if config.get("aromatic_rings_enabled", False):
            num_arom_rings = Descriptors.NumAromaticRings(mol)
            if num_arom_rings < config.get("aromatic_rings_min", 0):
                return False
            if num_arom_rings > config.get("aromatic_rings_max", 4):
                return False

        return True
