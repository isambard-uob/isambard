# The AMPAL Framework

This folder contains code for the AMPAL framwork, which is used to represent biomolecular structure in ISAMBARD. AMPAL classes for new biomolecules should be added here to a file specifically for that molecule type. If a new molecule class is added, a corresponding filter should be added to the PDBParser.

The `secondary_structure_elements` and `parameterisations` modules are also included in this directory as they directly inherit from the AMPAL framework.
