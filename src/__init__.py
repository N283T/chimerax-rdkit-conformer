"""ChimeraX-RDKitSMILES - Generate 3D molecules from SMILES using RDKit ETKDGv3"""

from chimerax.core.toolshed import BundleAPI


class _RdksmilesAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from . import cmd

        if command_info.name == "rdksmiles":
            from chimerax.core.commands import register

            register(command_info.name, cmd.rdksmiles_desc, cmd.rdksmiles)


bundle_api = _RdksmilesAPI()
