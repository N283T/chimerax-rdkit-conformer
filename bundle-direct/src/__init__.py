"""ChimeraX-RDKitConformer - Generate 3D conformers from molecular notations using RDKit ETKDGv3"""

from chimerax.core.toolshed import BundleAPI


class _RdkconfAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from . import cmd
        from chimerax.core.commands import register

        if command_info.name == "rdkconf":
            register(command_info.name, cmd.rdkconf_desc, cmd.rdkconf)


bundle_api = _RdkconfAPI()
