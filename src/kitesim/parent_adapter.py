from .initialisation.input_classes import InputVSM, InputBridleAero
from .initialisation.yaml_loader import setup_config
from .initialisation.mutable_variables import get_mutable_variables
from .solver import solver_main
from .post_processing import post_processing_main


def module_importer():
    return (
        InputVSM,
        InputBridleAero,
        setup_config,
        get_mutable_variables,
        solver_main,
        post_processing_main,
    )
