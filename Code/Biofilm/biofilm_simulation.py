import biofilm_utils as utils
import biofilm_config as config

def generate_plot(geometry_config: utils.general_param, solution_config: utils.solution_param, general_config: utils.general_param) -> None:

    geometry = utils.geometry_class(geometry_config)
    solution = utils.solution_class(geometry, solution_config)

    utils.biofilm_M_scheme(solution, general_config)

if __name__ == '__main__':

    # PDE-ODE simulation
    geometry_config = config.geometry_param_2D_PDEODE_biofilm
    solution_config = config.solution_param_2D_PDEODE_biofilm
    general_config = config.general_param_2D_PDEODE_biofilm
    generate_plot(geometry_config, solution_config, general_config)

    # PDE-PDE simulation
    geometry_config = config.geometry_param_2D_PDEPDE_biofilm
    solution_config = config.solution_param_2D_PDEPDE_biofilm
    general_config = config.general_param_2D_PDEPDE_biofilm
    generate_plot(geometry_config, solution_config, general_config)