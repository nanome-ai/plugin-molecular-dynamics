import nanome
import os
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

class AdvancedSettings:
    # general options
    general_forcefield_name_values = ['amber14-all', 'amber10', 'amber03', 'amber99sbnmr', 'amber99sbildn', 'amber99sb', 'amber96']
    general_water_model_name_values = ['spce', 'tip3pfb', 'tip3p', 'tip4pew', 'tip4pfb', 'tip5p', 'amber10_obc']
    general_platform_name_values = ['Reference', 'OpenCL', 'CPU', 'CUDA']
    general_platform_precisions = ['single', 'mixed', 'double']
    general_platform_precision_names = [None, 'OpenCLPrecision', None, 'CudaPrecision']
    general_platform_precision_values = [None, general_platform_precisions, None, general_platform_precisions]
    general_platform_device_index_names = [None, 'OpenCLDeviceIndex', None, 'CudaDeviceIndex']
    general_platform_device_index_values = [None, int, None, int]
    general_platform_index_names = [None, 'OpenCLPlatformIndex', None, None]
    general_platform_index_values= [None, int, None, None]

    general_options = {
        'Forcefield': {
            'value': 'general_forcefield_name_values',
            'depends': None,
            'default': 'amber14-all',
            'type': str
        },
        'Water Model': {
            'value': 'general_water_model_name_values',
            'depends': None,
            'default': 'tip3pfb',
            'type': str
        },
        'Platform': {
            'value': 'general_platform_name_values',
            'depends': None,
            'default': 'CPU',
            'type': str
        },
        'Precision': {
            'name': 'general_platform_precision_names',
            'value': 'general_platform_precision_values',
            'depends': 'general_platform_name',
            'default': 'single',
            'type': str
        },
        'Device index': {
            'name': 'general_platform_device_index_names',
            'value': 'general_platform_device_index_values',
            'depends': 'general_platform_name',
            'default': None,
            'type': int
        },
        'Platform index': {
            'name': 'general_platform_index_names',
            'value': 'general_platform_index_values',
            'depends': 'general_platform_name',
            'default': None,
            'type': int
        }
    }

    # system options
    system_nonbonded_method_names = 'nonbondedMethod'
    system_nonbonded_method_values = ['NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'Ewald', 'PME']
    system_ewald_error_tolerance_names = [None, None, None, 'ewaldErrorTolerance', 'ewaldErrorTolerance']
    system_ewald_error_tolerance_values = [None, None, None, float, float]
    system_contraint_names = 'constraints'
    system_constraint_values = ['None', 'HBonds', 'AllBonds', 'HAngles']
    integrator_contraint_error_tolerance_values = [None, float, float, float]
    system_rigid_water_names = 'rigidWater'
    system_rigid_water_values = [True, False]
    system_nonbonded_cutoff_names = 'nonbondedCutoff'
    system_nonbonded_cutoff_values = [None, float, float, float, float]
    system_barostats = ['None', 'Monte Carlo']
    system_barostat_values = [system_barostats, None, system_barostats, system_barostats, None]
    system_pressure_values = [None, float]
    system_barostat_interval_values = [None, int]
    system_thermostats = ['None', 'Andersen']
    system_thermostat_values = [None, system_thermostats, None, None, system_thermostats]
    system_random_init_vel_values = [True, False]
    system_generation_temp_values = [float, None]

    system_options = {
        'Nonbonded method' : {
            'name': 'system_nonbonded_method_names',
            'value': 'system_nonbonded_method_values',
            'depends': None,
            'default': 'CutoffPeriodic',
            'type': str
        },
        'Ewald error tolerance' : {
            'name': 'system_ewald_error_tolerance_names',
            'value': 'system_ewald_error_tolerance_values',
            'depends': 'system_nonbonded_method',
            'default': 0.0005,
            'type': float
        },
        'Constraints': {
            'name': 'system_contraint_names',
            'value': 'system_constraint_values',
            'depends': None,
            'default': 'HBonds',
            'type': str
        },
        'Constraint error tol.': {
            'value': 'integrator_contraint_error_tolerance_values',
            'depends': 'system_constraint',
            'default': 0.00001,
            'type': float
        },
        'Rigid water': {
            'name': 'system_rigid_water_names',
            'value': 'system_rigid_water_values',
            'depends': None,
            'default': True,
            'type': bool
        },
        'Nonbonded cutoff': {
            'name': 'system_nonbonded_cutoff_names',
            'value': 'system_nonbonded_cutoff_values',
            'depends': 'system_nonbonded_method',
            'default': 1.0,
            'unit': unit.nanometers,
            'type': float
        },
        'Random init vels.': { #TODO: Get dependency working
            'value': 'system_random_init_vel_values',
            'depends': None,
            'default': True,
            'type': bool
        },
        'Generation temp.': {
            'value': 'system_generation_temp_values',
            'depends': 'system_random_init_vel',
            'default': 300.0,
            'type': float
        }
    }

    # integrator options
    integrator_integrator_values = ['Langevin', 'Verlet', 'Brownian', 'VariableLangevin', 'VariableVerlet']
    integrator_timestep_values = [float, float, float, None, None]
    integrator_error_tolerance_values = [None, None, None, float, float]
    integrator_collision_rate_values = [float, None, float, None, None]
    integrator_temperature_values = [float, None, float, None, None]

    integrator_options = {
        'Integrator' : {
            'value': 'integrator_integrator_values',
            'depends': None,
            'default': 'Langevin',
            'type': str
        },
        'Timestep' : {
            'name': 'timestep',
            'value': 'integrator_timestep_values',
            'depends': 'integrator_integrator',
            'default': 2.0,
            'unit': unit.femtosecond,
            'type': float
        },
        'Error tolerance' : {
            'value': 'integrator_error_tolerance_values',
            'depends': 'integrator_integrator',
            'default': 0.0001,
            'type': float
        },
        'Collision rate' : {
            'value': 'integrator_collision_rate_values',
            'depends': 'integrator_integrator',
            'default': 1.0,
            'unit': unit.picosecond,
            'type': float
        },
        'Temperature' : {
            'value': 'integrator_temperature_values',
            'depends': 'integrator_integrator',
            'default': 300,
            'unit': unit.kelvin,
            'type': float
        },
        'Barostat': {
            'value': 'system_barostat_values',
            'depends': 'integrator_integrator',
            'default': 'None',
            'type': str
        },
        'Pressure': {
            'value': 'system_pressure_values',
            'depends': 'system_barostat',
            'default': 2.5,
            'unit': unit.atmosphere,
            'type': float
        },
        'Barostat Interval': {
            'value': 'system_barostat_interval_values',
            'depends': 'system_barostat',
            'default': 25,
            'type': int
        },
        'Thermostat': {
            'value': 'system_thermostat_values',
            'depends': 'integrator_integrator',
            'default': 'None',
            'type': str
        }
    }

    # simulation options
    simulation_reporter_interval_values = int
    simulation_reporter_options_values = {'Positions': bool, 'Velocities': bool, 'Forces': bool, 'Energies': bool}
    simulation_equilibrium_steps_values = int
    simulation_production_steps_values = int
    simulation_minimize_values = [True, False]
    simulation_max_minimization_steps_values = [int, None]

    simulation_options = {
        'Report Interval': {
            'value': 'simulation_reporter_interval_values',
            'depends': None,
            'default': 200,
            'type': int
        },
        'Reporter options': {
            'value': 'simulation_reporter_options_values',
            'depends': None,
            'default': {'Positions': True, 'Velocities': False, 'Forces': False, 'Energies': False},
            'type': dict
        },
        'Equilibration steps': {
            'value': 'simulation_equilibrium_steps_values',
            'depends': None,
            'default': 1000,
            'type': int
        },
        'Minimize?': {
            'value': 'simulation_minimize_values',
            'depends': None,
            'default': True,
            'type': bool
        },
        'Max minimize steps': {
            'value': 'simulation_max_minimization_steps_values',
            'depends': 'simulation_minimize',
            'default': None,
            'type': int
        }
    }

    Options = {'General' : general_options, 'System': system_options, 'Integrator' : integrator_options, 'Simulation': simulation_options}
    instance = None

    def __init__(self):
        AdvancedSettings.instance = self

        # general settings
        self.general_pdb_name = 'tmp.pdb'
        self.general_pdb = None
        self.general_forcefield_name = None
        self.general_water_model_name = None
        self.general_platform_name = None
        self.general_platform_precision = None
        self.general_platform_device_index = None
        self.general_platform_index = None
        self.general_platform = None
        self.topology = None
        self.platform = None
        self.forcefield = None
        # system settings
        self.system_nonbonded_method = None
        self.system_ewald_error_tolerance = None
        self.system_constraint = None
        self.system_contraint_error_tolerance = None
        self.system_rigid_water = None
        self.system_nonbonded_cutoff = None
        self.system_barostat = None
        self.system_barostat_interval = None
        self.system_thermostat = None
        self.system_random_init_vel = None
        self.system_generation_temp = None
        self.system = None
        # integrator settings
        self.integrator_integrator = None
        self.integrator_collision_rate = None
        self.integrator = None
        # simulation settings
        self.simulation_reporter_interval = None
        self.simulation_reporter_options = None
        self.simulation_equilibrium_steps = None
        self.simulation_production_steps = None
        self.simulation_minimize = None
        self.simulation_max_minimization_steps = None
        self.simulation = None
        # TODO: update
        self._reporter = None

        for category_name in AdvancedSettings.Options:
            self.set_to_defaults(AdvancedSettings.Options[category_name])

    def reset(self):
        self.platform = None
        self.forcefield = None
        self.system = None
        self.integrator = None
        self.simulation = None

    def get_setting(self, option, unit=False):
        setting_value = getattr(self, option['value'].replace('_values', ''), option['default'])
        if option.get('unit') and unit:
            setting_value = setting_value * option['unit']
        return setting_value

    def get_dependency_index(self, option, name):
        dependency = getattr(self, name, '')
        dependency_values = getattr(AdvancedSettings, f'{name}_values', None)
        first_dep_index, first_dep_value = next(((i, item) for i, item in enumerate(dependency_values) if item is not None), (0, None))
        if type(dependency) in dependency_values:
            return dependency_values.index(type(dependency))
        elif dependency in dependency_values:
            return dependency_values.index(dependency)
        else:
            return first_dep_index

    def get_option_display_type(self, option):
        dep_name = option['depends']
        display_values = getattr(AdvancedSettings, option['value'])
        if dep_name:
            display_values = next((item for item in display_values if item is not None), None)

        if type(display_values) is type:
            return display_values
        else:
            return type(display_values)

    def get_option_display(self, option):
        dependency_name = option['depends']
        display_values = getattr(AdvancedSettings, option['value'])
        display_type = self.get_option_display_type(option)
        disabled = False
        if dependency_name is None:
            if display_type is list:
                display_value = display_values
            else:
                display_value = self.get_setting(option)
        elif '.' in dependency_name:
            # get the dependency dictionary name and key
            (dep_dict_name, key) = dependency_name.split('.')
            # get the dependency dictionary
            dep_dict = getattr(self, dep_dict_name)
            bool_index = dep_dict[key]
            disabled = display_values[int(bool_index)] is None
            display_value = self.get_setting(option)
        else:
            dep_index = self.get_dependency_index(option, dependency_name)
            display_value = display_values[dep_index]
            disabled = display_value is None
            if disabled:
                display_value = next((val for val in display_values if val is not None), 0)

        return disabled, display_value, display_type

    def get_option_default(self, option, suboption_name=None):
        if suboption_name:
            return option['default'][suboption_name]
        else:
            return option['default']

    def get_option_type(self, option):
        return option['type']

    def get_option_by_name(self, options, option_name):
        option = options[option_name]
        return self.get_option(option)

    def get_option(self, option, unit=True):
        # ensure the option is named
        names = getattr(AdvancedSettings, option.get('name') or '', None)
        values = getattr(AdvancedSettings, option['value'], None) or option['default']
        name = None
        value = None

        dependency_name = option['depends']
        value = getattr(self, option['value'].replace('_values', ''), None)
        if dependency_name is None:
            # simply lookup the name and value
            name = names
        elif '.' in dependency_name:
            # if the dependency contains a subfield
            # assume its a bool-valued dictionary
            (dep_dict_name, key) = dependency_name.split('.')
            # get the dependency dictionary
            dep_dict = getattr(self, dep_dict_name)
            # lookup the dependency which indicates an index in our option's values array
            bool_index = dep_dict[key]
            # if values is None at the dependency index, assign value to None
            value = None if values[int(bool_index)] is None else value
        else:
            # lookup the name and value by dependency
            dep_index = self.get_dependency_index(option, dependency_name)
            name = names[dep_index] if names else None
            if type(values[dep_index]) not in [type, list]:
                value = values[dep_index]

        if unit and value is not None:
            value = value * option['unit'] if 'unit' in option else value

        return (name, value)

    def get_options(self, options_type, option_list, kind=dict):
        all_options = AdvancedSettings.Options[options_type]
        for display_name in option_list:
            name, value = self.get_option_by_name(all_options, display_name)
            if kind is dict:
                if name is not None and value is not None:
                    yield (name, value)
            elif kind is list:
                if value is not None:
                    yield value

    def get_named_options(self, option_type):
        options = AdvancedSettings.Options[option_type]
        named_options = {}
        if options is not None:
            for option in options.keys():
                (name, value) = self.get_option_by_name(options, option)
                # if ∃ property (defined or ∃ default) and should ∃ (according to dependency)
                if name is not None and value is not None:
                    # set it
                    named_options[name] = value

        return named_options

    def set_option(self, option, value):
        """Sets an option to a value
            options -- an options dictionary
            option  -- the the option to set
            value   -- the value to set the option to
        """
        var_name = option['value'].replace('_values', '')
        setattr(self, var_name, value)

    def set_to_default(self, options, option):
        """Sets an option to its default value
            options -- an options dictionary
            option  -- the option to set
        """
        var_name = option['value'].replace('_values', '')
        value = option['default']
        setattr(self, var_name, value)

    def set_to_defaults(self, category_options):
        """Sets all options in the category in an option dictionary to the default values
            options -- an options dictionary
        """
        for option_name in category_options.keys():
            self.set_to_default(category_options, category_options[option_name])

    def get_pdb(self):
        if self.pdb == None:
            self.pdb = app.PDBFile(self.pdb_name)
            self.topology = self.pdb.topology
            self.positions = self.pdb.positions
        return self.pdb

    def get_forcefield(self):
        if self.forcefield == None:
            self.forcefield = app.ForceField(f'{self.general_forcefield_name}.xml', f'{self.general_water_model_name}.xml')
        return self.forcefield

    def get_platform(self):
        if self.platform == None:
            self.platform = mm.Platform.getPlatformByName(self.general_platform_name)
        return self.platform

    def get_platform_properties(self):
        return self.get_named_options('General')

    def get_system(self, topology=None, residueTemplates=None):
        if topology:
            self.topology = topology

        if self.system is None:
            # get system args
            args = self.get_named_options('System')
            for name, arg in args.items():
                if type(arg) is str:
                    args[name] = getattr(app, arg, None)
            self.system = self.get_forcefield().createSystem(self.topology, residueTemplates=residueTemplates, **args)
        return self.system

    def get_integrator(self):
        if self.integrator == None:
            params = list(self.get_options('Integrator', ['Error tolerance', 'Temperature', 'Collision rate', 'Timestep'], list))
            self.integrator = getattr(mm, f'{self.integrator_integrator}Integrator', None)(*params)
            if self.integrator_integrator == 'VariableLangevin' or self.integrator_integrator == 'VariableVerlet':
                self.integrator.setConstraintTolerance(params[0])

        return self.integrator

    def get_simulation(self, positions=None):
        if positions:
            self.positions = positions

        if not self.simulation:
            self.simulation = app.Simulation(self.topology, self.get_system(), self.get_integrator(), self.get_platform(), self.get_platform_properties())
            self.simulation.context.setPositions(self.positions)
        return self.simulation

    def attach_reporter(self, reporter_class, result_callback):
        reporter = reporter_class(self, result_callback)
        self.get_simulation().reporters = [reporter]