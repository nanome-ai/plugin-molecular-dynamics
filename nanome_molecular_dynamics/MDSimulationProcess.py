from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.element import *
from simtk.openmm.app.internal.unitcell import computePeriodicBoxVectors
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat

import subprocess
from timeit import default_timer as timer

import nanome
from nanome.util import Logs

# TEMP
from nanome._internal._structure._bond import _Bond
from nanome._internal._structure._io._pdb.save import Options as PDBOptions

nanometer = nano * meter
picosecond = pico * second
nb_steps = 100
metalElements = ['Al','As','Ba','Ca','Cd','Ce','Co','Cs','Cu','Dy','Fe','Gd','Hg','Ho','In','Ir','K','Li','Mg',
        'Mn','Mo','Na','Ni','Pb','Pd','Pt','Rb','Rh','Sm','Sr','Te','Tl','V','W','Yb','Zn']

pdb_options = PDBOptions()
pdb_options.write_bonds = True

class MDSimulationProcess():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        self.__reporter = MDReporter(self)

    @staticmethod
    def get_bond_type(kind):
        if kind == _Bond.Kind.CovalentSingle:
            return Single
        if kind == _Bond.Kind.CovalentDouble:
            return Double
        if kind == _Bond.Kind.CovalentTriple:
            return Triple
        return None

    @staticmethod
    def get_atom_symbol(name, atoms_nb):
        upper = name.upper()
        if upper.startswith('CL'):
            return chlorine
        elif upper.startswith('NA'):
            return sodium
        elif upper.startswith('MG'):
            return magnesium
        elif upper.startswith('BE'):
            return beryllium
        elif upper.startswith('LI'):
            return lithium
        elif upper.startswith('K'):
            return potassium
        elif upper.startswith('ZN'):
            return zinc
        elif (atoms_nb == 1 and upper.startswith('CA')):
            return calcium
        else:
            return Element.getBySymbol(upper[0])

    def set_stream(self, stream):
        self.__stream = stream

    def delete_alternate_atoms(self, topology, positions):
        modeller = Modeller(topology, positions)
        delete_atoms = []
        for chain in topology.chains():
            for indexInChain, residue in enumerate(chain.residues()):
                atom_names = []
                for atom in residue.atoms():
                    if atom.name in atom_names:
                        delete_atoms.append(atom)
                    else:
                        atom_names.append(atom.name)

        modeller.delete(delete_atoms)
        return (modeller.getTopology(), modeller.getPositions())

    def fix_complexes(self, complex_list):
        fixed_complexes = []
        for complex in complex_list:
            for residue in complex.residues:
                atoms = residue._atoms
                for i in range(len(atoms) - 1, -1, -1):
                    if atoms[i].molecular.is_het == True:
                        del atoms[i]

            complex.io.to_pdb("tmp.pdb", pdb_options)

            fixer = PDBFixer(filename="tmp.pdb")
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.removeHeterogens(False)
            fixer.addMissingHydrogens(7.0)

            (topology, positions) = self.delete_alternate_atoms(fixer.topology, fixer.positions)
            with open('tmp2.pdb', 'w') as pdb_file:
                PDBFile.writeFile(topology, positions, pdb_file)

            fixed_complex = nanome.structure.Complex.io.from_pdb(path="tmp2.pdb")
            fixed_complex.index = complex.index
            fixed_complex.position = complex.position
            fixed_complex.rotation = complex.rotation
            fixed_complex.molecular.name = complex.molecular.name
            fixed_complex.rendering.visible = True
            fixed_complexes.append(fixed_complex)

        return fixed_complexes

    def init_simulation(self, complex_list):
        # Create topology
        topology = Topology()
        added_atoms = dict()
        positions = []
        PDBFile._loadNameReplacementTables()
        self.__complex_list = complex_list
        min_x = max_x = min_y = max_y = min_z = max_z = None
        Logs.debug("Create topology")
        for complex in complex_list:
            for molecule in complex.molecules:
                for chain in molecule.chains:
                    sim_chain = topology.addChain()
                    for residue in chain.residues:
                        residueName = residue.molecular.name
                        if residueName in PDBFile._atomNameReplacements:
                            atomReplacements = PDBFile._atomNameReplacements[residueName]
                        else:
                            atomReplacements = {}
                        sim_residue = topology.addResidue(residue.molecular.name, sim_chain)
                        for atom in residue.atoms:
                            molecular = atom.molecular
                            symbol = MDSimulationProcess.get_atom_symbol(molecular.name, len(residue._atoms))
                            atom_name = molecular.name
                            if atom_name in atomReplacements:
                                atom_name = atomReplacements[atom_name]
                            sim_atom = topology.addAtom(atom_name, symbol, sim_residue)
                            added_atoms[atom.index] = sim_atom
                            position = molecular.position
                            positions.append(Vec3(position.x * 0.1 * nanometer, position.y * 0.1 * nanometer, position.z * 0.1 * nanometer))
                            if min_x == None or position.x < min_x:
                                min_x = position.x
                            if max_x == None or position.x > max_x:
                                max_x = position.x
                            if min_y == None or position.y < min_y:
                                min_y = position.y
                            if max_y == None or position.y > max_y:
                                max_y = position.y
                            if min_z == None or position.z < min_z:
                                min_z = position.z
                            if max_z == None or position.z > max_z:
                                max_z = position.z

        topology.createStandardBonds()
        topology.createDisulfideBonds(positions)
        added_bonds = set(topology.bonds())
        for complex in complex_list:
            for molecule in complex.molecules:
                for chain in molecule.chains:
                    for residue in chain.residues:
                        for bond in residue.bonds:
                            if bond.index in added_bonds:
                                continue
                            atom1 = added_atoms[bond.atom1.index]
                            atom2 = added_atoms[bond.atom2.index]
                            type = MDSimulationProcess.get_bond_type(bond.molecular.kind)
                            topology.addBond(atom1, atom2, type)
                            added_bonds.add(bond.index)

        # topology.setPeriodicBoxVectors(computePeriodicBoxVectors(max_x - min_x, max_y - min_y, max_z - min_z, 90, 90, 90))
        topology.setPeriodicBoxVectors(computePeriodicBoxVectors(49.163, 45.981, 38.869, 90.00, 90.00, 90.00))

        # Create simulation parameters
        # nonbondedMethod = PME
        [templates, residues] = self.__forcefield.generateTemplatesForUnmatchedResidues(topology)
        for index, residue in enumerate(residues):
            template = templates[index]
            print("unmatched residue:", residue)
        for index, template in enumerate(templates):
            residue = residues[index]
            residue_bonds = list(residue.internal_bonds())+list(residue.external_bonds())
            print("residue bonds:", residue_bonds)
            (unique_res_bonds, unique_tmp_bonds) = ForceField.findMissingBonds(residue, template)
            print("UNIQUE RESIDUE BONDS:", unique_res_bonds)
            print("UNIQUE TEMPLATE BONDS:", unique_tmp_bonds)
            print(f"RESIDUE {residue.name}:")
            if len(unique_res_bonds) > 0:
                print("Missing bonds:")
                for res_bond in unique_res_bonds:
                    print(f"\n{res_bond}")
            print(f"TEMPLATE {template.name}:")
            if len(unique_tmp_bonds) > 0:
                print("Missing bonds:")
                for tmp_bond in unique_tmp_bonds:
                    print(f"\n{tmp_bond}")
            print("-------------------------")
            for bond in unique_res_bonds:
                template.name += '[' + ForceField.getAtomID(bond[0]) + '<-->' + ForceField.getAtomID(bond[1]) + ']'
                if bond[0] in residue.atoms() and bond[1] in residue.atoms():
                    template.addBondByName(bond[0].name, bond[1].name)
                else:
                    template.addExternalBondByName(bond[0].name)

            if template.name not in self.__forcefield._templates:
                self.__forcefield.registerResidueTemplate(template)
                print(f"registering template {template.name}")
            else:
                print(f"redundant template {template.name} ********************")

        system = self.__forcefield.createSystem(topology, nonbondedMethod = NoCutoff, nonbondedCutoff = 1 * nanometer, constraints = HBonds)

        # Set the simulation
        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picosecond)
        simulation = Simulation(topology, system, integrator)
        # Set reporting
        simulation.reporters.append(self.__reporter)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        self.__simulation = simulation

        self.__positions = [0.0] * ((len(positions) * 3) + 5000000)

    def simulate(self, complex_list):
        # positions = []
        # for complex in complex_list:
        #     for molecule in complex.molecules:
        #         for chain in molecule.chains:
        #             for residue in chain.residues:
        #                 for atom in residue.atoms:
        #                     position = atom.molecular.position
        #                     positions.append(Vec3(position.x * 0.1, position.y * 0.1, position.z * 0.1))

        # self.__simulation.context.setPositions(positions)

        self.__start = timer()
        self.__simulation.step(nb_steps)

    def simulation_result(self, positions):
        end = timer()
        Logs.debug("Simulation:", end - self.__start)
        self.__start = timer()
        new_positions = self.__positions
        for i in range(len(positions)):
            position = positions[i]
            x = position[0]._value * 10
            if math.isnan(x):
                Logs.warning("Got a NaN value, ignoring it")
                continue
            y = position[1]._value * 10
            if math.isnan(y):
                Logs.warning("Got a NaN value, ignoring it")
                continue
            z = position[2]._value * 10
            if math.isnan(z):
                Logs.warning("Got a NaN value, ignoring it")
                continue

            new_positions[i * 3] = x
            new_positions[i * 3 + 1] = y
            new_positions[i * 3 + 2] = z
        self.__stream.update(new_positions, self.__plugin.on_simulation_done)
        end = timer()
        Logs.debug("Sent new complexes:", end - self.__start)

# This class is a reporter for OpenMM Simulation class
class MDReporter(object):
    def __init__(self, process):
        self.__process = process

    def describeNextReport(self, simulation):
        return (nb_steps, True, False, False, False, None)

    def report(self, simulation, state):
        self.__process.simulation_result(state.getPositions())