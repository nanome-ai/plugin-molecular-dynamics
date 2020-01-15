import nanome
from nanome.util import Logs
from nanome.util.stream import StreamCreationError

import traceback

from timeit import default_timer as timer

from .MDSimulationMenu import MDSimulationMenu
from .MDSimulationProcess import MDSimulationProcess
from .MDAdvancedSettingsMenu import MDAdvancedSettingsMenu

class MDSimulation(nanome.PluginInstance):
    def start(self):
        self._selected_complexes = []
        self.__running = False
        self.__waiting_for_complexes = False
        self.__first_request = True
        self.__stream = None

        self.__menu = MDSimulationMenu(self)
        self.__menu.open()
        self.request_refresh()

        self.__advanced_settings_menu = MDAdvancedSettingsMenu(self, self.open_main_menu)
        self.__simulation = MDSimulationProcess(self)

    def update(self):
        if self.__running and not self.__waiting_for_complexes:
            self.__waiting_for_complexes = True
            if self.__first_request == True:
                self.request_complexes(self._selected_complexes, self.on_complexes_received)
                self.__first_request = False
            else:
                print('not first request...running simulation')
                self.__start = timer()
                # self.request_complex_list(self.on_complex_list_received)
                self.on_complex_list_received(None)

    def on_run(self):
        if self.__running == False:
            self.start_simulation()
        else:
            self.stop_simulation()

    def on_advanced_settings(self):
        self.__advanced_settings_menu.open()

    def open_main_menu(self, menu):
        self.__menu.open()

    def start_simulation(self):
        Logs.debug("Start Simulation")
        self.__start = timer()
        self.__running = True
        self.__waiting_for_complexes = False
        self.__first_request = True
        self.__stream = None
        self.__menu.change_state(True)
        # self.__menu.set_progress_bar("Loading", 10, "Converting input to mol2", update=True)

    def stop_simulation(self):
        Logs.debug("Stop Simulation")
        self.__running = False
        # self.__menu.remove_progress_bar()
        self.__menu.change_state(False)
        if self.__stream != None:
            self.__stream.destroy()

    def toggle_simulation(self):
        if (self.__running):
            self.stop_simulation()
        else:
            self.start_simulation()

    def on_complex_list_received(self, complex_list):
        if self.__running == False:
            self.__menu.change_complex_list(complex_list)
        elif self.__waiting_for_complexes == True:
            end = timer()
            Logs.debug("Request:", end - self.__start)
            self.__run_simulation(False, complex_list)

    def request_refresh(self):
        if self.__running == True:
            return
        self.request_complex_list(self.on_complex_list_received)

    def on_complex_added(self):
        self.request_refresh()

    def on_complex_removed(self):
        self.request_refresh()

    def on_complexes_received(self, complex_list):
        self.add_bonds(complex_list, self.bonds_added)

    def bonds_added(self, complex_list):
        complex_list = self.__simulation.fix_complexes(complex_list)
        self.__complex_list = complex_list

        # if settings.duplicateSystem:
            # for complex in complex_list:
            #     for molecule in complex.molecules:
            #         molecule.index = -1
            #         for chain in molecule.chains:
            #             chain.index = -1
            #             for residue in chain.residues:
            #                 residue.index = -1
            #                 for atom in residue.atoms:
            #                     atom.index = -1
            #                 for bond in residue.bonds:
            #                     bond.index = -1

        self.update_structures_deep(complex_list, self._complexes_updated)

    def _complexes_updated(self):
        end = timer()
        Logs.debug("First Request:", end - self.__start)
        self._start = timer()
        self.request_complexes(self._selected_complexes, self._updated_complexes_received)

    def _updated_complexes_received(self, complex_list):
        self.__complex_list = complex_list
        indices = []
        for complex in complex_list:
            for atom in complex.atoms:
                indices.append(atom.index)
        end = timer()
        Logs.debug("Second Request:", end - self.__start)
        self._start = timer()
        self.create_stream(indices, self.on_stream_creation)

    def on_stream_creation(self, stream, error):
        if error == StreamCreationError.AtomNotFound:
            Logs.error("Tried to create a stream with bad atoms")
            self.stop_simulation()
            return

        self.__stream = stream
        self.__simulation.set_stream(stream)
        end = timer()
        Logs.debug("Stream creation:", end - self.__start)
        self.__run_simulation(True, self.__complex_list)

    def __run_simulation(self, init, complex_list):
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                if init == True:
                    self.__start = timer()
                    self.__simulation.init_simulation(complex_list)
                    end = timer()
                    Logs.debug("Init simulation:", end - self.__start)
                self.__simulation.simulate(complex_list)
                return
            except:
                if attempt >= 3:
                    Logs.error("Got an error", attempt, "times, aborting simulation:")
                    Logs.error(traceback.format_exc())
                    self.stop_simulation()

    def on_simulation_done(self):
        self.__waiting_for_complexes = False

def main():
    plugin = nanome.Plugin("MD Simulation", "Run molecular dynamics on the selected complexes, using OpenMM", "MD", True)
    plugin.set_plugin_class(MDSimulation)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
