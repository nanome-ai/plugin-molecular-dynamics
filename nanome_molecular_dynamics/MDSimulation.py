import os
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
        nanome._internal._plugin_instance.UPDATE_RATE = 1.0/200

        self._selected_complexes = []
        self.running = False
        self.__stream = None

        self.__menu = MDSimulationMenu(self)
        self.__menu.open()
        self.request_refresh()

        self.__advanced_settings_menu = MDAdvancedSettingsMenu(self, self.open_main_menu)
        self.__process = MDSimulationProcess(self)

    def on_run(self):
        if self.running == False:
            self.start_simulation()
        else:
            self.stop_simulation()

    def on_stop(self):
        self.stop_simulation()

    def on_advanced_settings(self):
        self.__advanced_settings_menu.open()

    def open_main_menu(self, menu):
        self.__menu.open()

    def start_simulation(self):
        Logs.debug("Start Simulation")
        self.__start = timer()
        self.running = True
        self.__stream = None
        self.__menu.change_state(True)
        self.request_complexes(self._selected_complexes, self.on_complexes_received)
        # self.__menu.set_progress_bar("Loading", 10, "Converting input to mol2", update=True)

    def stop_simulation(self):
        Logs.debug("Stop Simulation")
        self.running = False
        # self.__menu.remove_progress_bar()
        self.__menu.change_state(False)
        if self.__stream != None:
            self.__stream.destroy()
        try:
            os.remove('tmp.pdb')
        except:
            pass

    def toggle_simulation(self):
        if self.running:
            self.stop_simulation()
        else:
            self.start_simulation()

    def on_complex_list_received(self, complex_list):
        self.__menu.change_complex_list(complex_list)

    def request_refresh(self):
        self.request_complex_list(self.on_complex_list_received)

    def on_complex_added(self):
        self.request_refresh()

    def on_complex_removed(self):
        self.request_refresh()

    def on_complexes_received(self, complex_list):
        # complex_list = self.__process.fix_complexes(complex_list)
        # self.add_bonds(complex_list, self.bonds_added)
        self.bonds_added(complex_list)

    def bonds_added(self, complex_list):
        self.__complex_list = complex_list
        self.update_structures_deep(complex_list, self._complexes_updated)

    def bonds_added_after_fixed(self, complex_list):
        self.__complex_list = complex_list

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
        self.create_writing_stream(indices, nanome.util.enums.StreamType.position, self.on_stream_creation)

    def on_stream_creation(self, stream, error):
        self.__waiting_for_complexes = False
        if error == StreamCreationError.AtomNotFound:
            Logs.error("Tried to create a stream with bad atoms")
            self.stop_simulation()
            return

        self.__stream = stream
        self.__process.set_stream(stream)
        end = timer()
        Logs.debug("Stream creation:", end - self.__start)
        self.__run_simulation(self.__complex_list)

    def __run_simulation(self, complex_list):
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                self.__start = timer()
                self.__process.init_simulation(complex_list)
                self.send_notification(nanome.util.enums.NotificationTypes.message, "Simulating...")
                self.__process.simulate(complex_list)
                return
            except:
                if attempt >= 3:
                    Logs.error("Got an error", attempt, "times, aborting simulation:")
                    Logs.error(traceback.format_exc())
                    self.stop_simulation()

def main():
    plugin = nanome.Plugin("MD Simulation", "Run molecular dynamics on the selected complexes, using OpenMM", "MD", True)
    plugin.set_plugin_class(MDSimulation)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
