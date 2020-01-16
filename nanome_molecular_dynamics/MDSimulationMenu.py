import nanome
from os import path

MENU_PATH = path.join(path.dirname(path.realpath(__file__)), "json", "menus", "main.json")

class MDSimulationMenu():
    def __init__(self, plugin):
        self.__plugin = plugin
        self.__btn_refresh = None
        self.__btn_run = None
        self.__complexes_list = None

        self.__menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self.__menu.index = 1
        self.setup_menu()

    def open(self):
        self.__plugin.menu = self.__menu
        self.__menu.enabled = True
        self.__plugin.update_menu(self.__menu)

    def change_state(self, state):
        self.__btn_run.set_all_text('Stop' if state else 'Start')
        self.__plugin.update_content(self.__btn_run)

    def change_complex_list(self, complex_list):
        def complex_pressed(button):
            complex = button.complex
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button")
                return

            if button.selected == False:
                button.selected = True
                self.__plugin._selected_complexes.append(complex.index)
            else:
                button.selected = False
                self.__plugin._selected_complexes.remove(complex.index)
            self.__plugin.update_content(button)

        self.__plugin._selected_complexes = []
        self.__complexes_list.items = []

        for complex in complex_list:
            clone = self.__complex_item_prefab.clone()
            btn = clone.get_content()
            btn.set_all_text(complex.molecular.name)
            btn.complex = complex
            btn.register_pressed_callback(complex_pressed)
            self.__complexes_list.items.append(clone)

        self.__plugin.update_menu(self.__menu)

    def remove_progress_bar(self, update=False):
        self.__progress_bar.enabled = False
        if update:
            self.__plugin.update_menu(self.__menu)

    def set_progress_bar(self, title, percentage, description, update=False):
        self.__progress_bar.enabled = True
        if title: self.__progress_bar.get_content().title = title
        if percentage: self.__progress_bar.get_content().percentage = percentage
        if description: self.__progress_bar.get_content().description = description
        if update: self.__plugin.update_menu(self.__menu)

    def disable_run(self, update=False):
        self.__btn_run.unusable = True
        if update: self.__plugin.update_menu(self.__menu)

    def enable_run(self, update=False):
        self.__btn_run.unusable = False
        if update: self.__plugin.update_menu(self.__menu)

    def setup_menu(self):
        def run_button_pressed_callback(button):
            self.__plugin.toggle_simulation()

        ln_list = self.__menu.root.find_node('List')

        # Set the complex list
        self.__complexes_list = ln_list.get_content()

        # Set the refresh and run buttons
        self.__btn_run = self.__menu.root.find_node("StartStop").get_content()
        self.__btn_run.register_pressed_callback(run_button_pressed_callback)

        # Set status bar
        self.__progress_bar = self.__menu.root.find_node("Progress Bar")
        self.__progress_bar.add_new_loading_bar()

        # Set the entry list prefab
        self.__complex_item_prefab = self.__menu.root.find_node('Entry Prefab')