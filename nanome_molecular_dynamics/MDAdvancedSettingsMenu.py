import re
from os import path
from functools import partial

import nanome
from .AdvancedSettings import AdvancedSettings
MENU_PATH = path.join(path.dirname(path.realpath(__file__)), "json/menus/advanced_settings.json")

class MDAdvancedSettingsMenu:
    def __init__(self, plugin, on_close):
        self.__plugin = plugin
        self.__settings = AdvancedSettings()
        self.__menu = nanome.ui.Menu.io.from_json(MENU_PATH)
        self.__menu.register_closed_callback(on_close)
        self.__tabs = {}
        self.setup_menu()

    def open(self):
        self.__plugin.menu = self.__menu
        self.__plugin.menu.enabled = True
        self.__plugin.update_menu(self.__menu)

    def switch_to_tab(self, button):
        my_category = button.name.replace(' Tab', '')
        view_name =  f'{my_category} Settings View'
        for child in self.__menu.root.find_node('Settings Views').get_children():
            settings_category = child.name.replace(' Settings View', '')
            child.enabled = child.name == f'{my_category} Settings View'
            tab = self.__tabs[settings_category]
            btn = tab.get_content()
            btn.selected = tab.name == f'{my_category} Tab'

        self.__plugin.update_menu(self.__menu)

    def render_category(self, category_name):
        category_view = self.__menu.root.find_node(f'{category_name} Settings View')
        category_view.remove_content()
        category_view.clear_children()
        category_options = AdvancedSettings.Options[category_name]
        # loop through each option
        # last_seperator = None # needed?
        for display_name in category_options:
            ln_option = nanome.ui.LayoutNode()
            category_view.add_child(ln_option)
            option = category_options[display_name]
            disabled, display_value, display_type = self.__settings.get_option_display(option)
            self.draw_option(category_name, ln_option, option, display_name, display_type, display_value, disabled)
            option['layout_node'] = ln_option

            ln_seperator = nanome.ui.LayoutNode()
            last_seperator = ln_seperator
            ln_seperator.sizing_type = nanome.util.enums.SizingTypes.ratio
            ln_seperator.sizing_value = 0.03
            ln_seperator.padding = (0.0, 0.0, 0.01, 0.01)
            seperator = ln_seperator.add_new_mesh()
            seperator.mesh_color = nanome.util.color.Color(127, 127, 127)
            category_view.add_child(ln_seperator)
        category_view.remove_child(last_seperator)

    def redraw_changed_options(self, category_name):
        category_view = self.__menu.root.find_node(f'{category_name} Settings View')
        category_options = AdvancedSettings.Options[category_name]
        for display_name in category_options:
            option = category_options[display_name]
            disabled, display_value, display_type = self.__settings.get_option_display(option)
            if disabled ^ option['layout_node'].disabled:
                option['layout_node'].disabled = disabled
                option['layout_node'].remove_content()
                option['layout_node'].clear_children()
                self.draw_option(category_name, option['layout_node'], option, display_name, display_type, display_value, disabled)

    def setup_menu(self):
        # Adds __tab_general, __tab_system, __tab_integrator, __tab_simulation buttons
        # registers callbacks for each to switch to their respective tabs
        for ln in self.__menu.root.find_node('Tabs').get_children():
            button = ln.get_content()
            tab_type = ln.name.replace(' Tab', '')
            self.__tabs[tab_type] = ln
            button.name = ln.name
            button.register_pressed_callback(self.switch_to_tab)

        # loop through each option category
        for category_name in AdvancedSettings.Options.keys():
            self.render_category(category_name)

    def update_display_value(self, choice_cell):
        content = choice_cell.get_content()
        if type(content) is nanome.ui.Button:
            content.set_all_text(str(choice_cell.display_value))
        elif type(content) is nanome.ui.TextInput:
            content.input_text = self.__settings.get_settings(option)
        elif type(content) is nanome.ui.UIList:
            for item in content.items:
                self.update_display_value(item.choice_cell)
        else:
            print('???')

        self.__plugin.update_content(content)

    def draw_option(self, category_name, ln, option, display_name, display_type, display_value, disabled, parent=None):
        # print(display_name, display_type, display_value)
        ln.layout_orientation = nanome.util.enums.LayoutTypes.horizontal
        ln.disabled = disabled

        if display_name:
            ln_label = nanome.ui.LayoutNode()
            ln_label.add_new_label(display_name)
            ln.add_child(ln_label)

        choice_cell = nanome.ui.LayoutNode()
        choice_cell.forward_dist = 0.002
        choice_cell.category_name = category_name
        choice_cell.option = option
        choice_cell.parent_cell = parent
        choice_cell.child_cells = []
        if parent: parent.child_cells.append(choice_cell)
        choice_cell.display_value = display_value
        choice_cell.display_type = display_type
        choice_cell.settings_name = display_name
        choice_cell.settings_value = self.__settings.get_setting(option)
        ln.add_child(choice_cell)
        ln.choice_cell = choice_cell

        if display_type is list or display_type is dict:
            choices_list = choice_cell.add_new_list()
            choices_list.display_rows = 2

            if display_type is list:
                for option_value in display_value:
                    ln_option = nanome.ui.LayoutNode()
                    choices_list.items.append(ln_option)
                    child_display_value = str(option_value)
                    child_display_type  = type(option_value)
                    self.draw_option(category_name, ln_option, option, '', child_display_type, child_display_value, disabled, choice_cell)
            elif display_type is dict:
                option_dict = display_value
                setting = self.__settings.get_setting(option)
                for suboption_name in option_dict:
                    ln_option = nanome.ui.LayoutNode()
                    choices_list.items.append(ln_option)
                    child_display_value = setting[suboption_name]
                    child_display_type  = type(child_display_value)
                    self.draw_option(category_name, ln_option, option, suboption_name, child_display_type, child_display_value, disabled, choice_cell)
        else:
            if display_name == 'Collision rate':
                print('Collision rate:')
                print('display value:', display_value)
            if display_type in [int, float]:
                content = choice_cell.add_new_text_input()
                content.register_changed_callback(partial(self.fix_input, choice_cell))
                content.register_submitted_callback(partial(self.handle_input_submitted, choice_cell))
                suboption_name = None if parent is None else display_name
                if disabled:
                    disabled_mesh = choice_cell.add_new_mesh()
                    disabled_mesh.mesh_color = nanome.util.color.Color(127, 127, 127)
                else:
                    content.placeholder_text = ''
                    content.input_text = self.__settings.get_setting(option) or ''
                    content.max_length = 120
            else:
                content = choice_cell.add_new_button(str(display_value))
                content.register_pressed_callback(partial(self.handle_button_pressed, choice_cell))
                # select current options
                if choice_cell.option['type'] is dict:
                    content.selected = display_value
                elif choice_cell.parent_cell:
                    content.selected = self.__settings.get_setting(option) == display_value
                content.unusable = disabled

    def fix_input(self, choice_cell, inp):
        text = inp.input_text
        data_type = choice_cell.option['type']
        if data_type is int:
            text = re.sub('[^\d]', '', text)
        elif data_type is float:
            if 'e-' in text:
                e_minus = text.index('e-') if 'e-' in text else None
                old_text = text
                text = re.sub('([^\d.])', '', text)
                text = text[:e_minus] + 'e-' + text[e_minus:]
            elif 'e' in text:
                e = text.index('e') if 'e' in text else None
                text = re.sub('([^\d.])', '', text)
                text = text[:e] + 'e' + text[e:]
            else:
                i = text.index('.') if '.' in text else None
                text = re.sub('([^\d])', '', text)
                text = text[:i] + '.' + text[i:] if i else text

        inp.input_text = text
        self.__plugin.update_content(choice_cell.get_content())

    def handle_input_submitted(self, choice_cell, inp):
        text = inp.input_text
        if text:
            setting_value = choice_cell.option['type'](text)
            # set the option value
            if choice_cell.parent_cell is None:
                choice_cell.settings_value = setting_value
            else:
                choice_cell.parent_cell.settings_value[choice_cell.settings_name] = setting_value
                setting_value = choice_cell.parent_cell.settings_value
        else:
            setting_value = None

        self.__settings.set_option(choice_cell.option, setting_value)

    def handle_button_pressed(self, choice_cell, button):
        obj_types = [list, dict]
        # if the choice_cell is for a simple option
        if choice_cell.parent_cell and choice_cell.parent_cell.display_type in obj_types:
            # if the parent is a list
            if choice_cell.parent_cell.display_type is list:
                # set the value of the parent to be the value of the child (single select)
                choice_cell.parent_cell.settings_value = choice_cell.display_value
                for ln_subcell in choice_cell.parent_cell.child_cells:
                    # None of its children are actually buttons, find the buttons another way!
                    ln_subcell.get_content().selected = False
                choice_cell.get_content().selected = True
            # otherwise if the parent is a dictionary
            elif choice_cell.parent_cell.display_type is dict:
                # toggle the value of the child in the parent (multi select)
                choice_cell.display_value = not choice_cell.display_value
                content = choice_cell.get_content()
                content.selected = choice_cell.display_value
                choice_cell.parent_cell.settings_value[choice_cell.settings_name] = choice_cell.display_value

        # identify the top choice_cell and use its value to set the option (caution: only goes up one level)
        top_cell = choice_cell.parent_cell if choice_cell.parent_cell else choice_cell
        value = top_cell.settings_value
        value = top_cell.option['type'](value) # cast the value
        self.__settings.set_option(choice_cell.option, value)
        self.update_display_value(choice_cell)
        self.redraw_changed_options(choice_cell.category_name)
        self.__plugin.update_menu(self.__menu)

