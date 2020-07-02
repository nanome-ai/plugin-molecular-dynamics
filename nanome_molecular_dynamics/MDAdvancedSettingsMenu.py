import re
from os import path
from functools import partial
from nanome.util import Logs
import nanome
from .AdvancedSettings import AdvancedSettings
from nanome.api.ui import Dropdown,DropdownItem

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
        self.__settings.reset()
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
        option_counter = 0
        for display_name in category_options:
            option = category_options[display_name]
            disabled, display_value, display_type = self.__settings.get_option_display(option)
            default_value = self.__settings.get_option_default(option)
            ln_option = nanome.ui.LayoutNode()
            ln_option.sizing_type = 1
            ln_option.set_size_fixed(0.001)
            ln_option.set_padding(0,0,0.01,0.01)
            category_view.add_child(ln_option)
            
            if display_type is dict:
                ln_option.set_size_fixed(0.001*len(display_value))
                option_counter += len(display_value)-1
            self.draw_option(category_name, ln_option, option, display_name, display_type, display_value, disabled, default = default_value)
            option['layout_node'] = ln_option
            option_counter += 1

            # ln_seperator = nanome.ui.LayoutNode()
            # last_seperator = ln_seperator
            # ln_seperator.sizing_type = nanome.util.enums.SizingTypes.ratio
            # ln_seperator.sizing_value = 0.03
            # ln_seperator.padding = (0.0, 0.0, 0.01, 0.01)
            # seperator = ln_seperator.add_new_mesh()
            # seperator.mesh_color = nanome.util.color.Color(127, 127, 127)
            # category_view.add_child(ln_seperator)
        # category_view.remove_child(last_seperator)
        while option_counter < 9:
            ln_option = nanome.ui.LayoutNode()
            ln_option.sizing_type = 1
            ln_option.set_size_fixed(0.001)
            ln_option.set_padding(0,0,0.01,0.01)
            category_view.add_child(ln_option)
            option_counter += 1
        

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
            #content.value.set_all(str(choice_cell.display_value))
        elif type(content) is nanome.ui.TextInput:
            content.input_text = self.__settings.get_settings(option)
        elif type(content) is nanome.ui.UIList:
            for item in content.items:
                self.update_display_value(item.choice_cell)
        elif type(content) is nanome.ui.Dropdown:
            #nothing to be displayed
            pass
        else:
            print('???')

        self.__plugin.update_content(content)

    def draw_option(self, category_name, ln, option, display_name, display_type, display_value, disabled, parent=None, default = None):
        ln.layout_orientation = nanome.util.enums.LayoutTypes.horizontal
        ln.disabled = disabled
        if display_name:
            ln_label_wrapper = nanome.ui.LayoutNode()
            ln_label = nanome.ui.LayoutNode()
            ln_label.add_new_label(display_name)
            label_text = ln_label.get_content()
            if parent == None:
                label_text._text_auto_size = False
                label_text._text_size = 0.3
            else:
                label_text._text_auto_size = False
                label_text._text_size = 0.25
            label_text._text_bold = True
            label_text._text_vertical_align = 1
            ln_label_wrapper.add_child(ln_label)
            ln.add_child(ln_label_wrapper)

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
            
            if display_type is list:
               
                # choices_list = choice_cell.add_new_list()
                # choices_list.display_rows = 2

                # for option_value in display_value:
                #     ln_option = nanome.ui.LayoutNode()
                #     choices_list.items.append(ln_option) 
                #     child_display_value = option_value
                #     child_display_type  = type(option_value)
                #     self.draw_option(category_name, ln_option, option, '', child_display_type, child_display_value, disabled, choice_cell)
                choices_dropdown = choice_cell.add_new_dropdown()
                choice_cell.forward_dist = 0.003
                dropdown_list = []
                for i,option_value in enumerate(display_value):
                    dd_item = DropdownItem()
                    dd_item._name = str(option_value)
                    #dd_item._selected = self.__settings.get_setting(option) == default
                    dd_item._selected = dd_item._name == str(default)
                    dropdown_list.append(dd_item)
                choices_dropdown.items = dropdown_list
                choices_dropdown.register_item_clicked_callback(partial(self.handle_dropdown_pressed,choice_cell))
                if display_name and display_name == 'Forcefield':
                    choices_dropdown._max_displayed_items = 4


            elif display_type is dict:
                choices_list = choice_cell.add_new_list()
                choices_list.display_rows = len(display_value)

                option_dict = display_value
                setting = self.__settings.get_setting(option)
                if display_type is dict and display_name:
                    label_text._text_vertical_align = 1
                    for _ in range(len(display_value)):
                        ln_label_wrapper.add_child(nanome.ui.LayoutNode())
                for suboption_name in option_dict:
                    ln_option = nanome.ui.LayoutNode()
                    choices_list.items.append(ln_option)
                    child_display_value = setting[suboption_name]
                    child_display_type  = type(child_display_value)
                    self.draw_option(category_name, ln_option, option, suboption_name, child_display_type, child_display_value, disabled, choice_cell)
        else:
            if display_type in [int, float]:
                text_input_wrapper = nanome.ui.LayoutNode()
                text_input_wrapper.set_padding(0,0,0.01,0.01)
                text_input_wrapper._forward_dist = 0.002
                content = text_input_wrapper.add_new_text_input()
                choice_cell.add_child(text_input_wrapper)
                content.register_changed_callback(partial(self.fix_input, choice_cell))
                content.register_submitted_callback(partial(self.handle_input_submitted, choice_cell))
                suboption_name = None if parent is None else display_name
                if disabled:
                    disabled_mesh = text_input_wrapper.add_new_mesh()
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
                else:
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

    def handle_dropdown_pressed(self,choice_cell,dropdown,item):
        
        value = item._name
        # bool("False") = True
        if value == 'False':
            value = False
        else:
            value = choice_cell.option['type'](value)
        self.__settings.set_option(choice_cell.option,value)
        self.update_display_value(choice_cell)
        self.redraw_changed_options(choice_cell.category_name)
        self.__plugin.update_menu(self.__menu)
        