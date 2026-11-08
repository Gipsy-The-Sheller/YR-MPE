from ..templates.base_plugin_ui import BasePlugin
from PyQt5.QtWidgets import QVBoxLayout, QLabel

class HelloPluginUI(BasePlugin):
    def __init__(self, parent=None):
        super().__init__(parent)

    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        for i, color in zip(range(100, 30, -12), 
                            ["#4CAF50", "#8BC34A",
                             "#CDDC39", "#FFEB3B",
                             "#FFC107"]):
            label = QLabel("Hello, MPE!")
            label.setStyleSheet(f"""font-size: {i}px; 
                                font-weight: bold; 
                                color: {color};""")
            layout.addWidget(label)
    def run(self):
        return self