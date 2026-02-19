"""
数据项按钮组件
支持颜色编码和单击/双击事件处理的QToolButton
"""
from PyQt5.QtWidgets import QToolButton, QMenu, QAction
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QSize
from PyQt5.QtGui import QIcon

from .dataset_models import (
    DatasetItem, DatasetInfo,
    SELECTION_STATE_NONE, SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE, SELECTION_STATE_RED
)


# 按钮样式定义（只改变背景颜色，保持其他属性不变）
BUTTON_STYLES = {
    SELECTION_STATE_NONE: """
        QToolButton {
            background-color: transparent;
        }
    """,
    SELECTION_STATE_GREEN: """
        QToolButton {
            background-color: #90EE90;
            border-radius: 6px;
        }
    """,
    SELECTION_STATE_BLUE: """
        QToolButton {
            background-color: #87CEEB;
            border-radius: 6px;
        }
    """,
    SELECTION_STATE_RED: """
        QToolButton {
            background-color: #FFB6C1;
            border-radius: 6px;
        }
    """
}

# 图标映射（根据数据类型）
ICON_MAPPING = {
    "sequence": "align.svg",
    "alignment": "align.svg",
    "model": "model.svg",
    "distance": "dist.svg",
    "phylogeny": "phylogeny.svg",
    "variant": "variants.svg",
    "coalescent": "coalescent.svg",
    "clock": "clock.svg"
}


class DatasetItemButton(QToolButton):
    """数据项按钮 - 支持颜色编码和单击/双击事件"""
    
    # 信号定义
    clicked_single = pyqtSignal(str)  # 单击事件，传递item_id
    clicked_double = pyqtSignal(str)  # 双击事件，传递item_id
    right_clicked = pyqtSignal(str)   # 右键点击事件，传递item_id
    
    def __init__(self, item: DatasetItem, parent=None):
        super().__init__(parent)
        
        self.item = item
        self.click_timer = QTimer()
        self.click_timer.setSingleShot(True)
        self.click_timer.timeout.connect(self._on_single_click_delayed)
        
        self.click_count = 0
        self.last_click_time = 0
        self.pending_single_click = False
        
        # 初始化按钮
        self._init_button()
        
    def _init_button(self):
        """初始化按钮"""
        # 不设置文本，只显示图标
        # self.setText(self._get_button_text())  # 移除文本
        
        # 设置工具按钮样式为只显示图标
        self.setToolButtonStyle(Qt.ToolButtonIconOnly)
        
        # 设置固定尺寸（与原有按钮一致）
        self.setFixedSize(45, 45)
        self.setIconSize(QSize(45, 45))
        
        # 应用样式
        self.update_style()
        
        # 连接点击事件
        self.clicked.connect(self._on_click)
        
        # 设置上下文菜单
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        
    def _get_button_text(self) -> str:
        """获取按钮文本"""
        # 获取类型标签
        type_label = self.item.item_type.capitalize()
        
        # 获取名称（截断过长的名称）
        name = self.item.get_name()
        if len(name) > 15:
            name = name[:12] + "..."
        
        # 添加额外信息
        extra_info = ""
        if self.item.item_type == "model":
            model_name = self.item.data.get("model", "")
            if model_name:
                extra_info = f"\n{model_name}"
        elif self.item.item_type == "distance":
            method = self.item.data.get("method", "ML")
            extra_info = f"\n{method}"
        
        return f"{type_label}\n{name}{extra_info}"
    
    def update_style(self):
        """更新按钮样式"""
        style = BUTTON_STYLES.get(
            self.item.selection_state,
            BUTTON_STYLES[SELECTION_STATE_NONE]
        )
        self.setStyleSheet(style)
    
    def update_item(self, item: DatasetItem):
        """更新数据项"""
        self.item = item
        self.setText(self._get_button_text())
        self.update_style()
    
    def _on_click(self):
        """处理点击事件"""
        self.click_count += 1
        
        if self.click_count == 1:
            # 第一次点击，等待可能的第二次点击
            self.pending_single_click = True
            self.click_timer.start(250)  # 等待250毫秒
        elif self.click_count == 2:
            # 第二次点击，是双击
            self.click_timer.stop()
            self.pending_single_click = False
            self.clicked_double.emit(self.item.id)
            self.click_count = 0
    
    def _on_single_click_delayed(self):
        """延迟的单击事件处理"""
        if self.pending_single_click and self.click_count == 1:
            # 确认是单击（没有第二次点击）
            self.clicked_single.emit(self.item.id)
            self.pending_single_click = False
            self.click_count = 0
    
    def _show_context_menu(self, position):
        """显示右键菜单"""
        menu = QMenu(self)
        
        # 查看详情
        view_action = QAction("View Details", self)
        view_action.triggered.connect(lambda: self.clicked_double.emit(self.item.id))
        menu.addAction(view_action)
        
        # 如果是选中状态，添加取消选择选项
        if self.item.selection_state == SELECTION_STATE_GREEN:
            menu.addSeparator()
            deselect_action = QAction("Deselect", self)
            deselect_action.triggered.connect(self._deselect)
            menu.addAction(deselect_action)
        
        # 添加删除选项
        menu.addSeparator()
        delete_action = QAction("Delete", self)
        delete_action.triggered.connect(self._delete_item)
        menu.addAction(delete_action)
        
        menu.exec_(self.mapToGlobal(position))
    
    def _deselect(self):
        """取消选择"""
        self.clicked_single.emit(self.item.id)
    
    def _delete_item(self):
        """删除数据项"""
        # 这个信号需要连接到管理器的删除方法
        # TODO: 实现删除逻辑
        pass
    
    # 移除了mouseDoubleClickEvent，因为双击已经在_on_click中处理
    
    def enterEvent(self, event):
        """鼠标进入事件"""
        super().enterEvent(event)
        # 可以在这里添加悬停效果
        self.setToolTip(self._get_tooltip())
    
    def _get_tooltip(self) -> str:
        """获取工具提示"""
        tooltip = f"<b>{self.item.item_type.capitalize()}</b><br>"
        tooltip += f"Name: {self.item.get_name()}<br>"
        tooltip += f"Created: {self.item.created_at.strftime('%Y-%m-%d %H:%M:%S')}<br>"
        
        if self.item.selection_state != SELECTION_STATE_NONE:
            tooltip += f"<br><b>Status:</b> {self.item.selection_state}<br>"
            if self.item.selection_reason:
                tooltip += f"<b>Reason:</b> {self.item.selection_reason}<br>"
        
        if not self.item.is_valid:
            tooltip += f"<br><span style='color: red'><b>Errors:</b><br>"
            tooltip += "<br>".join(self.item.validation_errors)
            tooltip += "</span>"
        
        return tooltip


class DatasetButton(QToolButton):
    """数据集按钮"""
    
    clicked_single = pyqtSignal(str)  # 单击事件
    clicked_double = pyqtSignal(str)  # 双击事件
    
    def __init__(self, dataset: DatasetInfo, parent=None):
        super().__init__(parent)
        
        self.dataset = dataset
        self.click_timer = QTimer()
        self.click_timer.setSingleShot(True)
        self.click_timer.timeout.connect(self._on_single_click_delayed)
        
        self.click_count = 0
        self.pending_single_click = False
        
        # 不设置文本，只显示图标
        # self.setText(self._get_button_text())  # 移除文本
        
        # 设置工具按钮样式为只显示图标
        self.setToolButtonStyle(Qt.ToolButtonIconOnly)
        
        # 设置固定尺寸（与原有按钮一致）
        self.setFixedSize(45, 45)
        self.setIconSize(QSize(45, 45))
        
        # 应用样式
        self.update_style()
        
        # 连接点击事件
        self.clicked.connect(self._on_click)
    
    def update_style(self):
        """更新按钮样式"""
        style = BUTTON_STYLES.get(
            self.dataset.selection_state,
            BUTTON_STYLES[SELECTION_STATE_NONE]
        )
        self.setStyleSheet(style)
    
    def _on_click(self):
        """处理点击事件"""
        self.click_count += 1
        
        if self.click_count == 1:
            # 第一次点击，等待可能的第二次点击
            self.pending_single_click = True
            self.click_timer.start(250)  # 等待250毫秒
        elif self.click_count == 2:
            # 第二次点击，是双击
            self.click_timer.stop()
            self.pending_single_click = False
            self.clicked_double.emit(self.dataset.id)
            self.click_count = 0
    
    def _on_single_click_delayed(self):
        """延迟的单击事件处理"""
        if self.pending_single_click and self.click_count == 1:
            # 确认是单击（没有第二次点击）
            self.clicked_single.emit(self.dataset.id)
            self.pending_single_click = False
            self.click_count = 0
    
    def _get_button_text(self) -> str:
        """获取按钮文本（已弃用，保留以兼容）"""
        text = f"{self.dataset.name}"
        
        if self.dataset.is_multigene:
            text += f"\n{self.dataset.loci_count} loci"
        
        return text
    
    # 移除了mouseDoubleClickEvent，因为双击已经在_on_click中处理
    
    def update_dataset(self, dataset: DatasetInfo):
        """更新数据集"""
        self.dataset = dataset
        self.setText(self._get_button_text())
