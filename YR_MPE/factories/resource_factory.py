import os
from PyQt5.QtGui import QIcon


class ResourceFactory:
    """统一资源管理工厂"""
    
    def __init__(self, base_path=None):
        if base_path is None:
            # 获取YR_MPE包的根目录
            self.base_path = os.path.dirname(os.path.dirname(__file__))
        else:
            self.base_path = base_path
    
    def get_icon(self, icon_path):
        """获取图标资源
        Args:
            icon_path: 相对于YR_MPE/icons/的路径
                    例如: "file/dataset.svg", "software/iqtree.svg"
        """
        full_path = os.path.join(self.base_path, "icons", icon_path)
        return QIcon(full_path)
    
    def get_icon_path(self, icon_path):
        """获取图标文件完整路径"""
        return os.path.join(self.base_path, "icons", icon_path)