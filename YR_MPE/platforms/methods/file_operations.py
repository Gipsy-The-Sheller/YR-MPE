"""
File Operations Module
统一的文件路径生成、目录创建、文件保存逻辑
"""
import os
import shutil
from pathlib import Path
from datetime import datetime
from typing import Union, Optional


class FileOperations:
    """文件操作工具类"""
    
    @staticmethod
    def create_workspace_subdirs(workspace_path: Union[str, Path]) -> bool:
        """
        创建workspace标准子目录结构
        
        Args:
            workspace_path: Workspace根目录路径
            
        Returns:
            bool: 创建是否成功
        """
        try:
            workspace_path = Path(workspace_path)
            subdirs = ['alignments', 'models', 'trees', 'distances', 'temp', 'history']
            for subdir in subdirs:
                (workspace_path / subdir).mkdir(exist_ok=True)
            return True
        except Exception as e:
            print(f"Failed to create workspace subdirs: {e}")
            return False
    
    @staticmethod
    def generate_filename(operation: str, extension: str, prefix: str = "") -> str:
        """
        生成标准化的文件名
        
        Args:
            operation: 操作类型（如 'alignment', 'model', 'tree'）
            extension: 文件扩展名（如 '.fasta', '.json', '.nwk'）
            prefix: 可选前缀
            
        Returns:
            str: 生成的文件名
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if prefix:
            filename = f"{prefix}_{operation}_{timestamp}{extension}"
        else:
            filename = f"{operation}_{timestamp}{extension}"
        return filename
    
    @staticmethod
    def save_file_to_workspace(workspace_path: Union[str, Path], 
                             file_type: str, 
                             filename: str, 
                             content: Union[str, bytes],
                             binary: bool = False) -> str:
        """
        将文件保存到workspace指定目录
        
        Args:
            workspace_path: Workspace根目录路径
            file_type: 文件类型目录 ('alignments', 'models', 'trees', 'distances')
            filename: 文件名
            content: 文件内容
            binary: 是否为二进制文件
            
        Returns:
            str: 完整文件路径
        """
        workspace_path = Path(workspace_path)
        file_path = workspace_path / file_type / filename
        
        mode = 'wb' if binary else 'w'
        encoding = None if binary else 'utf-8'
        
        with open(file_path, mode, encoding=encoding) as f:
            f.write(content)
        
        return str(file_path)
    
    @staticmethod
    def copy_file_to_workspace(workspace_path: Union[str, Path], 
                              file_type: str, 
                              source_path: Union[str, Path],
                              new_filename: Optional[str] = None) -> str:
        """
        复制文件到workspace指定目录
        
        Args:
            workspace_path: Workspace根目录路径
            file_type: 文件类型目录
            source_path: 源文件路径
            new_filename: 新文件名（可选，如果不提供则使用原文件名）
            
        Returns:
            str: 目标文件路径
        """
        workspace_path = Path(workspace_path)
        source_path = Path(source_path)
        
        if new_filename is None:
            new_filename = source_path.name
            
        target_path = workspace_path / file_type / new_filename
        shutil.copy2(source_path, target_path)
        
        return str(target_path)