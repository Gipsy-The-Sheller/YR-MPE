class DatasetItem:
    """Dataset数据项模型"""
    
    def __init__(self):
        self.selected = False          # 是否选中
        self.loci_name = ""           # 位点名称
        self.length = 0               # 序列长度
        self.sequence_count = 0       # 序列数量
        self.is_aligned = False       # 是否已比对
        self.file_path = ""           # 原始文件路径
        self.sequences = []           # 序列数据
        self.name = ""                # 新增name属性
        
    def __str__(self):
        return f"DatasetItem(loci_name={self.loci_name}, length={self.length}, count={self.sequence_count})"