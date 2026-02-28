def form_calibration_table(data):
    """
    form_calibration_table: Form a calibration file for LSD2 based on the input data
    
    :param data: list of calibrations: 
    {
    'name': str:name, 
    'set': list:taxon_set, 
    'type': str: 'fixed' / 'interval' / 'upper' / 'lower',
    'values': list:values
    }
    """

    switch_type = {
        'fixed': lambda x: f'{x[0]}',
        'interval': lambda x: f'b({x[0]},{x[1]})',
        'upper': lambda x: f'u({x[0]})',
        'lower': lambda x: f'l({x[0]})'
    }

    calibrations = f'{len(data)}'
    for calibration in data:
        if len(calibration['set']) > 1:
            taxa_str = ','.join(taxon for taxon in calibration['set'])
            calibrations = f'{calibrations}\nmrca({taxa_str})'
        elif len(calibration['set']) == 1:
            calibrations = f'{calibrations}\n{calibration["set"][0]}'
        else:
            raise ValueError('Empty taxon set')
        
        calibrations = f'{calibrations} {switch_type[calibration["type"]](calibration["values"])}'
    
    return calibrations

def form_lsd2_command(sequence_length = -1, findroot = False, tipdating = False):
    command = []
    if sequence_length > 0:
        command += ['-s', str(sequence_length)]
    
    if findroot:
        command += ['-r', 'a']
    else:
        command += ['-r', 'k']
    
    if not tipdating:
        command += ['-d', 0]
    
    return lambda tree_file, date_file: ['-i', tree_file, '-o', date_file] + command


def form_calibration_table_with_tips(node_calibrations, tip_calibrations):
    """
    生成包含 Node 和 Tip 标定的 LSD2 校准文件
    
    :param node_calibrations: Node 校准点（MRCA），列表格式：
        [{
            'name': str:name, 
            'set': list:taxon_set, 
            'type': str: 'fixed' / 'interval' / 'upper' / 'lower',
            'values': list:values
        }, ...]
    :param tip_calibrations: Tip 校准点（单个 OTU），字典格式：
        {
            otu_name: {
                'type': str: 'fixed' / 'interval' / 'upper' / 'lower',
                'values': list:values
            } or None,
            ...
        }
    :return: 校准文件内容字符串
    """
    
    switch_type = {
        'fixed': lambda x: f'{x[0]}',
        'interval': lambda x: f'b({x[0]},{x[1]})',
        'upper': lambda x: f'u({x[0]})',
        'lower': lambda x: f'l({x[0]})'
    }

    calibrations = []
    
    # 添加 Tip 校准（单个 OTU）
    if tip_calibrations:
        for otu_name, cal_data in tip_calibrations.items():
            if cal_data is not None:
                calibrations.append((otu_name, cal_data['type'], cal_data['values']))
    
    # 添加 Node 校准（MRCA）
    if node_calibrations:
        for calibration in node_calibrations:
            if not calibration:
                continue  # 跳过 None 或空校准点
            
            if len(calibration['set']) > 1:
                taxa_str = ','.join(taxon for taxon in calibration['set'])
                calibrations.append((f'mrca({taxa_str})', calibration['type'], calibration['values']))
            elif len(calibration['set']) == 1:
                calibrations.append((calibration['set'][0], calibration['type'], calibration['values']))
            else:
                raise ValueError('Empty taxon set')
    
    # 构建校准文件内容
    if not calibrations:
        return "0"  # 没有校准点
    
    result = f'{len(calibrations)}'
    for taxa, cal_type, values in calibrations:
        result += f'\n{taxa} {switch_type[cal_type](values)}'
    
    return result