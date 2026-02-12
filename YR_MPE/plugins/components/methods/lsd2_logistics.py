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