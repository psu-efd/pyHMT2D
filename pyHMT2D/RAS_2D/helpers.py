# adapted from pyras

"""
HEC-RAS Controller
==================
"""

import os

import win32api
import win32con


def get_supported_versions():
    """ """
    return ['RAS600', 'RAS500']  # Order gives the priority


def kill_all():
    """ """
    import os
    import subprocess

    ras_process_string = 'ras.exe'
    proc = subprocess.Popen('TASKLIST /FO "CSV"', stdout=subprocess.PIPE)
    tasklist = proc.stdout.read().split('\n')
    tasks = []
    pids = []
    for line in tasklist:
        l = line.lower()
        if ras_process_string in l:
            items = l.split(',')
            tasks.append(items)
            pids.append(int(eval(items[1])))

    for pid in pids:
        try:
            os.system('TASKKILL /PID {0} /F >nul'.format(pid))
        except Exception as e:
            print(e)


def get_available_versions():
    """ """
    ver = {'HEC-RAS\\6.0 Beta Update 1\\Ras.exe': 'RAS600',
           'HEC-RAS\\5.0.7\\Ras.exe': 'RAS500'}

    ldic = _get_registered_typelibs()

    # Check if files actually exist (another sanity check)
    available_versions = []

    for dic in ldic:
        fname = dic['filename']
        if os.path.isfile(fname):
            for k in ver:
                if k in fname:
                    available_versions.append(ver[k])
    return available_versions


def _get_typelib_info(keyid, version):
    """
    adapted from pywin32

    # Copyright (c) 1996-2008, Greg Stein and Mark Hammond.
    """
    collected = []
    help_path = ""
    key = win32api.RegOpenKey(win32con.HKEY_CLASSES_ROOT,
                              "TypeLib\\%s\\%s" % (keyid, version))
    try:
        num = 0
        while 1:
            try:
                sub_key = win32api.RegEnumKey(key, num)
            except win32api.error:
                break
            h_sub_key = win32api.RegOpenKey(key, sub_key)
            try:
                value, typ = win32api.RegQueryValueEx(h_sub_key, None)
                if typ == win32con.REG_EXPAND_SZ:
                    value = win32api.ExpandEnvironmentStrings(value)
            except win32api.error:
                value = ""
            if sub_key == "HELPDIR":
                help_path = value
            elif sub_key == "Flags":
                flags = value
            else:
                try:
                    lcid = int(sub_key)
                    lcidkey = win32api.RegOpenKey(key, sub_key)
                    # Enumerate the platforms
                    lcidnum = 0
                    while 1:
                        try:
                            platform = win32api.RegEnumKey(lcidkey, lcidnum)
                        except win32api.error:
                            break
                        try:
                            hplatform = win32api.RegOpenKey(lcidkey, platform)
                            fname, typ = win32api.RegQueryValueEx(hplatform, None)
                            if typ == win32con.REG_EXPAND_SZ:
                                fname = win32api.ExpandEnvironmentStrings(fname)
                        except win32api.error:
                            fname = ""
                        collected.append((lcid, platform, fname))
                        lcidnum = lcidnum + 1
                    win32api.RegCloseKey(lcidkey)
                except ValueError:
                    pass
            num = num + 1
    finally:
        win32api.RegCloseKey(key)

    return fname, lcid


def _get_registered_typelibs(match='HEC River Analysis System'):
    """
    adapted from pywin32

    # Copyright (c) 1996-2008, Greg Stein and Mark Hammond.
    """
    # Explicit lookup in the registry.
    result = []
    key = win32api.RegOpenKey(win32con.HKEY_CLASSES_ROOT, "TypeLib")
    try:
        num = 0
        while 1:
            try:
                key_name = win32api.RegEnumKey(key, num)
            except win32api.error:
                break
            # Enumerate all version info
            sub_key = win32api.RegOpenKey(key, key_name)
            name = None
            try:
                sub_num = 0
                best_version = 0.0
                while 1:
                    try:
                        version_str = win32api.RegEnumKey(sub_key, sub_num)
                    except win32api.error:
                        break
                    try:
                        version_flt = float(version_str)
                    except ValueError:
                        version_flt = 0  # ????
                    if version_flt > best_version:
                        best_version = version_flt
                        name = win32api.RegQueryValue(sub_key, version_str)
                    sub_num = sub_num + 1
            finally:
                win32api.RegCloseKey(sub_key)
            if name is not None and match in name:
                fname, lcid = _get_typelib_info(key_name, version_str)

                # Split version
                major, minor = version_str.split('.')

                result.append({'name': name,
                               'filename': fname,
                               'iid': key_name,
                               'lcid': lcid,
                               'major': int(major),
                               'minor': int(minor)})
            num = num + 1
    finally:
        win32api.RegCloseKey(key)
    result = sorted(result)
    return result


# %%

def RAS_2D_Initial_Cleaning():
    """ Initial cleaning of HEC-RAS 2D environment

    Returns
    -------

    """

    kill_all()

    __available_versions__ = get_available_versions()

    if len(__available_versions__) > 0:
        for ras_version in get_supported_versions():
            if ras_version in __available_versions__:
                os.environ['RAS_CONTROLLER_VERSION'] = ras_version
                break

        #from .hecrascontroller import HECRASController

        # Cleaning the namespace
        globals().pop('hecrascontroller')
        globals().pop('hecrasgeometry')
        globals().pop('win32api')
        globals().pop('win32con')
        globals().pop('runtime')
        globals().pop('os')
    else:
        error = '"HEC River Analysis System" type library not found. ' \
                'Please install HEC-RAS'
        raise Exception(error)
