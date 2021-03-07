# adapted from pyras

"""
"""
import os
import time

import win32con
import win32gui
import win32process


class Runtime(object):
    """ """
    def __init__(self, parent):
        self.window = None
        self.parent = parent
        self.parent_pid = None
        self.parent_window = None
        self.get_pid()

    def close(self):
        """ """
        kill_process(self.parent_pid)

    def get_pid(self):
        """ """
        self.parent.ShowRas()
        window_text = 'HEC-RAS '

        def enumHandler(hwnd, lParam):
            if window_text in win32gui.GetWindowText(hwnd):
                self.parent_window = hwnd
                return None

        win32gui.EnumWindows(enumHandler, None)
        _, pid = win32process.GetWindowThreadProcessId(self.parent_window)
        win32gui.ShowWindow(self.parent_window, win32con.SW_HIDE)
        self.parent_pid = pid

    # %% Handle GUI waiting for routines that do not stop runtime
    def pause_bc(self):
        """ """
        self._pause(window_text='Bridge Culvert Data')

    def pause_geo(self):
        """ """
        self._pause(window_text='Geometric Data')

    def pause_iw(self):
        """ """
        self._pause(window_text='Inline Structure Data')

    def pause_lw(self):
        """ """
        self._pause(window_text='Lateral Structure Editor')

    def pause_multiple(self):
        """ """
        self._pause(window_text='Run Multiple Plans')

    def pause_plan(self):
        """ """
        self._pause(window_text='Steady Flow Analysis')

    def pause_quasi(self):
        """ """
        self._pause(window_text='Quasi Unsteady Flow Editor')

    def pause_sediment(self):
        """ """
        self._pause(window_text='Sediment Data')

    def pause_steady(self):
        """ """
        self._pause(window_text='Steady Flow Data')

    def pause_unsteady(self):
        """ """
        self._pause(window_text='Unsteady Flow Data')

    def pause_quality(self):
        """ """
        self._pause(window_text='Water Quality Data')

    def pause_xs(self):
        """ """
        self._pause(window_text='Cross Section Data')

    def pause(self, time_seconds):
        """ """
        time.sleep(time_seconds)

    def _pause(self, window_text=None):
        """ """
        def enumHandler(hwnd, lParam):
            if window_text in win32gui.GetWindowText(hwnd):
                self.window = hwnd
        win32gui.EnumWindows(enumHandler, None)

        pause_check = True
        while pause_check:
            time.sleep(0.5)  # Prevent fan noise from CPU "over use"
            if not win32gui.IsWindowVisible(self.window):
                pause_check = False
                self.window = None


def kill_process(pid):
    """ """
    try:
        killed = os.system('TASKKILL /PID {} /F >nul'.format(pid))
    except Exception:
        killed = 0
    return killed
