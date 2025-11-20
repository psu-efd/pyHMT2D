# -*- coding: utf-8 -*-
""" A set of classes to handle SRH-2D data

SRH-2D data include information in SRHHydro, SRHGeom, SRHMat, and simulation results.

Author: Xiaofeng Liu, PhD, PE
Penn State University
"""

import sys
import os
import copy
import numpy as np
import glob
import h5py
from os import path
import shlex
import vtk
from vtk.util import numpy_support as VN
import re
import meshio


from pyHMT2D.Hydraulic_Models_Data import HydraulicData

from pyHMT2D.__common__ import *
from pyHMT2D.Misc.vtk_utilities import vtkCellTypeMap
from pyHMT2D.Misc import vtkHandler

class SRH_2D_SIF:
    """Class to parse and manage SRH-2D SIF (Simulation Input File) data"""
    
    VALID_MODULES = ['RIVER', 'WATERSHED', 'COAST']
    VALID_SOLVERS = ['FLOW', 'MOBILE', 'DIFF', 'SED_DIFF']
    VALID_TIME_TYPES = ['STEADY', 'UNSTEADY']
    VALID_MESH_UNITS = ['FOOT', 'METER', 'INCH', 'MM', 'MILE', 'KM', 'GSCALE']
    VALID_INIT_CONDITIONS = ['DRY', 'RST', 'AUTO', 'ZONAL', 'VARY_WSE', 'VARY_WD']
    VALID_MANNING_OPTIONS = ['SPATIAL', 'SPATIAL VEG', 'SPATIAL GRAIN']
    VALID_OUTPUT_FORMATS = ['SRHC', 'TEC', 'SRHN', 'XMDF', 'XMDFC', 'VTK']
    VALID_OUTPUT_UNITS = ['SI', 'EN']

    def __init__(self, srhsif_filename):
        """Initialize SIF parser with srhsif_filename"""
        self.srhsif_filename = srhsif_filename

        self.srhsif_content = {
            'ManningN': {},
            'BC': {},
            'IQParams': {},
            'EWSParamsC': {},
            'MONITORING': {},
            'wall_roughness': {},
            'pressurized_zones': {},
            'flow_obstructions': {},
            'output_max_dat': False,
            'intermediate_output': {}
        }
        self.comments = {}

        try:
            self._parse_file()
        except Exception as e:
            raise ValueError(f"Error parsing SIF file: {str(e)}")

    def _parse_file(self):
        """Parse the SIF file and store data in self.data dictionary"""
        try:
            with open(self.srhsif_filename, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"SIF file not found: {self.srhsif_filename}")

        #get the case name from the filename: assume the filename is like "Muncie_SIF.dat"
        case_name = self.srhsif_filename.split("_SIF.dat")[0]  #strip off "_SIF.dat" 
        self.srhsif_content["Case"] = case_name        

        i = 0
        while i < len(lines):

            line = lines[i].strip()
            if line.startswith('//'):
                comment = line[2:].strip()
                self.comments[i] = comment
                i += 1
                if i >= len(lines):
                    break
                value = lines[i].strip()
                
                try:
                    if "Simulation Description" in comment:
                        self.srhsif_content['simulation_description'] = value
                    elif "Module Selected" in comment:
                        self._validate_module(value.upper())
                        self.srhsif_content['module'] = value
                    elif "Solver Selection" in comment:
                        self._validate_solver(value.upper())
                        self.srhsif_content['solver'] = value
                    elif "Monitor-Point-Info" in comment:
                        self.srhsif_content['monitor_point_npoints'] = int(value)    #number of monitor points
                        if self.srhsif_content['monitor_point_npoints'] > 0:
                            self._parse_monitoring_points(lines, i)
                    elif "Steady-or-Unsteady" in comment:
                        self._validate_time_type(value.upper())
                        self.srhsif_content['time_type'] = value
                    elif "Time Parameters" in comment:
                        self.srhsif_content['time_params'] = self._parse_time_parameters(value)
                    elif "Turbulence-Model" in comment:
                        self.srhsif_content['turbulence_model'] = value
                    elif "A_TURB for the PARA Model" in comment:
                        self.srhsif_content['a_turb'] = float(value)
                    elif "Mesh-Unit" in comment:
                        self.srhsif_content['mesh_unit'] = value
                    elif "Mesh FILE_NAME and FORMAT" in comment:
                        self.srhsif_content['mesh_file_name'] = value.split()[0]
                        self.srhsif_content['mesh_format'] = value.split()[1]
                    elif "Initial Flow Condition Setup Option" in comment:
                        self._validate_init_condition(value.upper())
                        self.srhsif_content['initial_flow_condition'] = value

                        #if the initial flow condition is ZONAL, then the next lines are for the number of initial condition zones and the zone data
                        if value.upper() == 'ZONAL':
                            #dict for zonal IC data
                            res_IC_zonal = {}
                            
                            i += 1  #comment line
                            #check the current comment line contains "Constant Setup for Initial Condition"
                            if "Constant Setup for Initial Condition" in lines[i]:
                                #read the next line
                                i += 1
                                value = lines[i].strip()
                                self.srhsif_content['n_initial_condition_zones'] = int(value)

                                #read the next comment line and make sure it contains "Constant-Value Initial Condition"
                                i += 1
                                if "Constant-Value Initial Condition" in lines[i]:
                                    #read the next n_initial_condition_zones lines (one line per zone)
                                    for j in range(self.srhsif_content['n_initial_condition_zones']):
                                        i += 1
                                        values = lines[i].strip().split()   #values is a list of strings; the first three are floats, the rest are strings 
                                        values_float = [float(x) for x in values[:3]]
                                        values_str = values[3:]
                                        res_IC_zonal[j] = values_float + values_str
                                else:
                                    raise ValueError("The next line after the initial flow condition is not 'Constant-Value Initial Condition'.")

                            else:
                                raise ValueError("The next line after the initial flow condition is not 'Constant Setup for Initial Condition'.")
                            
                            self.srhsif_content['IC_zonal'] = res_IC_zonal
                            

                    elif "Manning Coefficient" in comment:
                        i = self._parse_manning_coefficients(lines, i, value)
                    elif "Any-Special-Modeling-Options" in comment:
                        self.srhsif_content['special_modeling_options'] = int(value)
                    elif "Boundary Type" in comment:
                        i = self._parse_boundary_conditions(lines, i)
                    elif "Wall-Roughess-Height-Specification" in comment:
                        i = self._parse_wall_roughness(lines, i)
                    elif "Pressurized Zone exists?" in comment:
                        i = self._parse_pressurized_zones(lines, i)
                    elif "Any In-Stream Flow Obstructions?" in comment:
                        i = self._parse_flow_obstructions(lines, i)       
                    elif "Results-Output-Format" in comment:
                        self._parse_output_format(value)
                    elif "Output File _MAX.dat" in comment:
                        self._parse_max_output(value)
                    elif "Intermediate Result Output Control" in comment:
                        i = self._parse_intermediate_output(lines, i)
                except ValueError as e:
                    raise ValueError(f"Error parsing line {i+1}: {str(e)}")
            i += 1

    def _parse_boundary_conditions(self, lines, start_index):
        """Parse all boundary conditions from the file
        
        Args:
            lines: List of all file lines
            start_index: Current line index where first boundary is found
        
        Returns:
            next_index: Index of the next line after all boundaries
        """

        res_BC = {}
        res_IQParams = {}
        res_EWSParamsC = {}

        i = start_index - 1  # Start from the line before the first boundary for the while loop below        


        index_BC = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Check if we've hit the next section (indicated by '//' without 'Boundary')
            if line.startswith('//') and 'Boundary' not in line:
                break
                
            # Parse boundary type
            if line.startswith('// Boundary Type'):
                index_BC += 1
                i += 1
                boundary_type = lines[i].strip().lower()

                res_BC[index_BC] = boundary_type

                if boundary_type == 'monitoring' or boundary_type == 'monitor' or boundary_type == 'symmetry' or boundary_type == 'wall':
                    #do nothing (there is no boundary values for monitoring, symmetry, or symmetry)
                    pass
                else:
                    # Get boundary values from next two lines
                    i += 1  # Skip to values comment line
                    i += 1  # Go to actual values line
                    boundary_values = lines[i].strip().split()

                    if boundary_type == 'inlet-q':
                        res_IQParams[index_BC] = boundary_values
                    elif boundary_type == 'exit-h':
                        res_EWSParamsC[index_BC] = boundary_values                    
                    else:
                        raise ValueError(f"Boundary type: {boundary_type} is not supported yet.")
            
            i += 1

        self.srhsif_content['BC'] = res_BC
        self.srhsif_content['IQParams'] = res_IQParams
        self.srhsif_content['EWSParamsC'] = res_EWSParamsC

        #print("res_BC = ", res_BC)
        #print("res_IQParams = ", res_IQParams)
        #print("res_EWSParamsC = ", res_EWSParamsC)
        
        return i  # Return the index of the next section


    def _parse_wall_roughness(self, lines, i):
        """Parse wall roughness section"""
        value = lines[i].strip()
        if value and value not in ['0', '']:
            self.srhsif_content['wall_roughness'] = value.split()
        return i

    def _parse_pressurized_zones(self, lines, i):
        """Parse pressurized zones section"""
        value = lines[i].strip()
        if value and value not in ['0', '']:
            self.srhsif_content['pressurized_zones'] = value.split()
        return i

    def _parse_flow_obstructions(self, lines, i):
        """Parse flow obstructions section"""
        value = lines[i].strip()
        if value and value not in ['0', '']:
            self.srhsif_content['flow_obstructions'] = value.split()
        return i

    def _parse_output_format(self, value):
        """Parse output format and units"""
        parts = value.split()
        format_parts = parts[0].split('/')
        output_format = format_parts[0].upper()
        if output_format not in self.VALID_OUTPUT_FORMATS:
            raise ValueError(f"Invalid output format: {output_format}")
        
        output_unit = parts[1].upper() if len(parts) > 1 else 'EN'
        if output_unit not in self.VALID_OUTPUT_UNITS:
            raise ValueError(f"Invalid output unit: {output_unit}")
        
        self.srhsif_content['OutputFormat'] = output_format
        self.srhsif_content['OutputUnit'] = output_unit
        
        # Check for optional STL FACE
        if len(parts) > 2 and parts[2:] == ['STL', 'FACE']:
            self.srhsif_content['output_stl_face'] = True

    def _parse_max_output(self, value):
        """Parse _MAX.dat output option"""
        self.srhsif_content['output_max_dat'] = bool(value.strip())

    def _parse_intermediate_output(self, lines, i):
        """Parse intermediate output control"""
        value = lines[i].strip()
        if value.lower() != 'empty':
            try:
                # Try parsing as a single interval
                interval = float(value)
                self.srhsif_content['intermediate_output'] = {'type': 'interval', 'value': interval}
            except ValueError:
                # Parse as list of times
                times = [float(t) for t in value.split()]
                self.srhsif_content['intermediate_output'] = {'type': 'list', 'value': times}
        return i

    def save_as(self, srhsif_filename=None):
        """Save the current data in self.srhsif_content back to a SIF file

        So to modify a SIF file, one needs to change the data in self.srhsif_content, and then call this function to save the data back to the SIF file.
        
        Parameters
        ----------
        srhsif_filename : str, optional
            The filename to save the SIF file to. If not provided, the filename will be the same as the original filename.

        Returns
        -------
        """
        if srhsif_filename is None:
            srhsif_filename = self.srhsif_filename

        try:
            with open(srhsif_filename, 'w') as f:
                # Simulation Description
                f.write("// Simulation Description (not used by SRH):\n")
                f.write(f"{self.srhsif_content.get('simulation_description', '')}\n")
                
                # Module
                f.write("// Module Selected (RIVER WATERSHED COAST)\n")
                f.write(f"{self.srhsif_content.get('module', '')}\n")
                
                # Solver
                f.write("// Solver Selection (FLOW MOBile DIFF SED_DIFF ...)\n")
                f.write(f"{self.srhsif_content.get('solver', '')}\n")
                
                # Monitor Points
                f.write("// Monitor-Point-Npoints: NPOINT\n")
                f.write(f"{self.srhsif_content.get('monitor_point_npoints', '0')}\n")

                # Monitoring Points
                if self.srhsif_content.get('monitor_point_npoints') > 0:
                    f.write("// Monitor Point Coordinates: x1 y1 x2 y2 ...\n")

                    #loop over the monitor points and write them to the file
                    for pointI in range(self.srhsif_content['monitor_point_npoints']):
                        point_x = self.srhsif_content['monitor_points'][pointI*2]
                        point_y = self.srhsif_content['monitor_points'][pointI*2+1]
                        f.write(f"{point_x} {point_y}\n")                 
                
                # Time Type
                f.write("// Steady-or-Unsteady (STEADY/UNS)\n")
                f.write(f"{self.srhsif_content.get('time_type', '')}\n")
                
                # Time Parameters
                f.write("// Time Parameters: T_Start(hr) T_end(hr) Dt(s) [Dt_max(s) CFL]\n")
                f.write(" ".join(map(str, self.srhsif_content['time_params'])) + "\n")
                    
                # Turbulence Model
                f.write("// Turbulence-Model\n")
                f.write(f"{self.srhsif_content.get('turbulence_model', '')}\n")
                
                # A_TURB Parameter
                f.write("// A_TURB for the PARA Model (0.05 to 1.0)\n")
                f.write(f"{self.srhsif_content.get('a_turb', '')}\n")
                
                # Mesh Unit
                f.write("// Mesh-Unit (FOOT METER INCH MM MILE KM GSCALE)\n")
                f.write(f"{self.srhsif_content.get('mesh_unit', '')}\n")
                
                # Mesh File
                f.write("// Mesh FILE_NAME and FORMAT(SMS...)\n")
                f.write(f"{self.srhsif_content.get('mesh_file_name', '')} {str(self.srhsif_content.get('mesh_format', '')).strip()}\n")
                
                # Initial Flow Condition
                f.write("// Initial Flow Condition Setup Option (DRY RST AUTO ZONAL Vary_WSE/Vary_WD)\n")
                f.write(f"{self.srhsif_content.get('initial_flow_condition', '')}\n")

                # If the initial flow condition is ZONAL, then the next lines are for the number of initial condition zones and the zone data
                if self.srhsif_content.get('initial_flow_condition', '').upper() == 'ZONAL':
                    f.write("// Constant Setup for Initial Condition: n_zone [MESH or 2DM_filename]\n")
                    f.write(f"{self.srhsif_content.get('n_initial_condition_zones', '0')}\n")
                    f.write("// Constant-Value Initial Condition for Mesh Zone: U V WSE [TK] [ED] [T]\n")
                    for j in range(self.srhsif_content['n_initial_condition_zones']):
                        #f.write(f"{self.srhsif_content['IC_zonal'][j]}\n")
                        f.write(f"{' '.join(map(str, self.srhsif_content['IC_zonal'][j]))}\n")
                
                # Manning Coefficients
                f.write("// Manning Coefficient n Input Options: SPATIAL or SPATIAL VEG GRAIN for 2D Model; SPATIAL for 3D model\n")
                f.write(f"{self.srhsif_content.get('manning_option', '')}\n")
                f.write("// Number of Material Types in 2D Mesh File\n")
                f.write(f"{str(self.srhsif_content.get('n_manning_zones', '')).strip()}\n")
                f.write("// Manning Coefficient in each mesh zone: a real value or a WD~n file name or Landuse\n")
                if 'ManningN' in self.srhsif_content:
                    # Sort the ManningN values by their keys before writing
                    sorted_values = sorted(self.srhsif_content['ManningN'].items(), key=lambda x: x[0])
                    for key, value in sorted_values:
                        if key != 0:    #Don't need the default Manning's n value
                            f.write(f"{value}\n")                        
                else:
                    #something is wrong
                    raise ValueError("ManningN is not found in the SIF file.")
                
                # Special Modeling Options
                f.write("// Any-Special-Modeling-Options? (0/1=no/yes)\n")
                f.write(f"{self.srhsif_content.get('special_modeling_options', '0')}\n")
                
                # Boundary Conditions
                bcDict = self.srhsif_content['BC']
                IQParams = self.srhsif_content['IQParams']
                EWSParamsC = self.srhsif_content['EWSParamsC']

                # loop over all the boundaries
                for index_BC, boundary_type in bcDict.items():
                    f.write("// Boundary Type (INLET-Q EXIT-H etc)\n")
                    f.write(f"{boundary_type}\n")                    
                    f.write("// Boundary Values (Q W QS TEM H_rough etc)\n")

                    if boundary_type == 'inlet-q':
                        f.write(" ".join(str(x) for x in IQParams[index_BC]) + "\n")
                    elif boundary_type == 'exit-h':
                        f.write(" ".join(str(x) for x in EWSParamsC[index_BC]) + "\n")
                    elif boundary_type == 'monitoring' or boundary_type == 'monitor' or boundary_type == 'symmetry':
                        #do nothing  
                        pass
                    else:
                        raise ValueError(f"Boundary type: {boundary_type} is not supported yet.")
              
                
                # Wall Roughness
                f.write("// Wall-Roughess-Height-Specification (empty-line=DONE)\n")
                if self.srhsif_content.get('wall_roughness'):
                    f.write(" ".join(map(str, self.srhsif_content['wall_roughness'])) + "\n")
                else:
                    f.write(" \n")
                
                # Pressurized Zones
                f.write("// Pressurized Zone exists? (empty-line or 0 == NO)\n")
                if self.srhsif_content.get('pressurized_zones'):
                    f.write(" ".join(map(str, self.srhsif_content['pressurized_zones'])) + "\n")
                else:
                    f.write(" \n")
                
                # Flow Obstructions
                f.write("// Any In-Stream Flow Obstructions? (empty-line or 0 = NO)\n")
                if self.srhsif_content.get('flow_obstructions'):
                    f.write(" ".join(map(str, self.srhsif_content['flow_obstructions'])) + "\n")
                else:
                    f.write(" \n")
                
                # Output Format
                f.write("// Results-Output-Format-and-Unit(SRHC/TEC/SRHN/XMDF/XMDFC/VTK;SI/EN) + Optional STL FACE\n")
                output_str = f"{self.srhsif_content.get('OutputFormat', '')} {self.srhsif_content.get('OutputUnit', '')}"
                if self.srhsif_content.get('output_stl_face'):
                    output_str += " STL FACE"
                f.write(output_str + "\n")
                
                # MAX.dat Output
                f.write("// Output File _MAX.dat is requested? (empty means NO)\n")
                f.write("1\n" if self.srhsif_content.get('output_max_dat') else " \n")
                
                # Intermediate Output
                f.write("// Intermediate Result Output Control: INTERVAL(hour) OR List of T1 T2 ...  EMPTY means the end\n")
                if 'intermediate_output' in self.srhsif_content:
                    output_data = self.srhsif_content['intermediate_output']
                    if output_data['type'] == 'interval':
                        f.write(f"{output_data['value']}\n")
                    else:
                        f.write(" ".join(map(str, output_data['value'])) + "\n")
                else:
                    f.write("EMPTY\n")

        except Exception as e:
            raise IOError(f"Error saving SIF file: {str(e)}")

    # Additional getter methods
    def get_simulation_start_end_time(self):
        """Return simulation start and end time"""
        return self.srhsif_content['time_params'][0], self.srhsif_content['time_params'][1]

    def get_simulation_time_step_size(self):
        """Return simulation time step size"""
        return self.srhsif_content['time_params'][2]

    def get_ManningN_dict(self):
        """Return Manning's n values"""
        return self.srhsif_content.get('ManningN', [])

    def get_wall_roughness(self):
        """Return wall roughness specifications"""
        return self.data.get('wall_roughness', [])


    def get_pressurized_zones(self):
        """Return pressurized zones specifications"""
        return self.data.get('pressurized_zones', [])

    def get_flow_obstructions(self):
        """Return flow obstructions specifications"""
        return self.data.get('flow_obstructions', [])

    def get_output_format(self):
        """Return output format and unit"""
        return {
            'format': self.data.get('OutputFormat'),
            'unit': self.data.get('OutputUnit'),
            'stl_face': self.data.get('output_stl_face', False)
        }

    def get_intermediate_output(self):
        """Return intermediate output control settings"""
        return self.data.get('intermediate_output', None)

    # Additional setter methods
    def set_wall_roughness(self, values):
        """Set wall roughness specifications"""
        self.data['wall_roughness'] = values

    def set_output_format(self, format_type, unit='EN', stl_face=False):
        """Set output format and unit"""
        if format_type not in self.VALID_OUTPUT_FORMATS:
            raise ValueError(f"Invalid output format: {format_type}")
        if unit not in self.VALID_OUTPUT_UNITS:
            raise ValueError(f"Invalid output unit: {unit}")
        
        self.data['output_format'] = format_type
        self.data['output_unit'] = unit
        self.data['output_stl_face'] = stl_face

    def set_intermediate_output(self, value, output_type='interval'):
        """Set intermediate output control
        
        Args:
            value: float or list of floats
            output_type: 'interval' or 'list'
        """
        if output_type not in ['interval', 'list']:
            raise ValueError("output_type must be 'interval' or 'list'")
            
        self.data['intermediate_output'] = {
            'type': output_type,
            'value': value
        }

    def _validate_module(self, value):
        """Validate module selection"""
        if value.upper() not in self.VALID_MODULES:
            raise ValueError(f"Invalid module: {value}. Must be one of {self.VALID_MODULES}")

    def _validate_solver(self, value):
        """Validate solver selection"""
        if value.upper() not in self.VALID_SOLVERS:
            raise ValueError(f"Invalid solver: {value}. Must be one of {self.VALID_SOLVERS}")

    def _validate_time_type(self, value):
        """Validate time type selection"""
        if value.upper() not in self.VALID_TIME_TYPES:
            raise ValueError(f"Invalid time type: {value}. Must be one of {self.VALID_TIME_TYPES}")

    def _validate_mesh_unit(self, value):
        """Validate mesh unit selection"""
        if value.upper() not in self.VALID_MESH_UNITS:
            raise ValueError(f"Invalid mesh unit: {value}. Must be one of {self.VALID_MESH_UNITS}")

    def _validate_init_condition(self, value):
        """Validate initial condition selection"""
        if value.upper() not in self.VALID_INIT_CONDITIONS:
            raise ValueError(f"Invalid initial condition: {value}. Must be one of {self.VALID_INIT_CONDITIONS}")

    def _validate_manning_option(self, value):
        """Validate Manning option selection"""
        if value.upper() not in self.VALID_MANNING_OPTIONS:
            raise ValueError(f"Invalid Manning option: {value}. Must be one of {self.VALID_MANNING_OPTIONS}")

    def _parse_time_parameters(self, value):
        """Parse time parameters line
        
        Format: T_Start(hr) T_end(hr) Dt(s) [Dt_max(s) CFL]
        Last two parameters are optional
        
        Args:
            value: String containing space-separated time parameters
            
        Returns:
            List of float values for time parameters
            
        Raises:
            ValueError: If parameters are invalid or missing required values
        """
        try:
            params = [float(x) for x in value.split()]
            if len(params) < 3:
                raise ValueError("Time parameters must include at least t_start, t_end, and dt")
            return params
        except ValueError as e:
            raise ValueError(f"Invalid time parameters format: {str(e)}")

    def _parse_monitoring_points(self, lines, i):
        """Parse monitoring points section
        
        Args:
            lines: List of all file lines
            i: Current line index (where NPOINT is specified)
            
        Returns:
            next_index: Index of the next line after monitoring points section
        """
        try:
            # Get number of monitoring points            
            n_points = int(lines[i].strip())
            
            # Skip comment line
            i += 1
            if not lines[i].strip().startswith('//'):
                raise ValueError("Expected comment line for monitoring point coordinates")
                
            # Get coordinates (x, y) and one point per line
            coords = []   #list of x1, y1, x2, y2, ...
            for j in range(n_points):
                i += 1
                #print("lines[i].strip().split(): ", lines[i].strip().split())
                #print("type(lines[i].strip().split()): ", type(lines[i].strip().split()))
                for item in lines[i].strip().split():
                    coords.append(float(item))

            #print(f"coords: {coords}")
            
            # Convert coordinates to float and pair them
            if len(coords) != 2 * n_points:
                raise ValueError(f"Expected {2*n_points} coordinates for {n_points} points: length of coords is {len(coords)}")
                
            #coords = [float(x) for x in coords]
                            
            # Store in data dictionary            
            self.srhsif_content['monitor_points'] = coords
            
            return i
            
        except ValueError as e:
            raise ValueError(f"Error parsing monitoring points: {str(e)}")
    
    def _parse_manning_coefficients(self, lines, i, value):
        """Parse Manning coefficients section
        
        Args:
            lines: List of all file lines
            i: Current line index where "VARY" is found
            value: String containing "VARY"
            
        Returns:
            next_index: Index of the next line after Manning coefficients
            
        Raises:
            ValueError: If Manning coefficient values are invalid
        """
        try:
            # First line should be the option, e.g., "VARY"
            manning_option = value.upper()
            self.srhsif_content['manning_option'] = manning_option
            
            # Next line should be the comment about number of materials
            i += 1
            if i >= len(lines) or not lines[i].strip().startswith('//'):
                raise ValueError("Expected comment line for number of materials")
                
            # Next line contains the actual number
            i += 1
            if i >= len(lines):
                raise ValueError("Unexpected end of file while reading number of materials")
            
            try:
                n_materials = int(lines[i].strip())
            except ValueError:
                raise ValueError(f"Invalid number of materials: {lines[i].strip()}")
                
            self.srhsif_content['n_manning_zones'] = n_materials

            # Read the comment line for the Manning values
            i += 1
            if i >= len(lines) or not lines[i].strip().startswith('//'):
                raise ValueError("Expected comment line for Manning values")
            
            # Read the Manning values
            res_ManningsN = {}  # dict for ManningsN (there cuould be multiple entries)

            #add a default value for the ManningN (not sure why SIF file does not have this)
            res_ManningsN[0] = 0.03

            for material_id in range(n_materials):
                i += 1

                if i >= len(lines):
                    raise ValueError("Unexpected end of file while reading Manning coefficients")

                value = lines[i].strip()
                # Handle both numeric values and file names
                try:
                    res_ManningsN[material_id+1] = float(value)
                except ValueError:
                    # If not a number, store as string (could be a file name)
                    res_ManningsN[material_id+1] = value
            
            self.srhsif_content['ManningN'] = res_ManningsN
            return i

        except ValueError as e:
            raise ValueError(f"Invalid Manning coefficient values: {str(e)}")
        
    def modify_one_manning_n_in_srhsif_content(self, material_id, new_manning_n):
        """Modify one Manning's n value for a specific material ID in self.srhsif_content"""

        #check material_id is an integer and is within the range of the Manning's n list
        if not isinstance(material_id, int) or material_id < 0 or material_id >= len(self.srhsif_content['ManningN']):
            raise ValueError(f"Invalid material ID: {material_id}. Must be an integer within the range of the Manning's n list.")

        #check new_manning_n is a float
        if not isinstance(new_manning_n, float):
            raise ValueError(f"Invalid Manning's n value: {new_manning_n}. Must be a float.")

        #modify the Manning's n value
        self.srhsif_content['ManningN'][material_id] = new_manning_n

        if gVerbose:
            print(f"Manning's n value for material ID {material_id} has been modified to {new_manning_n}")


    def modify_ManningsN(self, materialIDs, newManningsNValues, ManningN_MaterialNames):
        """Modify materialID's Manning's n values to new values for a given list of material IDs

        Parameters
        ----------
        materialIDs : list
            material ID list
        newManningsNValues : list
            new Manning's n value list
        ManningN_MaterialNames : list
            name of materials list

        Returns
        -------

        """

        if gVerbose: print("Modify Manning's n values ...")

        if not isinstance(materialIDs[0], int):
            raise Exception("Material ID has to be an integer. The type of materialID passed in is ", type(materialIDs[0]))

        if not isinstance(newManningsNValues[0], float):
            raise Exception("Manning's n has to be a float. The type of newManningsNValue passed in is ", type(newManningsNValues[0]))

        #get the Manning's n dictionary
        nDict = self.get_ManningN_dict()

        #loop over all the material IDs and modify the Manning's n values
        for i in range(len(materialIDs)):
            materialID = materialIDs[i]

            if materialID in nDict:
                if gVerbose: print("    Old Manning's n value =", nDict[materialID], "for material ID = ", materialID)
                nDict[materialID] = newManningsNValues[i]
                if gVerbose: print("    New Manning's n value =", nDict[materialID], "for material ID = ", materialID)

                #also modify the srhsif_content
                self.modify_one_manning_n_in_srhsif_content(materialID, newManningsNValues[i])
            else:
                print("The specified materialID", materialID, "is not in the Manning's n list. Please check.")

    def modify_InletQ(self, bcIDs, newInletQValues):
        """Modify the inlet flow rate for the specified boundary IDs

        Parameters
        ----------
        bcIDs : list
            list of boundary IDs
        newInletQValues : list
            list of new inlet flow rate values

        Returns
        -------

        """

        if gVerbose: print("Modifying inlet flow rate for the specified boundary IDs ...")

        if not isinstance(bcIDs[0], int):
            raise Exception("Boundary ID has to be an integer. The type of bcID passed in is ", type(bcIDs[0]))

        if not isinstance(newInletQValues[0], float):
            raise Exception("Inlet flow rate has to be a float. The type of newInletQValues passed in is ", type(newInletQValues[0]))
        
        #get the inlet flow rate dictionary
        bcDict = self.srhsif_content['BC']
        IQParamsDict = self.srhsif_content['IQParams']

        #loop over all the boundary IDs and modify the inlet flow rate values
        for i in range(len(bcIDs)):
            bcID = bcIDs[i]

            if bcID in bcDict:
                if gVerbose: print("    Old inlet flow rate =", IQParamsDict[bcID][0], "for boundary ID = ", bcID)
                IQParamsDict[bcID][0] = newInletQValues[i]
                if gVerbose: print("    New inlet flow rate =", IQParamsDict[bcID][0], "for boundary ID = ", bcID)
            else:
                print("The specified boundary ID", bcID, "is not in the inlet flow rate list. Please check.")

    def modify_ExitH(self, bcIDs, newExitHValues):
        """Modify the exit water surface elevation for the specified boundary IDs

        Parameters
        ----------
        bcIDs : list
            list of boundary IDs
        newExitHValues : list
            list of new exit water surface elevation values

        Returns
        -------

        """

        if gVerbose: print("Modifying exit water surface elevation for the specified boundary IDs ...")

        if not isinstance(bcIDs[0], int):
            raise Exception("Boundary ID has to be an integer. The type of bcID passed in is ",
                            type(bcIDs[0]))

        if not isinstance(newExitHValues[0], float):
            raise Exception("Exit water surface elevation has to be a float. The type of newExitHValues passed in is ",
                            type(newExitHValues[0]))

        #get the exit water surface elevation dictionary
        bcDict = self.srhsif_content['BC']
        EWSParamsC_Dict = self.srhsif_content['EWSParamsC']

        #loop over all the boundary IDs and modify the exit water surface elevation values
        for i in range(len(bcIDs)):
            bcID = bcIDs[i]

            if bcID in bcDict:
                if gVerbose: print("    Old exit water surface elevation =", EWSParamsC_Dict[bcID][0], "for boundary ID = ", bcID)
                EWSParamsC_Dict[bcID][0] = newExitHValues[i]
                if gVerbose: print("    New exit water surface elevation =", EWSParamsC_Dict[bcID][0], "for boundary ID = ", bcID)
            else:
                print("The specified boundary ID", bcID, "is not in the exit water surface elevation list. Please check.")


class SRH_2D_SRHHydro:

    """A class to handle srhhydro file for SRH-2D

    Attributes
    ----------
    srhhydro_filename : str
        name for the srhhydro file
    srhhydro_content : dictionary
        a dictionary to hold all information in the srhhydro file


    Methods


    """

    def __init__(self, srhhydro_filename):
        """ Constructor for SRH_2D_SRHHydro

        The SRH_2D_SRHHydro file contains all information for a SRH-2D run.

        Parameters
        ----------
        srhhydro_filename: str
            The name of the SRHHydro file
        """

        self.srhhydro_filename = srhhydro_filename #srhhydro file name

        #dict to hold the content of the srhhydro file
        self.srhhydro_content = {}

        #parse the srhhydro file and build srhhydro_content
        self.parse_srhhydro_file()

    def parse_srhhydro_file(self):
        """Parse the SRHHydro file"""
        
        res_all = {}
        res_ManningsN = {}
        res_InitZoneParams = {}
        
        res_BC = {}
        res_MONITORING = {}
        res_PRESSURE = {}
        res_INLET_Q = {}
        res_EXIT_H = {}
        res_WALL = {}

        res_IQParams = {}
        res_ISupCrParams = {}
        res_EWSParamsC = {}
        res_EWSParamsRC = {}
        res_EQParams = {}
        res_NDParams = {}
        res_PressureNodestrings = {}
        res_PressureParams = {}
        res_PressOvertop = {}  # Add dictionary for pressure overtop parameters
        res_PressWeirParams = {}  # Add dictionary for pressure weir parameters
        
        try:
            with open(self.srhhydro_filename, 'r') as f:
                lines = f.readlines()
                
            for i, line in enumerate(lines):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                
                # Handle different line types
                if line.startswith('SRHHYDRO'):
                    res_all['Version'] = parts[1]
                elif line.startswith('Case'):
                    res_all['Case'] = parts[1].strip('"')
                elif line.startswith('Description'):
                    res_all['Description'] = ' '.join(parts[1:]).strip('"')
                elif line.startswith('RunType'):
                    res_all['RunType'] = parts[1]
                elif line.startswith('ModelTemp'):
                    res_all['ModelTemp'] = parts[1]
                elif line.startswith('UnsteadyOutput'):
                    res_all['UnsteadyOutput'] = parts[1]
                elif line.startswith('SimTime'):
                    res_all['SimTime'] = [float(x) for x in parts[1:]]
                elif line.startswith('TurbulenceModel'):
                    res_all['TurbulenceModel'] = parts[1]
                elif line.startswith('ParabolicTurbulence'):
                    res_all['ParabolicTurbulence'] = float(parts[1])
                elif line.startswith('InitCondOption'):
                    res_all['InitCondOption'] = parts[1]
                elif line.startswith('NumInitZones'):
                    res_all['NumInitZones'] = int(parts[1])
                elif line.startswith('InitZoneParams'):
                    zone_id = int(parts[1])
                    p1_value = int(parts[2])
                    p2_value = int(parts[3])
                    p3_value = float(parts[4])
                    unit = parts[5]
                    res_InitZoneParams[zone_id] = {'p1': p1_value, 'p2': p2_value, 'p3': p3_value, 'unit': unit}
                elif line.startswith('Grid'):
                    res_all['Grid'] = parts[1].strip('"')
                elif line.startswith('HydroMat'):
                    res_all['HydroMat'] = parts[1].strip('"')
                elif line.startswith('MonitorPtFile'):
                    res_all['MonitorPtFile'] = parts[1].strip('"')
                elif line.startswith('OutputFormat'):
                    res_all['OutputFormat'] = parts[1]
                    res_all['OutputUnit'] = parts[2]
                elif line.startswith('OutputOption'):
                    res_all['OutputOption'] = int(parts[1])
                elif line.startswith('OutputInterval'):
                    res_all['OutputInterval'] = float(parts[1])
                elif line.startswith('ManningsN'):
                    zone_id = int(parts[1])
                    n_value = float(parts[2])
                    res_ManningsN[zone_id] = n_value
                elif line.startswith('BC'):
                    bc_id = int(parts[1])
                    bc_type = parts[2]
                    res_BC[bc_id] = bc_type
                    
                    # Track monitoring and pressure BCs
                    if bc_type == 'MONITORING':
                        res_MONITORING[int(parts[1])] = parts[2]
                    elif bc_type == 'PRESSURE':
                        res_PRESSURE[int(parts[1])] = parts[2]
                    elif bc_type == 'INLET-Q':
                        res_INLET_Q[int(parts[1])] = parts[2]
                    elif bc_type == 'EXIT-H':
                        res_EXIT_H[int(parts[1])] = parts[2]
                    elif bc_type == 'WALL':
                        res_WALL[int(parts[1])] = parts[2]
                elif line.startswith('IQParams'):
                    bc_id = int(parts[1])
                    res_IQParams[bc_id] = {
                        'discharge': float(parts[2]),
                        'unit': parts[3],
                        'type': parts[4]
                    }
                elif line.startswith('EWSParamsC'):
                    bc_id = int(parts[1])
                    res_EWSParamsC[bc_id] = {
                        'wse': float(parts[2]),
                        'unit': parts[3],
                        'type': parts[4]
                    }
                elif line.startswith('EWSParamsRC'):
                    bc_id = int(parts[1])
                    res_EWSParamsRC[bc_id] = {
                        'curve_file': parts[2].strip('"'),
                        'unit': parts[3],
                        'type': parts[4]
                    }
                elif line.startswith('PressureNodestrings'):
                    group_id = int(parts[1])
                    nodestring_ids = [int(x) for x in parts[2:]]
                    res_PressureNodestrings[group_id] = nodestring_ids
                elif line.startswith('PressureParams'):
                    nodestring_id = int(parts[1])
                    res_PressureParams[nodestring_id] = {
                        'type': parts[2],
                        'var1': float(parts[3]),    #not sure the meaning of var1, var2, var3
                        'var2': float(parts[4]),
                        'var3': float(parts[5]),
                        'unit': parts[6]
                    }

                    #If the current line is PressureParams, the next line should be PressOvertop
                    next_line = lines[i + 1].strip()
                    if next_line.startswith('PressOvertop'):
                        nodestring_id_overtop = int(parts[1])

                        if nodestring_id_overtop != nodestring_id:
                            raise ValueError(f"Mismatched nodestring IDs between PressureParams and PressOvertop: {nodestring_id} != {nodestring_id_overtop}")

                        overtop_flag = int(parts[2])

                        res_PressureParams[nodestring_id]['overtop_flag'] = overtop_flag                        
                        
                        # If overtop_flag is 1, next line should be PressWeirParams2
                        if overtop_flag == 1:
                            next_line = lines[i + 1].strip()
                            if next_line.startswith('PressWeirParams2'):
                                weir_parts = next_line.split()
                                weir_id = int(weir_parts[1])
                                if weir_id != nodestring_id:
                                    raise ValueError(f"Mismatched nodestring IDs between PressOvertop and PressWeirParams2: {nodestring_id} != {weir_id}")
                                
                                res_PressWeirParams[nodestring_id] = {
                                    'crest_elevation': float(weir_parts[2]),
                                    'weir_coefficient': float(weir_parts[3]),
                                    'unit': weir_parts[4],
                                    'type': weir_parts[5]
                                }
                            else:
                                raise ValueError(f"PressWeirParams2 not found for nodestring ID {nodestring_id}")
                    else:                        
                        raise ValueError(f"PressOvertop not found for nodestring ID {nodestring_id}")
                
            # Store all results in srhhydro_content
            self.srhhydro_content = {
                **res_all,
                'ManningsN': res_ManningsN,
                'InitZoneParams': res_InitZoneParams,
                'BC': res_BC,
                'MONITORING': res_MONITORING,
                'PRESSURE': res_PRESSURE,
                'INLET-Q': res_INLET_Q,
                'EXIT-H': res_EXIT_H,
                'WALL': res_WALL,
                'IQParams': res_IQParams,
                'ISupCrParams': res_ISupCrParams,
                'EWSParamsC': res_EWSParamsC,
                'EWSParamsRC': res_EWSParamsRC,
                'EQParams': res_EQParams,
                'NDParams': res_NDParams,
                'PressOvertop': res_PressOvertop,
                'PressWeirParams': res_PressWeirParams,
                'PressureNodestrings': res_PressureNodestrings,
                'PressureParams': res_PressureParams
            }
            
        except Exception as e:
            print(f"Error parsing SRHHYDRO file: {str(e)}")
            raise

    def get_ManningN_dict(self):
        """
        Get the dictionary for Manning's n

        Returns
        -------

        """

        return self.srhhydro_content["ManningsN"]

    def modify_ManningsN(self, materialIDs, newManningsNValues, ManningN_MaterialNames):
        """Modify materialID's Manning's n value to new value

        Parameters
        ----------
        materialIDs : list
            material ID list
        newManningsNValues : list
            new Manning's n value list
        ManningN_MaterialNames : list
            name of materials list

        Returns
        -------

        """

        if gVerbose: print("Modify Manning's n value ...")

        if not isinstance(materialIDs[0], int):
            raise Exception("Material ID has to be an integer. The type of materialID passed in is ", type(materialIDs[0]))

        if not isinstance(newManningsNValues[0], float):
            raise Exception("Manning's n has to be a float. The type of newManningsNValue passed in is ", type(newManningsNValues[0]))


        nDict = self.srhhydro_content["ManningsN"]

        for i in range(len(materialIDs)):
            materialID = materialIDs[i]

            if materialID in nDict:
                if gVerbose: print("    Old Manning's n value =", nDict[materialID], "for material ID = ", materialID)
                nDict[materialID] = newManningsNValues[i]
                if gVerbose: print("    New Manning's n value =", nDict[materialID], "for material ID = ", materialID)
            else:
                print("The specified materialID", materialID, "is not in the Manning's n list. Please check.")

    def get_BC(self):
        """
        Get boundary condition dictionary.

        Returns
        -------

        """

        BC_Dict = self.srhhydro_content["BC"]

        return BC_Dict

    def get_InletQ(self):
        """
        Get InletQ dictionary:

        Returns
        -------

        """

        if 'IQParams' not in self.srhhydro_content:
            IQParams_Dict = {}   #empty
        else:
            IQParams_Dict = self.srhhydro_content['IQParams']

        return IQParams_Dict

    def modify_InletQ(self, bcIDs, newInletQValues):
        """Modify bcID's InletQ value to new value.

        Currently only constant inlet flow is supported.

        Parameters
        ----------
        bcIDs : list
            bc ID list
        newInletQValues : list
            new InletQ value list

        Returns
        -------

        """

        if gVerbose: print("Modify InletQ value ...")

        if not isinstance(bcIDs[0], int):
            raise Exception("Boundary condition ID has to be an integer. The type of bcID passed in is ",
                            type(bcIDs[0]))

        if not isinstance(newInletQValues[0], float):
            raise Exception("IneltQ value has to be a float. The type of newInletQValues passed in is ",
                            type(newInletQValues[0]))

        BC_Dict = self.srhhydro_content["BC"]

        if 'IQParams' not in self.srhhydro_content:
            raise Exception("There is no INLET-Q boundary in the SRHHydro file.")

        IQParams_Dict = self.srhhydro_content['IQParams']    #key in IQParams is the boundary ID (1-based)

        for i in range(len(bcIDs)):
            bcID = bcIDs[i]

            if bcID in BC_Dict:
                # Make sure this bcID's corresponding BC type is InletQ
                if BC_Dict[bcID] != "INLET-Q":
                    raise Exception("The boundary condition for the specified bcID ", bcID, "is not INLET-Q.")

                # Also make sure bcID is in the "IQParams" dictionary
                if bcID not in IQParams_Dict:
                    raise Exception("The specified bcID ", bcID, "is not in the IQParams dictionary.")

                #print("bcID = ", bcID)
                #print("IQParams_Dict[bcID] = ", IQParams_Dict[bcID])
                #print("newInletQValues[i] = ", newInletQValues[i])
                
                if gVerbose: print("    Old InletQ value =", IQParams_Dict[bcID]['discharge'], "for boundary ID = ", bcID)
                IQParams_Dict[bcID]['discharge'] = newInletQValues[i]
                if gVerbose: print("    New InletQ value =", IQParams_Dict[bcID]['discharge'], "for boundary ID = ", bcID)
            else:
                print("The specified bcID", bcID, "is not in the boundary list. Please check.")

    def get_ExitH(self):
        """
        Get ExitH dictionary:

        Returns
        -------

        """

        if 'EWSParamsC' not in self.srhhydro_content:
            EWSParamsC_Dict = {}   #empty
        else:
            EWSParamsC_Dict = self.srhhydro_content['EWSParamsC']

        return EWSParamsC_Dict

    def modify_ExitH(self, bcIDs, newExitHValues):
        """Modify bcID's ExitH value to new value.

        Currently only constant WSE is supported.

        Parameters
        ----------
        bcIDs : list
            bc ID list
        newExitHValues : list
            new ExitH value list

        Returns
        -------

        """

        if gVerbose: print("Modify Exit-H value ...")

        if not isinstance(bcIDs[0], int):
            raise Exception("Boundary condition ID has to be an integer. The type of bcID passed in is ",
                            type(bcIDs[0]))

        if not isinstance(newExitHValues[0], float):
            raise Exception("Exit-H value has to be a float. The type of newExitHValues passed in is ",
                            type(newExitHValues[0]))

        BC_Dict = self.srhhydro_content["BC"]

        if 'EWSParamsC' not in self.srhhydro_content:
            raise Exception("There is no constant stage exit boundary in the SRHHydro file.")

        EWSParamsC_Dict = self.srhhydro_content['EWSParamsC']

        for i in range(len(bcIDs)):
            bcID = bcIDs[i]

            if bcID in BC_Dict:
                # Make sure this bcID's corresponding BC type is InletQ
                if BC_Dict[bcID] != "EXIT-H":
                    raise Exception("The boundary condition for the specified bcID ", bcID, "is not EXIT-H.")

                # Also make sure bcID is in the "EWSParamsC" dictionary
                if bcID not in EWSParamsC_Dict:
                    raise Exception("The specified bcID ", bcID, "is not in the EWSParamsC dictionary.")

                print("bcID = ", bcID)
                print("EWSParamsC_Dict[bcID] = ", EWSParamsC_Dict[bcID])
                print("newExitHValues[i] = ", newExitHValues[i])

                if gVerbose: print("    Old EXIT-H value =", EWSParamsC_Dict[bcID]['discharge'], "for boundary ID = ", bcID)
                EWSParamsC_Dict[bcID]['discharge'] = newExitHValues[i]
                if gVerbose: print("    New EXIT-H value =", EWSParamsC_Dict[bcID][0], "for boundary ID = ", bcID)
            else:
                print("The specified bcID", bcID, "is not in the boundary list. Please check.")

    def modify_Case_Name(self, newCaseName):
        """Modify grid file name

        Parameters
        ----------
        newCaseName : str
            new case name

        Returns
        -------

        """

        if gVerbose: print("Modify case name ...")

        self.srhhydro_content['Case'] = newCaseName

    def get_Case_Name(self):
        return self.srhhydro_content['Case']

    def modify_Grid_FileName(self, newGridFileName):
        """Modify grid file name

        Parameters
        ----------
        newGridFileName : str
            new grid file name

        Returns
        -------

        """

        if gVerbose: print("Modify grid file name ...")

        self.srhhydro_content['Grid'] = newGridFileName

    def get_Grid_FileName(self):
        return self.srhhydro_content['Grid']

    def modify_HydroMat_FileName(self, newHydroMatFileName):
        """Modify grid file name

        Parameters
        ----------
        newHydroMatFileName : str
            new HydroMat file name

        Returns
        -------

        """

        if gVerbose: print("Modify HydroMat file name ...")

        self.srhhydro_content['HydroMat'] = newHydroMatFileName

    def get_HydroMat_FileName(self):
        return self.srhhydro_content['HydroMat']

    def save_as(self, new_srhhydro_file_name=None):
        """Save as a new SRHHydro file (useful for modification of the SRHHydro file)

        Parameters
        -------
        new_srhhydro_file_name : str
            name of the new srhhydro file. If not provided, the original filename will be used.
        """

        if new_srhhydro_file_name is None:
            new_srhhydro_file_name = self.srhhydro_filename

        if gVerbose: print("Writing the SRHHYDRO file %s \n" % new_srhhydro_file_name)

        try:
            fid = open(new_srhhydro_file_name, 'w')
        except IOError:
            print('srhhydro file open error')
            sys.exit()

        # Write header and basic parameters in order
        fid.write(f"SRHHYDRO {self.srhhydro_content.get('Version', '')}\n")
        
        # Write quoted string parameters
        for param in ['Case', 'Description']:
            if param in self.srhhydro_content:
                fid.write(f'{param} "{self.srhhydro_content[param]}"\n')
            else:
                raise ValueError(f"Parameter {param} not found in the SRHHydro content.")

        # Write basic simulation parameters
        for param in ['RunType', 'ModelTemp', 'UnsteadyOutput']:
            if param in self.srhhydro_content:
                fid.write(f'{param} {self.srhhydro_content[param]}\n')
            else:
                raise ValueError(f"Parameter {param} not found in the SRHHydro content.")

        # Write SimTime
        if 'SimTime' in self.srhhydro_content:
            value = self.srhhydro_content['SimTime']
            fid.write(f'SimTime {value[0]} {value[1]} {value[2]}\n')
        else:
            raise ValueError("Parameter SimTime not found in the SRHHydro content.")

        # Write turbulence and initial condition parameters
        for param in ['TurbulenceModel', 'ParabolicTurbulence', 'InitCondOption']:
            if param in self.srhhydro_content:
                fid.write(f'{param} {self.srhhydro_content[param]}\n')

                #if the initial condition option is WSE, need to write the InitZoneParams
                if param == 'InitCondOption' and self.srhhydro_content['InitCondOption'] == 'WSE':
                    fid.write(f'NumInitZones {self.srhhydro_content["NumInitZones"]}\n')
                    if 'InitZoneParams' in self.srhhydro_content:
                        for zone_id, params in sorted(self.srhhydro_content['InitZoneParams'].items()):
                            fid.write(f'InitZoneParams {zone_id} {params["p1"]} {params["p2"]} {params["p3"]} {params["unit"]}\n')
                    else:
                        raise ValueError("Parameter InitZoneParams not found in the SRHHydro content.")
            else:
                raise ValueError(f"Parameter {param} not found in the SRHHydro content.")

        # Write grid and material files
        for param in ['Grid', 'HydroMat']:
            if param in self.srhhydro_content:
                fid.write(f'{param} "{self.srhhydro_content[param]}"\n')
            else:
                raise ValueError(f"Parameter {param} not found in the SRHHydro content.")

        # Write monitor point file (monitoring points are optional)
        for param in ['MonitorPtFile']:
            if param in self.srhhydro_content:
                fid.write(f'{param} "{self.srhhydro_content[param]}"\n')
            else:
                print(f"Parameter {param} not found in the SRHHydro content.")

        # Write output parameters
        if 'OutputFormat' in self.srhhydro_content:
            fid.write(f'OutputFormat {self.srhhydro_content["OutputFormat"]} {self.srhhydro_content.get("OutputUnit", "")}\n')
        else:
            raise ValueError("Parameter OutputFormat not found in the SRHHydro content.")

        for param in ['OutputOption', 'OutputInterval']:
            if param in self.srhhydro_content:
                fid.write(f'{param} {self.srhhydro_content[param]}\n')
            else:
                raise ValueError(f"Parameter {param} not found in the SRHHydro content.")

        # Write Manning's n values
        if 'ManningsN' in self.srhhydro_content:
            for material_id, value in sorted(self.srhhydro_content['ManningsN'].items()):
                fid.write(f'ManningsN {material_id} {value}\n')
        else:
            raise ValueError("Parameter ManningsN not found in the SRHHydro content.")

        # Write boundary conditions
        if 'BC' in self.srhhydro_content:
            for bc_id, bc_type in sorted(self.srhhydro_content['BC'].items()):
                fid.write(f'BC {bc_id} {bc_type}\n')
        else:
            raise ValueError("Parameter BC not found in the SRHHydro content.")

        # Write inlet and exit parameters
        if 'EWSParamsC' in self.srhhydro_content:
            #print("EWSParamsC = ", self.srhhydro_content['EWSParamsC'])
            for bc_id, params in sorted(self.srhhydro_content['EWSParamsC'].items()):
                fid.write(f'EWSParamsC {bc_id} {params["wse"]} {params["unit"]} {params["type"]}\n')

        if 'EWSParamsRC' in self.srhhydro_content:
            #print("EWSParamsRC = ", self.srhhydro_content['EWSParamsRC'])
            for bc_id, params in sorted(self.srhhydro_content['EWSParamsRC'].items()):
                print("EWSParamsRC = ", bc_id, params)
                #fid.write(f'EWSParamsRC {bc_id} "{params[0]}" {params[1]} {params[2]}\n')
                fid.write(f'EWSParamsRC {bc_id} "{params["curve_file"]}" {params["unit"]} {params["type"]}\n')

        if 'IQParams' in self.srhhydro_content:
            #print("IQParams = ", self.srhhydro_content['IQParams'])
            for bc_id, params in sorted(self.srhhydro_content['IQParams'].items()):
                fid.write(f'IQParams {bc_id} {params["discharge"]} {params["unit"]} {params["type"]}\n')
                
        # Write pressure nodestrings
        if 'PressureNodestrings' in self.srhhydro_content:
            for group_id, nodestrings in sorted(self.srhhydro_content['PressureNodestrings'].items()):
                fid.write(f'PressureNodestrings {group_id} {" ".join(map(str, nodestrings))}\n')

        # Write pressure parameters and overtop settings
        if 'PressureParams' in self.srhhydro_content:
            for group_id, params in sorted(self.srhhydro_content['PressureParams'].items()):
                fid.write(f'PressureParams {group_id} {" ".join(map(str, params))}\n')
                
                # Write PressOvertop immediately after its corresponding PressureParams
                if 'PressOvertop' in self.srhhydro_content and group_id in self.srhhydro_content['PressOvertop']:
                    overtop_value = self.srhhydro_content['PressOvertop'][group_id]
                    fid.write(f'PressOvertop {group_id} {overtop_value}\n')
                    
                    # If overtop is 1, write corresponding weir parameters
                    if overtop_value == 1 and 'PressWeirParams' in self.srhhydro_content and group_id in self.srhhydro_content['PressWeirParams']:
                        weir_params = self.srhhydro_content['PressWeirParams'][group_id]
                        fid.write(f'PressWeirParams2 {group_id} {" ".join(map(str, weir_params))}\n')

        fid.close()

    def get_simulation_start_end_time(self):
        """Get the start and end time of simulation (in hours, from the "SimTime" entry)

        Returns
        -------
        startTime : float
            start time in hours
        endTime : float
            end time in hours

        """

        value = self.srhhydro_content["SimTime"]

        return value[0],value[2]

    def get_simulation_time_step_size(self):
        """Get the time step size of simulation (in seconds, from the "SimTime" entry)

        Returns
        -------
        deltaT : float
            time step size in seconds

        """

        value = self.srhhydro_content["SimTime"]

        return value[1]

    def get_grid_file_name(self):
        """Get grid file name (should be like "xxx.srhgeom")

        Returns
        -------
        value : str
            Name of the grid file

        """

        value = self.srhhydro_content["Grid"]

        return value

    def get_mat_file_name(self):
        """Get material file name (should be like "xxx.srhmat")

        Returns
        -------
        value : str
            Name of the material file

        """

        value = self.srhhydro_content["HydroMat"]

        return value

    # Add getter methods
    def get_pressure_overtop(self, nodestring_id=None):
        """Get pressure overtop parameters
        
        Args:
            nodestring_id: Optional nodestring ID. If None, returns all parameters
            
        Returns:
            dict: Pressure overtop parameters
        """
        if nodestring_id is None:
            return self.srhhydro_content.get('PressOvertop', {})
        return self.srhhydro_content.get('PressOvertop', {}).get(nodestring_id)

    def get_pressure_weir_params(self, nodestring_id=None):
        """Get pressure weir parameters
        
        Args:
            nodestring_id: Optional nodestring ID. If None, returns all parameters
            
        Returns:
            dict: Pressure weir parameters
        """
        if nodestring_id is None:
            return self.srhhydro_content.get('PressWeirParams', {})
        return self.srhhydro_content.get('PressWeirParams', {}).get(nodestring_id)

class SRH_2D_SRHGeom:
    """A class to handle srhgeom file for SRH-2D

    Notes:
        1. In SMS (SRH-2D), monitoring lines (ML) and monitoring points (MP) are treated separately and differently.
           ML is created in the same way as boundary lines. Nodes closest to the line are recorded as a NodeString
           in SRHGEOM file. MP is created as a feature point and its coordinates are recorded in the "srhmpoint" file.
           Interpolation is done to probe the results at MPs.

           This potentially have some drawbacks:

           a. The MLs and MPs have to be defined before the simulation is run, not afterwards. One can use plot over a
              line or point, but it is not convenient and it is limited. For example, to calculate flow over a line,
              one has to cut a line and then integrate.
           b. The way MLs are defined if slightly limited by the mesh (node locations). To be exactly at a ML, one has
              to interpolate the simulated solution to the points on MLs.

        2. Because BC and ML are mixed together, we need to separate them. However, this can only be done with the BC
           information in the srhhydro file. That is why we need to pass in bcDict in the constructor.

    Attributes
    ----------

    Methods
    ----------

    """

    def __init__(self, srhgeom_filename, bcDict):
        """Constructor for SRH_2D_SRHGeom

        Parameters
        ----------
        srhgeom_filename : str
            Name of the SRHGeom file
        bcDict : dict
            Boundary condtion dictionary
        """

        self.srhgeom_filename = srhgeom_filename

        self.Name = ''
        self.GridUnit = ''

        #boundary condition dictionary. In srhgeom file, not all nodeStrings are boundary conditions. Some of them
        #could be monitoring lines. This bcDict is passed in to distinguish the two (srhgeom does not have this
        #information; only srhhydro has). An SRH_2D_SRHHydro object has to be previously created, and then
        #bcDict <- srhhydro_obj.srhhydro_content["BC"].
        self.bcDict = bcDict

        #mesh information:
        #   number of elements
        self.numOfElements = -1

        #   number of nodes
        self.numOfNodes = -1

        # number of NodeString
        self.numOfNodeStrings = -1

        # First read of the SRHGEOM file to get number of elements and nodes
        self.getNumOfElementsNodes()

        # list of nodes for all elements
        self.elementNodesList = np.zeros([self.numOfElements, gMax_Nodes_per_Element], dtype=int)

        # number of nodes for each element (3,4,...,gMax_Nodes_per_Element)
        self.elementNodesCount = np.zeros(self.numOfElements, dtype=int)

        # list of edges for all elements (For each element, the number of nodes and the number of edges are the same.
        # Thus, there is no need for elementEdgesCount.
        self.elementEdgesList = np.zeros([self.numOfElements, gMax_Nodes_per_Element], dtype=int)

        # each element's vtk cell type
        self.vtkCellTypeCode = np.zeros(self.numOfElements, dtype=int)

        # each node's 3D coordinates
        self.nodeCoordinates = np.zeros([self.numOfNodes, 3], dtype=np.float64)

        # twoD mesh's boundingBox: [xmin, ymin, zmin, xmax, ymax, zmax]
        self.twoDMeshBoundingbox = []

        # each element (cell)'s bed elevation (z)
        self.elementBedElevation = np.zeros(self.numOfElements, dtype=np.float64)

        # center of each element
        self.elementCenters = np.zeros([self.numOfElements, 3], dtype=np.float64)

        # bed slope of each element at cell center
        self.elementBedSlope = np.zeros([self.numOfElements, 2], dtype=np.float64)

        # bed slope of each node
        self.nodeBedSlope = np.zeros([self.numOfNodes, 2], dtype=np.float64)

        # each NodeString's list of nodes (stored in a dictionary)
        self.nodeStringsDict = {}

        # get the mesh information (elementNodesList, elementNodesCount, vtkCellTypeCode, and nodeCoordinates)
        # from reading the SRHGEOM file again
        self.readSRHGEOMFile()

        # list of elements for all nodes
        self.nodeElementsList = np.zeros([self.numOfNodes, gMax_Elements_per_Node], dtype=int)

        # number of elements for each node (1,2,...,gMax_Elements_per_Node)
        self.nodeElementsCount = np.zeros(self.numOfNodes, dtype=int)

        # build the element list for each node
        self.buildNodeElements()

        #edges (connecting two nodes). Not using "faces" because they are lines. When 2D mesh is extruded to 3D,
        #these lines will become faces.
        self.edges = {} #dictionary: {[node1, node2]: edgeID}
        self.edges_r = {} #dictionary: {edgeID: [node1, node2]}  #revise of self.edges

        self.edgeElements = {} #dictionary: {edgeID: [element list]}. Here since one edge can be shared by two
                               #elements at the most, the element list's length is either 1 or 2

        #edges of each boundary. Only contains real boundaries, not monitoring lines.
        #boundaryID = nodeString ID read in from srhgeom file
        self.boundaryEdges = {} #dictionary: {boundaryID: [list of edge IDs]}
        self.boundaryEdges_normal = {} #dictionary: {boundaryID: [list of edge normal vectors]}. The normal vector is pointing outwards.

        #list of all boundary edge IDs (all lumped to one list)
        self.allBoundaryEdgeIDs = []

        #build "lines" and "boundaryLines" dictionaries. SRH-2D has no such information in SRHGEOM. We need to build it.
        self.buildEdgesAndBoundaryEdges()

        #calculate the bed slope of each element and then interpolate the bed slope to each node
        self.calculate_bed_slope()


    def getNumOfElementsNodes(self):
        """ Get the number of elements and nodes in srhgeom mesh file

        Returns


        """

        if gVerbose: print("Getting numbers of elements and nodes from the SRHGEOM file ...")

        # read the "srhgeom" mesh file
        try:
            #print("Opening the SRHGEOM file %s" % self.srhgeom_filename)
            srhgeomfile = open(self.srhgeom_filename, 'r')
        except:
            raise Exception('Failed openning the SRHGEOM file %s' % self.srhgeom_filename)

        count = 0
        elemCount = 0
        nodeCount = 0
        nodeStringCount = 0

        while True:
            count += 1

            # Get next line from file
            line = srhgeomfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

            if len(search) != 0:
                if search[0] == "Elem":
                    elemCount += 1
                    # print("Elem # %d: %s" % (elemCount, line))
                elif search[0] == "Node":
                    nodeCount += 1
                    # print("Node # %d: %s" % (nodeCount, line))
                elif search[0] == "NodeString":
                    nodeStringCount += 1

        srhgeomfile.close()

        self.numOfElements = elemCount
        self.numOfNodes = nodeCount
        self.numOfNodeStrings = nodeStringCount

        if gVerbose: print("There are %d elements, %d nodes, and %d node strings in the mesh." % (self.numOfElements,
                                                                                    self.numOfNodes, self.numOfNodeStrings))


    def readSRHGEOMFile(self):
        """ Get mesh information by reading the srhgeom file

        Parameters
        ----------

        Returns
        -------

        """

        if gVerbose: print("Reading the SRHGEOM file ...")

        # read the "srhgeom" mesh file
        try:
            srhgeomfile = open(self.srhgeom_filename, 'r')
        except:
            raise Exception('Failed openning srhgeom file %s' % self.srhgeom_filename)

        count = 0
        elemCount = 0
        nodeCount = 0
        nodeStringCount = 0

        # nodeString could be long enough to have more than one line. Use this to record the current nodeString during reading.
        currentNodeStringID = -1

        while True:
            count += 1

            # Get next line from file
            line = srhgeomfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

            if len(search) != 0:
                if search[0] == "Elem":
                    elemCount += 1
                    # print("Elem # %d: %s" % (elemCount, line))
                    self.elementNodesList[elemCount - 1][0:len(search[2:])] = [int(i) for i in search[2:]]
                    self.elementNodesCount[elemCount - 1] = len(search[2:])
                    if len(search[2:]) < 1 or len(search[2:]) > gMax_Nodes_per_Element:
                        sys.exit("Number of nodes for element %d is less than 1 or larger than the max of %d." % (elemCount,gMax_Nodes_per_Element))
                    self.vtkCellTypeCode[elemCount - 1] = vtkCellTypeMap[len(search[2:])]

                    #print("elem ID = ", elemCount-1)
                    #print("elementNodesList = ", self.elementNodesList[elemCount - 1])
                elif search[0] == "Node":
                    nodeCount += 1
                    # print("Node # %d: %s" % (nodeCount, line))
                    self.nodeCoordinates[nodeCount - 1] = [float(i) for i in search[2:]]
                elif search[0] == "NodeString":
                    nodeStringCount += 1
                    self.nodeStringsDict[int(search[1])] = [int(i) for i in search[2:]]
                    currentNodeStringID = int(search[1])
                elif search[0].lower() == "Name".lower():
                    self.Name = search[1]
                elif search[0].lower() == "GridUnit".lower():
                    self.GridUnit = search[1]
                elif ("SRHGEOM".lower() not in search[0].lower()) : #assume this is still nodeString
                    node_list = self.nodeStringsDict[currentNodeStringID]
                    for i in search:
                        node_list.append(int(i))
                    self.nodeStringsDict[currentNodeStringID] = node_list

        srhgeomfile.close()

        #calculate cell center and the bed elevation at cell center
        for cellI in range(self.numOfElements):

            cell_center = np.array([0.0, 0.0, 0.0])
            elev_temp = 0.0

            #loop over all nodes of current element
            for nodeI in range(self.elementNodesCount[cellI]):
                cell_center += self.nodeCoordinates[self.elementNodesList[cellI][nodeI]-1]
                elev_temp += self.nodeCoordinates[self.elementNodesList[cellI][nodeI]-1,2]  #z elevation; nodeI-1 because node number is 1-based for SRH-2D

            if self.elementNodesCount[cellI] > 0:
                self.elementBedElevation[cellI] = elev_temp / self.elementNodesCount[cellI]
                self.elementCenters[cellI] = cell_center / self.elementNodesCount[cellI]
            else:
                print("element ", cellI, " has no nodes. Exiting...")
                sys.exit()

        #calculate the 2D mesh's bounding box
        xmin = np.min(self.nodeCoordinates[:,0])
        ymin = np.min(self.nodeCoordinates[:,1])
        zmin = np.min(self.nodeCoordinates[:,2])
        xmax = np.max(self.nodeCoordinates[:,0])
        ymax = np.max(self.nodeCoordinates[:,1])
        zmax = np.max(self.nodeCoordinates[:,2])

        self.twoDMeshBoundingbox = [xmin, ymin, zmin, xmax, ymax, zmax]

        if gVerbose: print("2D mesh's bounding box = ", self.twoDMeshBoundingbox)

        if False:
            print("elementNodesList = ", self.elementNodesList)
            print("elementNodesCount = ", self.elementNodesCount)
            print("vtkCellTypeCode = ", self.vtkCellTypeCode)
            print("nodeCoordinates = ", self.nodeCoordinates)
            print("elementBedElevation = ", self.elementBedElevation)
            print("nodeStrings = ", self.nodeStringsDict)

    def calculate_bed_slope(self):
        """ 
        Calculate the bed slope of each element and then interpolate the bed slope to each node. The slope of each element (e.g., triangle and quadrilateral) is calculated
        by using the Gauss theorem. For each element, the bed slope is calculated by:
        1. compute the outward normal vector of each edge of the element.
        2. interpolate the bed elevation at the center of each edge of the element from the the two points of the edge.
        3. use the Gauss theorem to calculate the bed slope of the element: S = -(1/A) * \int_{\partial A} z * n * ds, where A is the area of the element,
           z is the bed elevation, n is the outward normal vector, and ds is the length of the edge.
        4. interpolate the bed slope of the element to each node of the element.

        Returns
        -------
        
        """
        import numpy as np
        import sys

        # Initialize array to store bed slopes [Sx, Sy] for each element
        num_elements = self.numOfElements
        element_bed_slopes = np.zeros((num_elements, 2))

        #print("self.elementNodesCount = ", self.elementNodesCount)

        for i in range(num_elements):
            # Get nodes for this element
            element_nodes = self.elementNodesList[i]
            num_nodes = self.elementNodesCount[i]
            
            # Get coordinates of element nodes
            node_coords = self.nodeCoordinates[element_nodes-1]  #-1 because node number is 1-based for SRH-2D
            
            # Calculate element area and determine node ordering
            if num_nodes == 3:  # Triangle
                # Area of triangle using cross product
                v1 = node_coords[1] - node_coords[0]
                v2 = node_coords[2] - node_coords[0]
                area = 0.5 * np.cross(v1[:2], v2[:2])  # Remove abs() to preserve sign
            elif num_nodes == 4:  # Quadrilateral
                # Area of quadrilateral using shoelace formula
                x = node_coords[:num_nodes, 0]
                y = node_coords[:num_nodes, 1]
                area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)  # Remove abs() to preserve sign

                # Check for degenerate quads
                #if np.any(np.diff(x) == 0) and np.any(np.diff(y) == 0):
                #    print(f"Warning: Element {i} may be degenerate")

            else:  # General polygon (pentagon, hexagon, etc.)
                # Area of general polygon using shoelace formula
                x = node_coords[:num_nodes, 0]
                y = node_coords[:num_nodes, 1]
                area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)  # Remove abs() to preserve sign
            #else:
            #    print("Element %d is not a triangle or quadrilateral. Other shapes are not supported yet. Exiting..." % i)
            #    sys.exit()

            # Take absolute value of area for size, but keep sign for orientation
            area_abs = np.abs(area)
            is_clockwise = area < 0

            # Check for degenerate elements (zero or very small area)            
            # if area_abs < 1e-12:  # Very small threshold for numerical precision
            #     print(f"ERROR: Element {i} has zero area ({area})")
            #     print(f"Element nodes: {element_nodes}")
            #     print(f"Node coordinates:")
            #     for j, coord in enumerate(node_coords):
            #         print(f"  Node {j}: {coord}")
            #     print(f"Number of nodes: {num_nodes}")
                
            #     # For now, assign a very small area to avoid division by zero
            #     # This is a temporary fix - the real issue should be investigated
            #     area_abs = 1e-12
            #     area = 1e-12 if area >= 0 else -1e-12
            #     print(f"Assigned temporary area: {area}")

            #debug
            #print("element id = ", i)
            #print("element nodes = ", element_nodes)
            #print("num_nodes = ", num_nodes)
            #print("area = ", area)
            #print("is_clockwise = ", is_clockwise)

            # Initialize slope components
            Sx = 0.0
            Sy = 0.0

            # Process each edge of the element (number of edges = number of nodes)
            for j in range(num_nodes):
                # Get nodes of current edge
                node1 = element_nodes[j]
                node2 = element_nodes[(j + 1) % num_nodes]
                
                # Get coordinates of edge nodes
                p1 = self.nodeCoordinates[node1-1]   #-1 because node number is 1-based for SRH-2D
                p2 = self.nodeCoordinates[node2-1]   #-1 because node number is 1-based for SRH-2D
                
                # Calculate edge length
                edge_length = np.sqrt(np.sum((p2[:2] - p1[:2])**2))
                
                # Calculate outward normal vector (normalized)
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                # For counterclockwise elements, normal = [dy, -dx]
                # For clockwise elements, normal = [-dy, dx]
                if is_clockwise:
                    normal = np.array([-dy, dx]) / edge_length
                else:
                    normal = np.array([dy, -dx]) / edge_length

                #debug
                #print("edge id = ", j)
                #print("normal = ", normal)
                #print("edge center = ", (p1 + p2) / 2)
                
                # Interpolate bed elevation at edge center
                z_center = 0.5 * (p1[2] + p2[2])
                
                # Add contribution to slope components using Gauss theorem
                Sx += z_center * normal[0] * edge_length
                Sy += z_center * normal[1] * edge_length

            # Final slope calculation (note the negative sign from the Gauss theorem)
            element_bed_slopes[i] = [-Sx / area_abs, -Sy / area_abs]

        self.elementBedSlope = element_bed_slopes
        
        # Interpolate bed slope to each node in the mesh
        for i in range(self.numOfNodes):
            # Get elements for this node
            node_elements = self.nodeElementsList[i]
            num_elements = self.nodeElementsCount[i]

            # Get slopes of the elements that share this node
            node_slopes = element_bed_slopes[node_elements]

            # Interpolate the slope components separately
            node_slope_x = np.sum(node_slopes[:, 0]) / num_elements
            node_slope_y = np.sum(node_slopes[:, 1]) / num_elements

            # Store the slope vector
            self.nodeBedSlope[i] = [node_slope_x, node_slope_y]

    def save_as(self, newSrhgeomFileName, dir=''):
        """ Save as a srhgeom file.

        It can be used to save to a new srhgeom file after things have been modified, e.g., bathymetry

        Parameters
        ----------

        Returns
        -------

        """

        if gVerbose: print("Wring the SRHGEOM file %s \n" % newSrhgeomFileName)

        if len(dir) != 0:
            newSrhgeomFileName = dir + "/" + newSrhgeomFileName

        try:
            fid = open(newSrhgeomFileName, 'w')
        except IOError:
            print('srhgeom file open error')
            sys.exit()

        fid.write('SRHGEOM 30\n')
        fid.write('Name \"Saved from pyHMT2D \"\n')

        fid.write('\n')

        fid.write('GridUnit %s \n' % self.GridUnit)

        #write elements
        for elemI in range(1, self.numOfElements + 1):
            fid.write('Elem %d' % elemI)
            for nodeI in range(1, self.elementNodesCount[elemI-1] + 1):
                fid.write(' %d' % self.elementNodesList[elemI-1, nodeI-1])
            fid.write('\n')

        #write nodes
        for nodeI in range(1, self.numOfNodes + 1):
            fid.write('Node %d %f %f %f\n' % (nodeI,
                                            self.nodeCoordinates[nodeI-1, 0],
                                            self.nodeCoordinates[nodeI-1, 1],
                                            self.nodeCoordinates[nodeI-1, 2]))

        #write nodeStrings
        for nodeStringID, nodeStringNodeList in self.nodeStringsDict.items():
            fid.write('NodeString %d' % nodeStringID)

            # line break counter (start a new line every 10 nodes)
            line_break_counter = 0

            for nodeI in nodeStringNodeList:
                fid.write(' %d' % nodeI)

                line_break_counter += 1

                # 10 numbers per line
                if ((line_break_counter % 10) == 0):
                    fid.write("\n")

                    line_break_counter = 0

            fid.write('\n')

        fid.close()


    def buildNodeElements(self):
        """ Build node's element list for all nodes

        Returns
        -------

        """

        if gVerbose: print("Building mesh's node, elements, and topology ...")

        #create an empty list of lists for all nodes
        self.nodeElementsList = [[] for _ in range(self.numOfNodes)]

        #loop over all cells
        for cellI in range(self.numOfElements):
            #loop over all nodes of current cell
            for i in range(self.elementNodesCount[cellI]):
                #get the node ID of current node
                nodeID = self.elementNodesList[cellI][i]

                self.nodeElementsList[nodeID-1].append(cellI)  #nodeID-1 because SRH-2D is 1-based

        #count beans in each node's element list
        for nodeI in range(self.numOfNodes):
            self.nodeElementsCount[nodeI] = len(self.nodeElementsList[nodeI])

        #print("nodeElementsCount = ", self.nodeElementsCount)
        #print("nodeElementsList = ", self.nodeElementsList)

    def buildEdgesAndBoundaryEdges(self):
        """ Build edges and boundaryEdges dictionaries

        Returns
        -------

        """

        current_edgeID = 1  #note: SRH-2D is 1-based.

        #loop over all elements
        for cellI in range(self.numOfElements):
            # loop over all edges of current element
            for i in range(self.elementNodesCount[cellI]):
                # get the node ID of current and next nodes
                # connecting current edge

                nodeID_1 = 0
                nodeID_2 = 0

                if i != (self.elementNodesCount[cellI]-1):  #if not the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][i+1]
                else:                                      #if the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][0]

                curr_edge_node_IDs = tuple(sorted([nodeID_1, nodeID_2]))

                #check whether the current edge is in the edges dictionary or not
                if curr_edge_node_IDs not in self.edges:  #if no, added a new edge
                    self.edges[curr_edge_node_IDs] = current_edgeID
                    self.edges_r[current_edgeID] = curr_edge_node_IDs

                    #add the new edge to edgeElements dictionary too
                    self.edgeElements[current_edgeID] = [cellI+1]  #+1 because SRH-2D is 1-based

                    current_edgeID += 1
                else: # if yes, neighbor element already has that edge and no need to do anything, but need to add
                      # the current element to edgeElements list

                    #get the existing edge's ID
                    existingEdgeID = self.edges[curr_edge_node_IDs]

                    #in this case, there should be one and only one element in the list
                    assert len(self.edgeElements[existingEdgeID]) == 1

                    self.edgeElements[existingEdgeID].append(cellI+1) #+1 because SRH-2D is 1-based


        #build the allBoundaryEdgeIDs list: boundary edges are those with only one element
        #loop over all edges in the mesh
        for edge in self.edgeElements:
            if len(self.edgeElements[edge]) == 1:
                self.allBoundaryEdgeIDs.append(edge)

        #list to flag whether a boundary edge in self.allBoundaryEdgeIDs is used or not
        #if not, it is a default boundary (wall)
        allBoundaryEdgeUsageFlag = {}
        for boundaryEdgeID in self.allBoundaryEdgeIDs:
            allBoundaryEdgeUsageFlag[boundaryEdgeID] = False

        #now all edges are built, we check the boundary edges
        #loop through all boundaries.
        for nodeString in self.nodeStringsDict:

            #not all NodeStrings are boundaries. Could be monitoring lines. Need to exclude them.
            if nodeString not in self.bcDict.keys():
                continue

            #also exclude internal BC nodeStrings such as weir and pressure
            if 'WEIR' in self.bcDict[nodeString] or 'PRESSURE' in self.bcDict[nodeString]:
                continue

            #list of nodes in current nodeString
            nodeString_nodeList = self.nodeStringsDict[nodeString]

            current_boundary_edge_list = []

            #loop through all edges in current boundary
            for i in range(len(nodeString_nodeList)-1):
                nodeID_1 = nodeString_nodeList[i]
                nodeID_2 = nodeString_nodeList[i + 1]

                #curr_edge_node_IDs = tuple(sorted([nodeID_1, nodeID_2]))
                curr_edge_node_IDs = tuple([nodeID_1, nodeID_2])
                curr_edge_node_IDs_reverse = tuple([nodeID_2, nodeID_1])

                #check whether the edge is in the edge list. If not, we got a problem.
                if curr_edge_node_IDs in self.edges:
                    current_boundary_edge_list.append(self.edges[curr_edge_node_IDs])
                    allBoundaryEdgeUsageFlag[self.edges[curr_edge_node_IDs]] = True
                elif curr_edge_node_IDs_reverse in self.edges:
                    current_boundary_edge_list.append(-self.edges[curr_edge_node_IDs_reverse])
                    allBoundaryEdgeUsageFlag[self.edges[curr_edge_node_IDs_reverse]] = True
                else:
                    print("Boundary edge ", curr_edge_node_IDs, "in NodeString", nodeString, "can not be found in edge list. Mesh is wrong. Exiting...")
                    sys.exit()

            self.boundaryEdges[nodeString] = current_boundary_edge_list

        #build default boundary as wall (in SMS, the default boundary is not exported in SRHGEOM)
        #boundary edges are those who is used only by one element.
        unusedBoundaryEdgeList = []
        for boundaryEdgeID in self.allBoundaryEdgeIDs:
            if not allBoundaryEdgeUsageFlag[boundaryEdgeID]: #if not used
                unusedBoundaryEdgeList.append(boundaryEdgeID)

        if len(unusedBoundaryEdgeList) > 0:  #if there are unused boundary edges
            defaultWallBoundaryID = len(self.boundaryEdges)+1   #boundary ID for default boundary
            self.boundaryEdges[defaultWallBoundaryID] = unusedBoundaryEdgeList

        #compute the outward normal vectors of all boundary edges
        for boundaryID in self.boundaryEdges:
            self.boundaryEdges_normal[boundaryID] = []
            for edgeID in self.boundaryEdges[boundaryID]:
                #get the nodes of the current edge
                edge_nodes = self.edges_r[abs(edgeID)]   #edgeID is negative for the reverse edge
                node1 = edge_nodes[0]
                node2 = edge_nodes[1]
                p1 = self.nodeCoordinates[node1-1]
                p2 = self.nodeCoordinates[node2-1]

                edge_center = (p1 + p2) / 2

                edge_vector = p2 - p1
                edge_vector_magnitude = np.linalg.norm(edge_vector)
                edge_normal_vector = edge_vector / edge_vector_magnitude

                #make sure the normal vector is pointing outwards: get the cell/element that the edge belongs to, 
                #get the cell center, 
                element_ID = self.edgeElements[abs(edgeID)][0]
                cell_center = self.elementCenters[element_ID-1]

                #compute the vector from the cell center to the edge node
                cell_center_vector = edge_center - cell_center

                #compute the dot product of the cell center vector and the edge normal vector
                dot_product = np.dot(cell_center_vector, edge_normal_vector)

                #if the dot product is positive, the normal vector is pointing outwards
                if dot_product < 0:
                    edge_normal_vector = -edge_normal_vector

                self.boundaryEdges_normal[boundaryID].append(edge_normal_vector)


        #build elementEdgesList
        # loop over all elements
        for cellI in range(self.numOfElements):
            # loop over all edges of current element
            for i in range(self.elementNodesCount[cellI]):    #i is the current edge index for this element
                # get the node ID of current and next nodes
                # connecting current edge

                nodeID_1 = 0
                nodeID_2 = 0

                if i != (self.elementNodesCount[cellI] - 1):  # if not the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][i + 1]
                else:  # if the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][0]

                curr_edge_node_IDs = tuple([nodeID_1, nodeID_2])           #the original order of the edge nodes from mesh file
                curr_edge_node_IDs_reverse = tuple([nodeID_2, nodeID_1])   #the order of the edge nodes reversed

                # check whether the current edge is in the edges dictionary or not
                if curr_edge_node_IDs in self.edges:
                    self.elementEdgesList[cellI, i] = self.edges[curr_edge_node_IDs]
                elif curr_edge_node_IDs_reverse in self.edges:
                    self.elementEdgesList[cellI, i] = -self.edges[curr_edge_node_IDs_reverse]   #negative means the order is opposite.
                else:  # something is wrong if the edge is not found in the edge list.
                    print("Cell: ", cellI, " , edge: ", i, "  edge node IDs: ", curr_edge_node_IDs, "can not be found in edge list. Mesh is wrong. Exiting...")
                    sys.exit()

        #debug
        if False:
            print("edges", self.edges)
            print("edgeElements",self.edgeElements)
            print("elementEdgesList", self.elementEdgesList)
            print("allBoundaryEdgeIDs",self.allBoundaryEdgeIDs)
            print("boundaryEdges", self.boundaryEdges)

    def extrude_to_3D(self, layerHeights, mshFileName, bTerrainFollowing=False, dir=''):
        """ Extrude the 2D mesh to 3D with layers

        When extruded: node -> line, edge -> face, 2D element-> volume, boundary line -> boundary face



        The extruded mesh is written to a GMSH MSH file. Tried MSH version 4. But it seems
        too complicated (entities etc.). They are not necessary for our purpose.
        Therefore, we will stick with version 2 and OpenFOAM takes this version well.

        Ref: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029

        Parameters
        ----------
        layerHeights : list
            a list with strictly increasing heights for the layers. If bTerrainFollowing
            is True, only the last value is used as the top elevation.
        mshFileName : str
            name of the GMSH MSH file to be written
        bTerrainFollowing : bool
            whether to follow the terrain. If False, the bottom will be at z = 0, and the layer heights
            are specified in layerHeights. If True, the bottom will key the terrain elevation, the top
            will be at layerHeights[-1] and interpolation inbetween.

        Returns
        -------

        """

        #GMSH element type code
        MSHTRI = 2  # 3-node triangle
        MSHQUAD = 3  # 4-node quadrilateral

        MSHHEX = 5  # 8-node hexahedron
        MSHPRISM = 6  # 6-node prism

        allNodes = {}   #dict: {nodeID: [x,y,z]}
        allCells = {}   #dict: {cellID: [cell_type, [list of nodes]]}
        number_of_prism = 0 #number of prism in the extruded 3D mesh
        number_of_hex = 0 #number of hex in the extruded 3D mesh

        #some sanity check nlayers and layerHeights
        nlayers = len(layerHeights)

        if nlayers == 0:
            raise Exception("layerHeight is empty. Please check. Exiting ...")

        #check the values in the list layerHeights are strictly increasing
        bHeighIncreasing = all(i < j for i, j in zip(layerHeights, layerHeights[1:]))

        if not bHeighIncreasing:
            raise Exception("Values in layerHeight are not strictly increasing. Please check. Exiting ...")

        #a list for all side boundaries (not including top and bottom)
        allSideBoundaryFaces = {} #dict: {boundaryID: [list of faces]}. Here face is a list of nodes.

        #create all nodes
        #1. add the original nodes in 2D mesh
        for nodeI in range(self.nodeCoordinates.shape[0]):
            if bTerrainFollowing:
                allNodes[nodeI] = self.nodeCoordinates[nodeI,:]
            else:
                allNodes[nodeI] = self.nodeCoordinates[nodeI, :]
                allNodes[nodeI][2] = 0.0  # make the 2D mesh flat

        #2. add the extruded nodes
        for layerI in range (1, nlayers+1):
            for nodeI in range(self.nodeCoordinates.shape[0]):
                if bTerrainFollowing: #interpolation between bottom elevation and the top elevation in layerHeights[-1]
                    elev = allNodes[nodeI][2] + (layerHeights[-1] - allNodes[nodeI][2])/nlayers*layerI
                    allNodes[nodeI + layerI * self.numOfNodes] = np.array([self.nodeCoordinates[nodeI, 0],
                                                                           self.nodeCoordinates[nodeI, 1],
                                                                           elev])
                else:
                    # make the original 2D mesh flat; only add the layer height
                    allNodes[nodeI + layerI*self.numOfNodes] = np.array([self.nodeCoordinates[nodeI,0],
                                                                         self.nodeCoordinates[nodeI,1],
                                                                         layerHeights[layerI-1]])

        #create all cells

        #bottom and top boundary faces list
        bottomBoundaryFaces = [None] * self.numOfElements
        topBoundaryFaces = [None] * self.numOfElements

        cellID = 0

        #layer by layer
        for layerI in range(1, nlayers+1):
            #loop through each element in 2D mesh
            for elementI in range(self.numOfElements):

                if self.elementNodesCount[elementI] == 3: #triangle -> prism
                    cell_type = MSHPRISM
                    number_of_prism += 1
                elif self.elementNodesCount[elementI] == 4: #quad -> hex
                    cell_type = MSHHEX
                    number_of_hex += 1

                cellID += 1

                cell_list = []

                top_face_node_list = []
                bottom_face_node_list = []

                #loop over all nodes of the current elment in 2D mesh (at current layer's bottom)
                for nodeI in range(self.elementNodesCount[elementI]):
                    cell_list.append(self.elementNodesList[elementI, nodeI] + (layerI-1)*self.numOfNodes)

                    #if this is the first layer, record the bottom face to bottomBoundaryFaces list
                    if layerI == 1:
                        bottom_face_node_list.append(self.elementNodesList[elementI, nodeI])

                #now add the nodes at current layer's top.
                for nodeI in range(self.elementNodesCount[elementI]):
                    cell_list.append(self.elementNodesList[elementI, nodeI] + layerI*self.numOfNodes)

                    #if this is the last layer, record the top face to topBoundaryFaces list
                    if layerI == nlayers:
                        top_face_node_list.append(self.elementNodesList[elementI, nodeI] + layerI*self.numOfNodes)

                allCells[cellID] = [cell_type, cell_list]

                bottomBoundaryFaces[elementI] = bottom_face_node_list
                topBoundaryFaces[elementI] = top_face_node_list


        #create all boundaries

        #allSideBoundaryFaces = {} #dict: {boundaryID: [list of faces]}. Here face is a list of nodes.

        #loop through all boundary lines (nodeStrings) in 2D mesh. self.boundaryEdges only
        #contains real boundaries, not monitoring lines.
        for boundaryID in self.boundaryEdges:

            #list of all faces for current boundary
            #each face is a list four nodes (quad because we extrude a line in vertical direction)
            faceList = []

            #loop through all edges of current boundary
            for edgeID in self.boundaryEdges[boundaryID]:
                nodeID_1 = self.edges_r[edgeID][0]
                nodeID_2 = self.edges_r[edgeID][1]

                #loop through all layers
                for layerI in range(1, nlayers + 1):
                    #add the four nodes to the current face (ordered)
                    curFace = [nodeID_1+(layerI-1)*self.numOfNodes, nodeID_1+layerI*self.numOfNodes,
                               nodeID_2+ layerI*self.numOfNodes, nodeID_2+(layerI-1)*self.numOfNodes]

                    #add the current face to the faceList for current boundary
                    faceList.append(curFace)

            #now we have all faces in the current boundary
            allSideBoundaryFaces[boundaryID] = faceList


        #write to GMSH MSH file

        if dir!='':
            mshFileName_base = dir + '/' + mshFileName
        else:
            mshFileName_base = mshFileName
        if gVerbose: print("Write to GMESH MSH file with name = ", mshFileName_base)


        fid = open(mshFileName_base, 'w')

        #write MeshFormat
        fid.write("$MeshFormat\n")
        fid.write("2.1 0 8\n")
        fid.write("$EndMeshFormat\n")

        #write PhysicalNames:
        fid.write("$PhysicalNames\n")
        fid.write("%d\n" % (len(allSideBoundaryFaces)+2+1))   #need to output all side boundaries,
                                                              #top + bottom (2), and the volume (1)

        # maximum side boundary ID. This is the starting point for IDs of top, bottom, and volume
        maxSideBoundaryID = -999

        #loop through all side boundaries
        for boundaryID in allSideBoundaryFaces:
            if boundaryID in self.bcDict.keys(): #if the current BC is defined in srhhydro
                boundary_name = "boundary_"+str(boundaryID)+"_"+self.bcDict[boundaryID]  #name = boundary_ID_BCType
            else:                                #else, this should be default wall boundary
                boundary_name = "boundary_"+str(boundaryID)+"_"+"wall"                   #name = boundary_ID_wall

            fid.write("2 %d \"%s\"\n" % (boundaryID, boundary_name))  #dimension, physicalTag, name

            maxSideBoundaryID = max(maxSideBoundaryID, boundaryID)

        #output the top and bottom boundaries
        boundary_name = "top"
        fid.write("2 %d \"%s\"\n" % (maxSideBoundaryID+1, boundary_name))  # dimension, physicalTag, name

        boundary_name = "bottom"
        fid.write("2 %d \"%s\"\n" % (maxSideBoundaryID+2, boundary_name))  # dimension, physicalTag, name


        #output the volume
        volume_name = "channel"  #just a generic name for the volume
        fid.write("3 %d \"%s\"\n" % (maxSideBoundaryID+2+1, volume_name))

        fid.write("$EndPhysicalNames\n")

        #write Nodes
        fid.write("$Nodes\n")
        #total_num_nodes = (nlayers+1)*self.numOfNodes
        total_num_nodes = len(allNodes)

        fid.write("%d\n" % total_num_nodes)

        #output node coordinates
        for nodeI in allNodes:
            #GMSH is 1-based
            fid.write("%d %f %f %f\n" %(nodeI+1, allNodes[nodeI][0], allNodes[nodeI][1], allNodes[nodeI][2]))

        fid.write("$EndNodes\n")

        #write Elements (boundary conditions and volume)
        #  numElements = number of all faces in all boundaries + number of cells in all volumes

        elementTag_counter = 0 #count the number of elements

        fid.write("$Elements\n")

        if (number_of_hex == 0) and (number_of_prism == 0):
            print("There is on hex or prism, the only supported cell types, in the extruded 3D mesh. Check the mesh. Exiting ...")
            sys.exit()

        #output total number of elements
        #                    num. of side boundary faces + num. of faces on top and bottom + num. of hex + num. of prism
        fid.write("%d \n" % (len(self.allBoundaryEdgeIDs)*nlayers + self.numOfElements*2 + number_of_hex + number_of_prism))

        #loop through all side boundaries
        for boundaryID in allSideBoundaryFaces:
            #output each quad face's
            for faceID in range(len(allSideBoundaryFaces[boundaryID])):
                elementTag_counter += 1
                fid.write("%d 3 1 %d %d %d %d %d\n" % (elementTag_counter, boundaryID, allSideBoundaryFaces[boundaryID][faceID][0],
                                    allSideBoundaryFaces[boundaryID][faceID][1], allSideBoundaryFaces[boundaryID][faceID][2],
                                    allSideBoundaryFaces[boundaryID][faceID][3]))

        #output top and bottom boundaries
        for faceI in topBoundaryFaces:
            elementTag_counter += 1

            if len(faceI) == 3: #triangle
                fid.write("%d 2 1 %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID+1,
                                                     faceI[0],faceI[1],faceI[2]))
            elif len(faceI) == 4: #quad
                fid.write("%d 3 1 %d %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID+1,
                                                     faceI[0], faceI[1], faceI[2], faceI[3]))

        for faceI in bottomBoundaryFaces:
            elementTag_counter += 1

            if len(faceI) == 3:  # triangle
                fid.write("%d 2 1 %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID + 2,
                                                     faceI[0], faceI[1], faceI[2]))
            elif len(faceI) == 4:  # quad
                fid.write("%d 3 1 %d %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID + 2,
                                                        faceI[0], faceI[1], faceI[2], faceI[3]))

        #output the volume's tag and node list. Hex and prism need to be seperate
        #1.Hex
        if number_of_hex > 0:
            #loop through all cells in the volume
            for cellID in allCells:
                if allCells[cellID][0] == MSHHEX:
                    elementTag_counter += 1
                    fid.write("%d 5 1 %d %d %d %d %d %d %d %d %d\n" % (elementTag_counter, maxSideBoundaryID+2+1,
                                                                       allCells[cellID][1][0], allCells[cellID][1][1],
                                            allCells[cellID][1][2], allCells[cellID][1][3], allCells[cellID][1][4],
                                            allCells[cellID][1][5], allCells[cellID][1][6], allCells[cellID][1][7]))

        #2. prism
        if number_of_prism > 0:
            #loop through all cells in the volume
            for cellID in allCells:
                if allCells[cellID][0] == MSHPRISM:
                    elementTag_counter += 1
                    fid.write("%d 6 1 %d %d %d %d %d %d %d\n" % (elementTag_counter, maxSideBoundaryID+2+1,
                                                                 allCells[cellID][1][0], allCells[cellID][1][1],
                                                                 allCells[cellID][1][2], allCells[cellID][1][3],
                                                                 allCells[cellID][1][4], allCells[cellID][1][5]))

        fid.write("$EndElements")

        fid.close()


    def output_2d_mesh_to_vtk(self,meshVTKFileName, bFlat=False, dir=''):
        """output the 2D mesh to vtk


        Parameters
        --------
        meshVTKFileName : str
            VTK file name for the mesh
        bFlat : bool, optional
            whether to make a flat mesh (z=0)
        dir : str, optional
            directory to write the VTK file to

        Returns
        -------

        """

        vtkFileName = ''
        if len(dir) == 0:
            vtkFileName = meshVTKFileName
        else:
            vtkFileName = dir + "/" + meshVTKFileName

        # build VTK object:
        # points
        pointsVTK = vtk.vtkPoints()
        flatCoordinates = np.copy(self.nodeCoordinates)

        if bFlat:
            flatCoordinates[:,2] = 0.0

        pointsVTK.SetData(VN.numpy_to_vtk(flatCoordinates))

        # cell topology information list: [num. of nodes, node0, node1, .., num. of nodes, nodexxx]
        # the list start with the number of nodes for a cell and then the list of node indexes
        connectivity_list = []
        # type of cells (contains the number of face points
        cellFPCounts = np.zeros(self.elementNodesList.shape[0], dtype=np.int64)

        # loop over all elements
        for k in range(self.elementNodesList.shape[0]):
            connectivity_list.append(self.elementNodesCount[k])

            for nodeI in range(self.elementNodesCount[k]):
                connectivity_list.append(
                    self.elementNodesList[k][nodeI] - 1)  # -1 becasue SRH-2D is 1-based.

            cellFPCounts[k] = self.elementNodesCount[k]

        connectivity = np.array(connectivity_list, dtype=np.int64)

        # convert cell's number of face points to VTK cell type
        vtkHandler_obj = vtkHandler()
        cell_types = vtkHandler_obj.number_of_nodes_to_vtk_celltypes(cellFPCounts)

        cellsVTK = vtk.vtkCellArray()
        cellsVTK.SetCells(self.elementNodesList.shape[0], VN.numpy_to_vtkIdTypeArray(connectivity))

        uGrid = vtk.vtkUnstructuredGrid()
        uGrid.SetPoints(pointsVTK)
        uGrid.SetCells(cell_types.tolist(), cellsVTK)

        # write to vtk file
        unstr_writer = vtk.vtkUnstructuredGridWriter()
        unstr_writer.SetFileName(vtkFileName)
        unstr_writer.SetInputData(uGrid)
        unstr_writer.Write()


    def output_nodeString_line_coordinates(self, nodeStringID, nodeStringFileName, dir=''):
        """ Output the nodeString line coordinates into a text file


        Parameters
        ----------
        nodeStringID : int
            the ID of the nodeString to be written out
        nodeStringFileName : str
            the name of the file to be written to
        dir : str, optional
            directory to write the output file to

        Returns
        -------

        """

        nodeStringFileName_final = ''
        if len(dir) == 0:
            nodeStringFileName_final = nodeStringFileName
        else:
            nodeStringFileName_final = dir + "/" + nodeStringFileName

        #check whether the given nodeStringID is in the nodeStringDict
        if nodeStringID not in self.nodeStringsDict.keys():
            print("The given nodeStringID", nodeStringID, "is not valid. Valid nodeString IDs are: ",
                  self.nodeStringsDict.keys())

        try:
            fid = open(nodeStringFileName_final, 'w')
        except IOError:
            print('NodeString file open error')
            sys.exit()

        #write the header
        fid.write("station, x, y, z\n")

        # list of nodes in current nodeString
        nodeString_nodeList = self.nodeStringsDict[nodeStringID]

        station = [0] * len(nodeString_nodeList)
        x = [0] * len(nodeString_nodeList)
        y = [0] * len(nodeString_nodeList)
        z = [0] * len(nodeString_nodeList)

        #loop over all nodes in the list
        for nodeI in range(len(nodeString_nodeList)):
            x[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 0]
            y[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 1]
            z[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 2]

            if nodeI!=0: #if not the first node in the list, calculate the station
                station[nodeI] = station[nodeI-1] + np.sqrt( (x[nodeI]-x[nodeI-1])**2 + (y[nodeI]-y[nodeI-1])**2 )

        #write the data
        for nodeI in range(len(nodeString_nodeList)):
            fid.write("%f, %f, %f, %f\n" % (station[nodeI], x[nodeI], y[nodeI], z[nodeI]))

        fid.close()

    def output_as_gmsh(self, gmshFileName, dir=''):
        """ Output SRH-2D mesh in GMESH format

        Not implemented yet.

        Parameters
        ----------
        gmshFileName : str
            GMSH MSH file name to write to
        dir : str, optional
            directory name to write to

        Returns
        ---------

        """

        return NotImplemented

        if dir!='':
            gmshFileName_base = dir + '/' + gmshFileName
        else:
            gmshFileName_base = gmshFileName
        if gVerbose: print("GMESH file name = ", gmshFileName_base)


class SRH_2D_SRHMat:
    """A class to handle srhmat file for SRH-2D

    Attributes
    ----------
    srhmat_filename : str
        name of the srhmat file
    numOfMaterials : int
        number of materials defined in the srhmat file
    matZoneCells : dict
        each material zone's cell list (in a dictionary: material zone ID and cell list)


    Methods

    """
    def __init__(self, srhmat_filename):
        """SRH_2D_SRHMat constructor from a srhmat file

        Parameters
        ----------
        srhmat_filename : str
            name of the srhmat file
        """

        self.srhmat_filename = srhmat_filename

        #number of materials (Manning's n zones)
        self.numOfMaterials = -1

        #list of Material names (in a dictionary: ID and name)
        self.matNameList = {}

        #each material zone's cell list (in a dictionary: material zone ID and cell list)
        self.matZoneCells = {}

        #build material zones data
        self.buildMaterialZonesData()

    def buildMaterialZonesData(self):
        """Build the data for material zones from reading the srhmat file

        Returns
        -------

        """

        if gVerbose: print("Reading the SRHMAT file ...")

        # read the "srhmat" material file
        try:
            srhmatfile = open(self.srhmat_filename, 'r')
        except:
            print('Failed openning SRHMAT file', self.srhmat_filename)
            sys.exit()

        res_MatNames = {} #dict to hold "MatName" entries for material names list: ID and name
        res_Materials = {} #dict to hold "Material" information: ID and list of cells

        #add default material name
        res_MatNames["-1"] = "Default"

        current_MaterialID = -1
        current_MaterialCellList = []

        while True:
            # Get a line from file
            line = srhmatfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

#            print("search[0] == Material", (search[0] == "Material"))

            if len(search) != 0: #if it is not just an empty line
                if search[0] == "SRHMAT":
                    continue
                elif search[0] == "NMaterials":
                    self.numOfMaterials = int(search[1])
                    continue
                elif search[0] == "MatName":
                    res_MatNames[search[1]] = search[2]
                    continue
                elif search[0] == "Material": # a new material zone starts
                    #if this is not the first Material zone; save the previous Material zone
                    if ((current_MaterialID != -1) and (len(current_MaterialCellList)!=0)):
                        res_Materials[current_MaterialID] = copy.deepcopy(current_MaterialCellList)

                        #clear up current_MaterialID and current_MaterialCellList
                        current_MaterialID = -1
                        current_MaterialCellList.clear()

                    current_MaterialID = int(search[1])
                    current_MaterialCellList.extend(search[2:])
                else: #still in a material zone (assume there is no other things other than those listed above in the
                      #srhmat file)
                    current_MaterialCellList.extend(search)

            #add the last material zone
            res_Materials[current_MaterialID] = copy.deepcopy(current_MaterialCellList)

        srhmatfile.close()

        #within SRHMat file, we have no idea whether there are any cells in the default material zone. Add a default
        #wiht an empty cell list. The user of SRH_2D_SRHMat should take care of that.
        res_Materials[-1] = []

        #convert material ID from str to int
        for zoneID, cellList in res_Materials.items():
            #print("cellList ", cellList)
            temp_list = [int(i) for i in cellList]
            res_Materials[zoneID] = temp_list

        self.matNameList = res_MatNames
        self.matZoneCells = res_Materials

        #print("matNameList = ", self.matNameList)
        #print("matZoneCells = ", self.matZoneCells)


    def find_cell_material_ID(self, cellID):
        """Given a cell ID, find its material (Manning's n) ID

        Parameters
        ----------
        cellID : int
            cell ID to find material ID for

        Returns
        -------

        """

        if not self.matZoneCells:
            self.buildMaterialZonesData()

        #whether the cell is found in the list
        bFound = False

        for matID, cellList in self.matZoneCells.items():
            if (cellID+1) in cellList:  #cellID+1 because SRH-2D is 1-based
                bFound = True
                return matID

        #if cellID is not found, report problem.
        #This could be due to several reasons. Most likely, the material coverage in SMS does not
        #cover certain area of the mesh. In this case, SMS will not report the material for the missed
        #cells in srhmat (no warning given either). Here, if we can't find the cell in material list,
        #simply use the default Manning's n ID.
        if not bFound:
            print("pyHMT2D: In find_cell_material_ID(cellID), cellID =", cellID, "is not found. Check mesh and material coverage. Default is used.")
            return 0 #return the default (?)
            #sys.exit()

    def save_as(self, new_srhmat_file_name):
        """Save as a new SRHMat file (useful for modification of the SRHMAT file)

        Parameters
        -------
        new_srhmat_file_name : str
            name of the new srhmat file

        """

        if gVerbose: print("Wring the SRHMAT file %s \n" % new_srhmat_file_name)

        try:
            fid = open(new_srhmat_file_name, 'w')
        except IOError:
            print('srhmat file open error')
            sys.exit()

        fid.write('SRHMAT 30\n')
        fid.write('NMaterials %d\n' % (self.numOfMaterials))

        for matID, matName in self.matNameList.items():
            fid.write('MatName %s %s\n' % (matID, matName))

        element_counter = 0

        for matID, cellList in self.matZoneCells.items():

            #don't write out the default
            if matID==-1:
                continue

            fid.write('Material %d' % (matID))

            for cellI in cellList:
                element_counter += 1
                fid.write(' %d' % cellI)

                if element_counter==10:
                    element_counter = 0
                    fid.write('\n')

        fid.close()


class SRH_2D_Data(HydraulicData):
    """
    A class for SRH-2D data I/O, manipulation, and format conversion
    
    This class is designed to read SRH-2D results in format other than VTK. It can 
    save SRH-2D results into VTK format for visualization in Paraview, parse
    SRH-2D mesh information and convert/save to other formats, query SRH-2D
    results (e.g., for calibration), etc.

    Between srhdryo and srhsif, only one of them is needed depending on the
    type of the input file.
    
    Attributes
    ----------
    hdf_filename : str
        The name for the HDF result file generated by HEC-RAS
    srhhydro_obj : SRH_2D_SRHHydro
        An object to hold information in the SRHHYDRO file
    srhsif_obj : SRH_2D_SIF
        An object to hold information in the SIF file
    srhgeom_obj : SRH_2D_SRHGeom
        An object to hold information in the SRHGEOM file
    srhmat_obj : SRH_2D_SRHMat
        An object to hold information in the SRHMAT file

    Methods
    -------

    
    """
    
    def __init__(self, srhcontrol_filename):
        """Constructor for SRH_2D_Data

        srhgeom_filename and srhmat_filename are contained in the SRHHydro or SRH-2D SIF file.

        Parameters
        ----------
        srhcontrol_filename : str
            Name of the SRH-2D control file: either SRHHydro or SIF file

        """

        HydraulicData.__init__(self, "SRH-2D")

        #make sure the control file is either SRHHydro or SIF
        if not srhcontrol_filename.endswith(".srhhydro") and \
           not ("sif" in srhcontrol_filename.lower()):  #contains "sif" or "SIF"
            raise ValueError("SRH-2D control file must be either SRHHydro or SIF file")

        #check if the control file is a SRHHydro or SIF file
        self.control_type = "SRHHydro" if srhcontrol_filename.endswith(".srhhydro") else "SIF"

        self.srhcontrol_filename = srhcontrol_filename

        #extract path in srhcontrol_filename. We assume the SRHGeom and SRHMat files
        #are located in the same directory as the srhcontrol_filename file, which can be
        #either a srhhydro or sif file. In the control file,
        #the file names for SRHGeom and SRHMat do not contain the directory.
        file_path, _, = os.path.split(srhcontrol_filename)

        #read SRH_2D_SRHHydro file and build SRH_2D_SRHHydro object
        self.srhhydro_obj = None
        self.srhsif_obj = None
        if self.control_type == "SRHHydro":
            self.srhhydro_obj = SRH_2D_SRHHydro(self.srhcontrol_filename)
        else:
            self.srhsif_obj = SRH_2D_SIF(self.srhcontrol_filename)

        #get the srhgeom_filename and srhmat_filename from srhhydro_obj
        self.srhgeom_filename = None
        self.srhmat_filename = None    #Note: it seems SRH-2D preprocess assumes the file name for srhmat is the same as srhgeom. For example, if srhgeom is "Muncie.srhgeom", then srhmat is "Muncie.srhmat".
        if self.control_type == "SRHHydro":
            self.srhgeom_filename =      self.srhhydro_obj.get_grid_file_name() if len(file_path) == 0 \
                                    else file_path+"/"+self.srhhydro_obj.get_grid_file_name()

            self.srhmat_filename =       self.srhhydro_obj.get_mat_file_name() if len(file_path) == 0 \
                                    else file_path+"/"+self.srhhydro_obj.get_mat_file_name()
        else:
            self.srhgeom_filename = self.srhsif_obj.srhsif_content["mesh_file_name"]
            self.srhmat_filename = self.srhsif_obj.srhsif_content["mesh_file_name"].replace(".srhgeom", ".srhmat")


        if gVerbose:
            print("self.srhgeom_filename = ", self.srhgeom_filename)
            print("self.srhmat_filename = ", self.srhmat_filename)

        #read and build SRH_2D_Geom, and SRH_2D_Mat objects
        self.srhgeom_obj = None
        self.srhmat_obj = None
        if self.control_type == "SRHHydro":
            self.srhgeom_obj = SRH_2D_SRHGeom(self.srhgeom_filename, self.srhhydro_obj.srhhydro_content["BC"])
            self.srhmat_obj = SRH_2D_SRHMat(self.srhmat_filename)
        else:
            self.srhgeom_obj = SRH_2D_SRHGeom(self.srhgeom_filename, self.srhsif_obj.srhsif_content["BC"])
            self.srhmat_obj = SRH_2D_SRHMat(self.srhmat_filename)

        #Manning's n value at cell centers and nodes
        self.ManningN_cell = np.zeros(self.srhgeom_obj.numOfElements, dtype=float)
        self.ManningN_node = np.zeros(self.srhgeom_obj.numOfNodes, dtype=float)

        #build Manning's n value at cell centers and nodes
        self.buildManningN_cells_nodes()

        #XMDF data: if desired, the result data in a XMDF file
        #can be read in and stored in the following arrays.
        #Depending on whether the XMDF is nodal or at cell centers, use different arrays.
        self.xmdfTimeArray_Nodal = None   #numpy array to store all time values
        self.xmdfAllData_Nodal = {}     #dict to store {varName: varValues (numpy array)}

        self.xmdfTimeArray_Cell = None   #numpy array to store all time values
        self.xmdfAllData_Cell = {}     #dict to store {varName: varValues (numpy array)}

    def get_case_name(self):
        """Get the "Case" name for this current case ("Case" in srhhyro file)

        Returns
        -------
        case_name : str
            Case name

        """
        if self.control_type == "SRHHydro":
            return self.srhhydro_obj.srhhydro_content["Case"]
        else:
            return self.srhsif_obj.srhsif_content["Case"]

    def get_output_format_unit(self):
        """ Get the format and unit of the output

        Returns
        -------
        
        """

        if self.control_type == "SRHHydro":
            return self.srhhydro_obj.srhhydro_content["OutputFormat"], self.srhhydro_obj.srhhydro_content["OutputUnit"]
        else:
            return self.srhsif_obj.srhsif_content["OutputFormat"], self.srhsif_obj.srhsif_content["OutputUnit"]

    def buildManningN_cells_nodes(self):
        """ Build Manning's n values at cell centers and nodes


        This calculation is based on the srhhydro and srhgeom files, not interpolation from a GeoTiff file. This is the
        Manning's n value used by SRH-2D.

        Returns
        -------

        """

        if gVerbose: print("Building Manning's n values for cells and nodes in SRH-2D mesh ...")

        if self.control_type == "SRHHydro":
            #Manning's n dictionary in srhhydro
            nDict = self.srhhydro_obj.srhhydro_content["ManningsN"]
        elif self.control_type == "SIF":
            nDict = self.srhsif_obj.srhsif_content["ManningN"]

        if gVerbose: print("nDict = ", nDict)


        #loop over all cells in the mesh
        for cellI in range(self.srhgeom_obj.numOfElements):
            #get the material ID of current cell
            matID = self.srhmat_obj.find_cell_material_ID(cellI)

            #print("matID = ", matID)
            #print("nDict[matID] = ", nDict[matID])

            #get the Manning's n value for current cell
            self.ManningN_cell[cellI] = nDict[matID]

        #interpolate Manning's n from cell centers to nodes
        #loop over all nodes in the mesh
        for nodeI in range(self.srhgeom_obj.numOfNodes):
            n_temp = 0.0

            #loop over all cells that share the current node
            for i in range(self.srhgeom_obj.nodeElementsCount[nodeI]):
                cellI = self.srhgeom_obj.nodeElementsList[nodeI][i]
                n_temp += self.ManningN_cell[cellI]

            #take the average
            self.ManningN_node[nodeI] = n_temp/self.srhgeom_obj.nodeElementsCount[nodeI]
        

    def output_boundary_manning_n_profile(self, nodeStringID, nodeStringFileName, dir=''):
        """ Output Manning's n profile for a specified boundary

        Parameters
        ----------
        nodeStringID : int
            nodeString ID of specified boundary
        nodeStringFileName : str
            name of output file
        dir : str, optional
            directory for the output file

        Returns
        -------

        """

        nodeStringFileName_final = ''
        if len(dir) == 0:
            nodeStringFileName_final = nodeStringFileName
        else:
            nodeStringFileName_final = dir + "/" + nodeStringFileName

        #check whether the given nodeStringID is in the nodeStringDict
        if nodeStringID not in self.srhgeom_obj.nodeStringsDict.keys():
            print("The given nodeStringID", nodeStringID, "is not valid. Valid nodeString IDs are: ",
                  self.srhgeom_obj.nodeStringsDict.keys())
            sys.exit()

        try:
            fid = open(nodeStringFileName_final, 'w')
        except IOError:
            print('Boundary\'s Manning n file open error')
            sys.exit()

        #write the header
        fid.write("station, n\n")

        # list of node ID (1-based) in current nodeString
        nodeString_nodeList = self.srhgeom_obj.nodeStringsDict[nodeStringID]

        station = [0] * len(nodeString_nodeList)
        x = [0] * len(nodeString_nodeList)
        y = [0] * len(nodeString_nodeList)
        z = [0] * len(nodeString_nodeList)

        #loop over all nodes in the list
        for nodeI in range(len(nodeString_nodeList)):
            x[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 0]
            y[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 1]
            z[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 2]

            if nodeI!=0: #if not the first node in the list, calculate the station
                station[nodeI] = station[nodeI-1] + np.sqrt( (x[nodeI]-x[nodeI-1])**2 + (y[nodeI]-y[nodeI-1])**2 )

        #write the data
        for nodeI in range(len(nodeString_nodeList)):
            fid.write("%f, %f\n" % (station[nodeI], self.ManningN_node[nodeString_nodeList[nodeI] - 1])) #-1 because SRH-2D is 1-based

        fid.close()

    def get_ManningN_dict(self):
        """ Get the Manning's n dictionary from the SRH-2D control file

        Returns
        -------
        
        """

        if self.control_type == "SRHHydro":
            return self.srhhydro_obj.srhhydro_content["ManningsN"]
        else:
            return self.srhsif_obj.srhsif_content["ManningN"]
        
    def modify_ManningsNs(self, materialIDs, newManningsNValues, materialNames):
        """ Modify the Manning's n value for the specified material IDs

        Parameters
        ----------
        materialIDs : list
            list of material IDs
        newManningsNValues : list
            list of new Manning's n values
        materialNames : list
            list of material names

        Returns
        -------

        """

        if gVerbose: print("Modifying Manning's n values for the specified material IDs ...")

        if self.control_type == "SRHHydro":
            self.srhhydro_obj.modify_ManningsN(materialIDs, newManningsNValues, materialNames)
        else:
            self.srhsif_obj.modify_ManningsN(materialIDs, newManningsNValues, materialNames)


        #update the ManningN at cells and nodes after the modification
        self.buildManningN_cells_nodes()


    def modify_InletQ(self, inlet_q_bc_IDs, new_inlet_q_values):
        """ Modify the inlet flow rate for the specified boundary IDs


        Parameters
        ----------
        inlet_q_bc_IDs : list   
            list of inlet flow rate boundary IDs
        new_inlet_q_values : list
            list of new inlet flow rate values

        Returns
        -------

        """

        if gVerbose: print("Modifying inlet flow rate for the specified boundary IDs ...")

        if self.control_type == "SRHHydro":
            self.srhhydro_obj.modify_InletQ(inlet_q_bc_IDs, new_inlet_q_values)
        else:
            self.srhsif_obj.modify_InletQ(inlet_q_bc_IDs, new_inlet_q_values)

    def modify_ExitH(self, exit_h_bc_IDs, new_exit_h_values):
        """ Modify the exit water surface elevation for the specified boundary IDs

        Parameters
        ----------
        exit_h_bc_IDs : list
            list of exit water surface elevation boundary IDs
        new_exit_h_values : list
            list of new exit water surface elevation values

        Returns
        -------

        """

        if gVerbose: print("Modifying exit water surface elevation for the specified boundary IDs ...")

        if self.control_type == "SRHHydro":
            self.srhhydro_obj.modify_ExitH(exit_h_bc_IDs, new_exit_h_values)
        else:
            self.srhsif_obj.modify_ExitH(exit_h_bc_IDs, new_exit_h_values)

    def readSRHXMDFFile(self, xmdfFileName, bNodal):


        """ Read SRH-2D result file in XMDF format. The XMDF file can be either nodal or cell center. It also contains results from multiple time steps.

        Parameters
        ----------
        xmdfFileName : str
            file name for the XMDF data
        bNodal : bool
            whether it is nodal (True) or at cell center (False)

        Returns
        -------


        """

        #for debug
        #dumpXMDFFileItems(xmdfFileName)

        if bNodal:
            if "XMDFC" in xmdfFileName:
                if gVerbose: print("Warning: bNodal indices the result is for nodal values. However, the XMDF file name contains "
                      "\"XMDFC\", which indicates it is for cell center values. Please double check.")
        else:
            if "XMDFC" not in xmdfFileName:
                if gVerbose: print("Warning: bNodal indices the result is for cell centers. However, the XMDF file name "
                      "does not contain \"XMDFC\", which indicates it is for nodal values. Please double check.")

        if gVerbose: print("Reading the XMDF file ...\n")

        xmdfFile = h5py.File(xmdfFileName, "r")

        #build the list of solution variables
        varNameList = []

        #copy of the original velocity vector variable name
        varNameVelocity = ''

        for ds in xmdfFile.keys():
            #print(ds)

            if ds != "File Type" and ds != "File Version":
                if "Velocity" in ds:
                    varNameVelocity = '%s' % ds
                    vel_x = ds.replace("Velocity","Vel_X",1)
                    vel_y = ds.replace("Velocity","Vel_Y",1)
                    varNameList.append(vel_x)
                    varNameList.append(vel_y)
                else:
                    varNameList.append(ds)

        if not varNameList: #empty list
            print("There is no solution varialbes in the XMDF file. Exiting ...")
            sys.exit()

        if gVerbose: print("Variables in the XMDF file: ", varNameList)

        #build the 1d array of time for the solution
        #It seems all solution variables share the same time (although
        #each has their own recrod of "Time"). So only get one is sufficient.
        if bNodal:
            self.xmdfTimeArray_Nodal = np.array(xmdfFile[varNameList[0]]['Times'])
            if gVerbose: print("Time values for the solutions: ", self.xmdfTimeArray_Nodal)
        else:
            self.xmdfTimeArray_Cell = np.array(xmdfFile[varNameList[0]]['Times'])
            if gVerbose: print("Time values for the solutions: ", self.xmdfTimeArray_Cell)

        #result for all varialbes and at all times in a dictionary (varName : numpy array)
        for varName in varNameList:
            if gVerbose: print("varName =", varName)

            varName_temp = varName

            if ("Vel_X" in varName) or ("Vel_Y" in varName):  # if it is velocity components
                varName = varNameVelocity                     # use the original velocity variable name

            if bNodal: #for nodal data
                if "Velocity" in varName:
                    if "Vel_X" in varName_temp: #vel-x
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Nodal[varName_temp] = np.array(xmdfFile[varName]['Values'][:,:,0])
                    else:
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Nodal[varName_temp] = np.array(xmdfFile[varName]['Values'][:,:,1])
                else:
                    self.xmdfAllData_Nodal[varName] = np.array(xmdfFile[varName]['Values'])

                    #fixe the water elevation = -999 in SRH-2D
                    if varName == "Water_Elev_ft" or varName == "Water_Elev_m":
                        for nodeI in range(len(self.xmdfAllData_Nodal[varName])):
                            #loop over time steps
                            for timeI in range(self.xmdfAllData_Nodal[varName].shape[0]):
                                if self.xmdfAllData_Nodal[varName][timeI, nodeI] == -999:
                                    self.xmdfAllData_Nodal[varName][timeI, nodeI] = self.srhgeom_obj.nodeCoordinates[nodeI, 2] #use node elevation


                if np.array(xmdfFile[varName]['Values']).shape[1] != self.srhgeom_obj.numOfNodes:
                    print("The number of nodes in the XMDF file (%d) is different from that in the mesh (%d). "
                          "Abort readSRHXMDFFile(...) function."
                          % (np.array(xmdfFile[varName]['Values']).shape[1], self.srhgeom_obj.numOfNodes))
                    return

                #print(self.xmdfAllData_Nodal)
            else: #for cell center data
                if "Velocity" in varName:
                    if "Vel_X" in varName_temp: #vel-x
                        #print("varName_temp, varName = ", varName_temp, varName)
                        self.xmdfAllData_Cell[varName_temp] = np.array(xmdfFile[varName]['Values'])[:,:,0]
                    else: #vel-y
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Cell[varName_temp] = np.array(xmdfFile[varName]['Values'])[:,:,1]
                else:
                    self.xmdfAllData_Cell[varName] = np.array(xmdfFile[varName]['Values'])

                    #fixe the water elevation = -999 in SRH-2D
                    if varName == "Water_Elev_ft" or varName == "Water_Elev_m":
                        #loop through time
                        for timeI in range(self.xmdfAllData_Cell[varName].shape[0]):
                            for cellI in range(self.xmdfAllData_Cell[varName].shape[1]):
                                if self.xmdfAllData_Cell[varName][timeI, cellI] == -999:
                                    self.xmdfAllData_Cell[varName][timeI, cellI] = self.srhgeom_obj.elementBedElevation[cellI]

                if np.array(xmdfFile[varName]['Values']).shape[1] != self.srhgeom_obj.numOfElements:
                    print("The number of elements in the XMDF file (%d) is different from that in the mesh (%d). "
                          "Abort readSRHXMDFFile(...) function."
                          % (np.array(xmdfFile[varName]['Values']).shape[1], self.srhgeom_obj.numOfElements))

                #print(self.xmdfAllData_Cell)

        xmdfFile.close()

    def outputXMDFDataToVTK(self, bNodal, timeStep=-1,lastTimeStep=False, dir=''):
        """Output XMDF result data to VTK

        This has to be called after the XMDF data have been loaded by calling readSRHXMDFFile(...).


        Parameters
        ----------
            bNodal : bool
                whether export nodal data or cell center data. Currently, it can't output both.
            timeStep : int
                optionly only the specified time step will be saved
            lastTimeStep : bool
                optionally specify only the last time step
            dir : str, optional
                directory to write to

        Returns
        -------
            vtkFileNameList : list
                a list of vtkFileName (for each time step)

        """
        if gVerbose: print("Output all data in the XMDF file to VTK ...")

        if (bNodal and (not self.xmdfAllData_Nodal)) or ((not bNodal) and (not self.xmdfAllData_Cell)):
            print("Empty XMDF data arrays. Call readSRHXMDFFile() function first. Exiting ...")
            sys.exit()

        #check the sanity of timeStep
        if (timeStep != -1):
            if bNodal:
                if (not timeStep in range(len(self.xmdfTimeArray_Nodal))):
                    message = "Specified timeStep = %d not in range (0 to %d)." % (timeStep, len(self.xmdfTimeArray_Nodal))
                    sys.exit(message)
            else:
                if (not timeStep in range(len(self.xmdfTimeArray_Cell))):
                    message = "Specified timeStep = %d not in range (0 to %d)." % (timeStep, len(self.xmdfTimeArray_Cell))
                    sys.exit(message)

        vtkFileName_base = ''
        #print("self.srhhydro_obj.srhhydro_content[\"Case\"] = ", self.srhhydro_obj.srhhydro_content["Case"])

        if self.control_type == "SRHHydro":

            if dir!='':
                vtkFileName_base = dir + '/' + 'SRH2D_' + self.srhhydro_obj.srhhydro_content["Case"] #use the case name of the
                # SRH-2D case
            else:
                vtkFileName_base = 'SRH2D_' + self.srhhydro_obj.srhhydro_content["Case"]  # use the case name of the SRH-2D case

        else:
            if dir!='':
                vtkFileName_base = dir + '/' + 'SRH2D_' + self.srhsif_obj.srhsif_content["Case"] #use the case name of the
                # SRH-2D case
            else:
                vtkFileName_base = 'SRH2D_' + self.srhsif_obj.srhsif_content["Case"]  # use the case name of the SRH-2D case

        if gVerbose: print("vtkFileName_base = ", vtkFileName_base)

        #build result variable names
        resultVarNames = list(self.xmdfAllData_Nodal.keys()) if bNodal else list(self.xmdfAllData_Cell.keys())

        #add Manning's n (at cell center regardless bNodal)
        ManningNVarName = "ManningN"

        if bNodal:
            if gVerbose: print("All nodal solution variable names: ", resultVarNames[:-1], "and cell center ManningN")
        else:
            if gVerbose: print("All cell center solution variable names: ", resultVarNames)

        #add bed elevation at nodes regardless whether bNodal is True or False. Nodal elevation is more accurate representation
        #of the terrain because cell center elevation is averaged from nodal elevations.
        #get the units of the results
        if self.control_type == "SRHHydro":
            units = self.srhhydro_obj.srhhydro_content['OutputUnit']
        elif self.control_type == "SIF":
            units = self.srhsif_obj.srhsif_content['OutputUnit']
        else:
            raise NotImplementedError("TODO: implement the outputXMDFDataToVTK for SIF file")
        

        if gVerbose: print("SRH-2D result units:", units)
        bedElevationVarName = ''
        velocityVarName = ''
        if units == "SI":
            bedElevationVarName = "Bed_Elev_m"
            velocityVarName = "Velocity_m_p_s"
        else:
            bedElevationVarName = "Bed_Elev_ft"
            velocityVarName = "Velocity_ft_p_s"

        if gVerbose: print("Nodal bed elevation name: ", bedElevationVarName)

        #list of vtkFileName (to be returned to caller)
        vtkFileNameList = []

        #loop through each time step
        timeArray = self.xmdfTimeArray_Nodal if bNodal else self.xmdfTimeArray_Cell
        for timeI in range(timeArray.shape[0]):

            if lastTimeStep:
                if timeI < (timeArray.shape[0] - 1):
                    continue

            if (timeStep != -1) and (timeI != timeStep):
                continue

            if gVerbose: print("timeI and time = ", timeI, timeArray[timeI])

            #numpy array for all solution variables at one time
            resultData =      np.zeros((self.srhgeom_obj.numOfNodes, len(resultVarNames)), dtype="float32") if bNodal \
                         else np.zeros((self.srhgeom_obj.numOfElements, len(resultVarNames)), dtype="float32")

            str_extra = "_N_" if bNodal else "_C_"

            vtkFileName = vtkFileName_base + str_extra + str(timeI+1).zfill(4) + ".vtk"
            if gVerbose: print("vtkFileName = ", vtkFileName)

            #loop through each solution variable (except Bed_Elev and ManningN, which will be added seperately)
            for varName, varI in zip(resultVarNames, range(len(resultVarNames))):
                if gVerbose: print("varName = ", varName)
                #get the values of current solution varialbe at current time
                resultData[:,varI] =      self.xmdfAllData_Nodal[varName][timeI,:] if bNodal \
                                     else self.xmdfAllData_Cell[varName][timeI,:]

            #add Mannng's n to cell center
            #self.ManningN_cell

            #add nodal bed elevation to VTK's point data
            nodalElevation = self.srhgeom_obj.nodeCoordinates[:,2]

            #build VTK object:
            # points
            pointsVTK = vtk.vtkPoints()
            pointsVTK.SetData(VN.numpy_to_vtk(self.srhgeom_obj.nodeCoordinates))

            # cell topology information list: [num. of nodes, node0, node1, .., num. of nodes, nodexxx]
            # the list start with the number of nodes for a cell and then the list of node indexes
            connectivity_list = []

            # type of cells (contains the number of face points
            cellFPCounts = np.zeros(self.srhgeom_obj.elementNodesList.shape[0], dtype=np.int64)

            #loop over all elements
            for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
                connectivity_list.append(self.srhgeom_obj.elementNodesCount[k])

                for nodeI in range(self.srhgeom_obj.elementNodesCount[k]):
                    connectivity_list.append(self.srhgeom_obj.elementNodesList[k][nodeI]-1) #-1 becasue SRH-2D is 1-based.

                cellFPCounts[k] = self.srhgeom_obj.elementNodesCount[k]

            connectivity = np.array(connectivity_list, dtype=np.int64)

            # convert cell's number of face points to VTK cell type
            vtkHandler_obj = vtkHandler()
            cell_types = vtkHandler_obj.number_of_nodes_to_vtk_celltypes(cellFPCounts)

            cellsVTK = vtk.vtkCellArray()
            cellsVTK.SetCells(self.srhgeom_obj.elementNodesList.shape[0], VN.numpy_to_vtkIdTypeArray(connectivity))

            uGrid = vtk.vtkUnstructuredGrid()
            uGrid.SetPoints(pointsVTK)
            uGrid.SetCells(cell_types.tolist(), cellsVTK)

            cell_data = uGrid.GetCellData()  # This holds cell data
            point_data = uGrid.GetPointData()  # This holds point data

            # add solutions

            # column numbers for Vel_X and Vel_Y for vector assemble
            nColVel_X = -1
            nColVel_Y = -1

            # First output all solution variables as scalars
            if gVerbose: print('The following solution variables are processed and saved to VTK file: \n')
            for k in range(len(resultVarNames)):
                if gVerbose: print('     %s\n' % resultVarNames[k])

                #if it is a veloctiy component, only record its location
                #not output to VTK by itself, will be assembed to vector.
                if resultVarNames[k].find('Vel_X') != -1:
                    nColVel_X = k
                    continue
                elif resultVarNames[k].find('Vel_Y') != -1:
                    nColVel_Y = k
                    continue

                if bNodal:
                    temp_point_data_array = VN.numpy_to_vtk(resultData[:,k])
                    temp_point_data_array.SetName(resultVarNames[k])
                    point_data.AddArray(temp_point_data_array)
                else:
                    temp_cell_data_array = VN.numpy_to_vtk(resultData[:,k])
                    temp_cell_data_array.SetName(resultVarNames[k])
                    cell_data.AddArray(temp_cell_data_array)

            #add veloctiy by combining components
            currentTimePointV = np.zeros((self.srhgeom_obj.numOfNodes, 3)) if bNodal else \
                                np.zeros((self.srhgeom_obj.numOfElements,3))
            currentTimePointV[:, 0] = resultData[:,nColVel_X]
            currentTimePointV[:, 1] = resultData[:,nColVel_Y]
            currentTimePointV[:, 2] = 0.0

            if bNodal:
                temp_point_data_array = VN.numpy_to_vtk(currentTimePointV)
                temp_point_data_array.SetName(velocityVarName)
                point_data.AddArray(temp_point_data_array)
            else:
                temp_cell_data_array = VN.numpy_to_vtk(currentTimePointV)
                temp_cell_data_array.SetName(velocityVarName)
                cell_data.AddArray(temp_cell_data_array)

            #add Manning's n
            if bNodal:
                temp_point_data_array = VN.numpy_to_vtk(self.ManningN_node)
                temp_point_data_array.SetName(ManningNVarName)
                point_data.AddArray(temp_point_data_array)
            else:
                temp_cell_data_array = VN.numpy_to_vtk(self.ManningN_cell)
                temp_cell_data_array.SetName(ManningNVarName)
                cell_data.AddArray(temp_cell_data_array)

            #add nodal bed elevation
            currentBedElev = np.zeros(self.srhgeom_obj.numOfNodes) if bNodal else np.zeros(self.srhgeom_obj.numOfElements)

            currentBedElev = self.srhgeom_obj.nodeCoordinates[:,2]

            temp_point_data_array = VN.numpy_to_vtk(currentBedElev)
            temp_point_data_array.SetName(bedElevationVarName)
            point_data.AddArray(temp_point_data_array)

            # write to vtk file
            unstr_writer = vtk.vtkUnstructuredGridWriter()     #ASCII
            #unstr_writer = vtk.vtkXMLUnstructuredGridWriter()   #Binary
            unstr_writer.SetFileName(vtkFileName)
            unstr_writer.SetInputData(uGrid)
            unstr_writer.Write()

            #add the vtkFileName to vtkFileNameList
            vtkFileNameList.append(vtkFileName)

        #vtkFileNameList
        return vtkFileNameList

    
    def outputXMDFDataToPINNData(self, bNodal, bBoundary=False, dir=''):
        """Output XMDF result data to PINN data files
        - data_points.npy: rows of (x, y) coordinates of the points for steady state data (currently only support one time step)
        - data_values.npy: rows of (h, u, v) values of the solution variables at the points in data_points.npy
        - data_flags.npy: rows of flags for the points in data_points.npy (h_flag, u_flag, v_flag)
        - all_data_points_stats.pny: statistics of the data points (min, max, mean, std, etc.)

        This has to be called after the XMDF data have been loaded by calling readSRHXMDFFile(...).

        Parameters
        ----------
            bNodal : bool
                whether export nodal data or cell center data. Currently, it can't output both.
            bBoundary : bool
                whether export boundary data            
            dir : str, optional
                directory to write to

        Returns
        -------
            None

        """
        if gVerbose: print("Output all data in the XMDF file to PINN data ...")

        if (bNodal and (not self.xmdfAllData_Nodal)) or ((not bNodal) and (not self.xmdfAllData_Cell)):
            print("Empty XMDF data arrays. Call readSRHXMDFFile() function first. Exiting ...")
            sys.exit()

       
        #build result variable names
        resultVarNames = list(self.xmdfAllData_Nodal.keys()) if bNodal else list(self.xmdfAllData_Cell.keys())

        if bNodal:
            if gVerbose: print("All nodal solution variable names: ", resultVarNames[:-1])
        else:
            if gVerbose: print("All cell center solution variable names: ", resultVarNames)

        #add bed elevation at nodes regardless whether bNodal is True or False. Nodal elevation is more accurate representation
        #of the terrain because cell center elevation is averaged from nodal elevations.
        #get the units of the results
        if self.control_type == "SRHHydro":
            units = self.srhhydro_obj.srhhydro_content['OutputUnit']
        elif self.control_type == "SIF":
            units = self.srhsif_obj.srhsif_content['OutputUnit']
        else:
            raise NotImplementedError("TODO: implement the outputXMDFDataToVTK for SIF file")        

        if gVerbose: print("SRH-2D result units:", units)
        
        velocityVarName = ''
        if units == "SI":            
            velocityVarName = "Velocity_m_p_s"
        else:
            velocityVarName = "Velocity_ft_p_s"

        if gVerbose: print("Nodal velocity name: ", velocityVarName)

  
        #loop through each time step
        timeArray = self.xmdfTimeArray_Nodal if bNodal else self.xmdfTimeArray_Cell

        #get the last time step
        timeI = timeArray.shape[0] - 1

        #get the last time step data
        resultData =      np.zeros((self.srhgeom_obj.numOfNodes, len(resultVarNames)), dtype="float32") if bNodal \
                         else np.zeros((self.srhgeom_obj.numOfElements, len(resultVarNames)), dtype="float32")
  
        if gVerbose: print("timeI and time = ", timeI, timeArray[timeI])

        #numpy array for all solution variables at one time
        resultData =      np.zeros((self.srhgeom_obj.numOfNodes, len(resultVarNames)), dtype="float32") if bNodal \
                         else np.zeros((self.srhgeom_obj.numOfElements, len(resultVarNames)), dtype="float32")
        
        #loop through each solution variable 
        for varName, varI in zip(resultVarNames, range(len(resultVarNames))):
            if gVerbose: print("varName = ", varName)
            #get the values of current solution varialbe at current time
            resultData[:,varI] =      self.xmdfAllData_Nodal[varName][timeI,:] if bNodal \
                                     else self.xmdfAllData_Cell[varName][timeI,:]
        
        # Get coordinates (x,y) for all points
        if bNodal:
            data_points = self.srhgeom_obj.nodeCoordinates[:, :2]  # Only x,y coordinates
        else:
            data_points = self.srhgeom_obj.elementCenters[:, :2]  # Only x,y coordinates

        # Get solution variables (h, u, v)
        # Find indices for water elevation and velocity components
        water_depth_idx = -1
        vel_x_idx = -1
        vel_y_idx = -1
        
        for i, var_name in enumerate(resultVarNames):
            if "Water_Depth" in var_name:
                water_depth_idx = i
            elif "Vel_X" in var_name:
                vel_x_idx = i
            elif "Vel_Y" in var_name:
                vel_y_idx = i

        if water_depth_idx == -1 or vel_x_idx == -1 or vel_y_idx == -1:
            print("Error: Could not find all required variables (water depth and velocity components)")
            sys.exit()

        # Extract h, u, v values
        h_values = resultData[:, water_depth_idx]
        u_values = resultData[:, vel_x_idx]
        v_values = resultData[:, vel_y_idx]

        #number of points so far (before we add wall boundary points)
        num_points_excluding_walls = len(data_points)

        # if bBoundary is True, we include points on the boundary (currenlty only the wall boundaries are supported)
        if bBoundary:
            #for wall boundary points
            wall_node_IDs = []
            x_values_wall = []
            y_values_wall = []
            h_values_wall = []
            u_values_wall = []
            v_values_wall = []

            #loop through all boundaries.
            for nodeString in self.srhgeom_obj.nodeStringsDict:

                #not all NodeStrings are boundaries. Could be monitoring lines. Need to exclude them.
                if nodeString not in self.srhgeom_obj.bcDict.keys():
                    continue

                #we only include wall boundaries (for now)
                if 'wall' in self.srhgeom_obj.bcDict[nodeString].lower():
                    #list of nodes in current nodeString
                    nodeString_nodeList = self.srhgeom_obj.nodeStringsDict[nodeString]                    

                    #loop through all edges in current boundary
                    for i in range(len(nodeString_nodeList)-1):
                        nodeID_1 = nodeString_nodeList[i]
                        nodeID_2 = nodeString_nodeList[i + 1]

                        # Get coordinates of the two nodes
                        node1_coords = self.srhgeom_obj.nodeCoordinates[nodeID_1-1]   #SRH-2D node IDs are 1-based
                        node2_coords = self.srhgeom_obj.nodeCoordinates[nodeID_2-1]   #SRH-2D node IDs are 1-based

                        if nodeID_1 not in wall_node_IDs:
                            wall_node_IDs.append(nodeID_1)

                            x_values_wall.append(node1_coords[0])
                            y_values_wall.append(node1_coords[1])
                            h_values_wall.append(0.0)   #h does not matter for wall boundary because it is not used in training 
                            u_values_wall.append(0.0)   #no-slip condition
                            v_values_wall.append(0.0)   #no-slip condition

                        if nodeID_2 not in wall_node_IDs:
                            wall_node_IDs.append(nodeID_2)

                            x_values_wall.append(node2_coords[0])
                            y_values_wall.append(node2_coords[1])
                            h_values_wall.append(0.0)   #h does not matter for wall boundary because it is not used in training 
                            u_values_wall.append(0.0)   #no-slip condition
                            v_values_wall.append(0.0)   #no-slip condition

            print("Found %d wall boundary points" % len(wall_node_IDs))

            #number of wall points
            num_wall_points = len(wall_node_IDs)

            #append the wall boundary points to the data_points
            if len(x_values_wall) > 0:  
                # First stack x and y values into a 2D array
                wall_points = np.column_stack((np.array(x_values_wall), np.array(y_values_wall)))
                # Then concatenate with the original data_points
                data_points = np.concatenate((data_points, wall_points))

                h_values = np.concatenate((h_values, np.array(h_values_wall)))
                u_values = np.concatenate((u_values, np.array(u_values_wall)))
                v_values = np.concatenate((v_values, np.array(v_values_wall)))                
                        
            
        # Compute U magnitude
        Umag = np.sqrt(u_values**2 + v_values**2)

        # Combine into data_values array
        data_values = np.column_stack((h_values, u_values, v_values))

        # Create flags array (1 for valid data, 0 for invalid)
        # For now, we'll mark all points as valid (1)
        data_flags = np.ones_like(data_values)

        # the data_flags for h on walls should be set to 0. 
        if bBoundary:
            data_flags[num_points_excluding_walls:num_points_excluding_walls + num_wall_points, 0] = 0

        #compute the statistics of the data points
        #min, max, mean, std, median, etc.
        x_min = np.min(data_points[:, 0])  #stats of x, y, and t are computed here, but not used in training. The normalization should use the statistics of the pde/boundary points from mesh_points.json.
        x_max = np.max(data_points[:, 0])
        x_mean = np.mean(data_points[:, 0])
        x_std = np.std(data_points[:, 0])
        y_min = np.min(data_points[:, 1])
        y_max = np.max(data_points[:, 1])
        y_mean = np.mean(data_points[:, 1])
        y_std = np.std(data_points[:, 1])
        t_min = 0.0
        t_max = 0.0
        t_mean = 0.0
        t_std = 0.0
        h_min = np.min(data_values[:, 0])
        h_max = np.max(data_values[:, 0])
        h_mean = np.mean(data_values[:, 0])
        h_std = np.std(data_values[:, 0])
        u_min = np.min(data_values[:, 1])
        u_max = np.max(data_values[:, 1])
        u_mean = np.mean(data_values[:, 1])
        u_std = np.std(data_values[:, 1])
        v_min = np.min(data_values[:, 2])
        v_max = np.max(data_values[:, 2])
        v_mean = np.mean(data_values[:, 2])
        v_std = np.std(data_values[:, 2])
        Umag_min = np.min(Umag)
        Umag_max = np.max(Umag)
        Umag_mean = np.mean(Umag)
        Umag_std = np.std(Umag)

        all_data_points_stats = np.array([x_min, x_max, x_mean, x_std, y_min, y_max, y_mean, y_std, t_min, t_max, t_mean, t_std, h_min, h_max, h_mean, h_std, u_min, u_max, u_mean, u_std, v_min, v_max, v_mean, v_std, Umag_min, Umag_max, Umag_mean, Umag_std])
        print(f"x_min: {x_min}, x_max: {x_max}, x_mean: {x_mean}, x_std: {x_std}")
        print(f"y_min: {y_min}, y_max: {y_max}, y_mean: {y_mean}, y_std: {y_std}")
        print(f"t_min: {t_min}, t_max: {t_max}, t_mean: {t_mean}, t_std: {t_std}")
        print(f"h_min: {h_min}, h_max: {h_max}, h_mean: {h_mean}, h_std: {h_std}")
        print(f"u_min: {u_min}, u_max: {u_max}, u_mean: {u_mean}, u_std: {u_std}")
        print(f"v_min: {v_min}, v_max: {v_max}, v_mean: {v_mean}, v_std: {v_std}")
        print(f"Umag_min: {Umag_min}, Umag_max: {Umag_max}, Umag_mean: {Umag_mean}, Umag_std: {Umag_std}")    

        # Save files
        if dir:
            os.makedirs(dir, exist_ok=True)
            prefix = os.path.join(dir, '/')
        else:
            os.makedirs('pinn_points', exist_ok=True)
            prefix = 'pinn_points/'

        np.save(f'{prefix}data_points.npy', data_points)
        np.save(f'{prefix}data_values.npy', data_values)
        np.save(f'{prefix}data_flags.npy', data_flags)
        np.save(f'{prefix}all_data_points_stats.npy', all_data_points_stats)

        print(f"Saved data points shape: {data_points.shape}")
        print(f"Saved data values shape: {data_values.shape}")
        print(f"Saved data flags shape: {data_flags.shape}")
        print(f"Files saved in: {dir}")

        #print the last 10 data points
        print(f"Last 10 data points: {data_points[-10:]}")
        print(f"Last 10 data values: {data_values[-10:]}")
        print(f"Last 10 data flags: {data_flags[-10:]}")

        #save the data_points, data_values, and data_flags to a vtk file (for visualization)
        vtkFileName = os.path.join(dir, 'data_points.vtk')
        
        # Create points array with z=0 for 2D points
        points = np.column_stack((data_points, np.zeros(len(data_points))))
        
        # Create point data dictionary
        point_data = {
            'h': data_values[:, 0],
            'u': data_values[:, 1],
            'v': data_values[:, 2],
            'h_flag': data_flags[:, 0],
            'u_flag': data_flags[:, 1],
            'v_flag': data_flags[:, 2]
        }
        
        # Create meshio mesh object with point cloud
        mesh = meshio.Mesh(
            points=points,
            cells=[("vertex", np.arange(len(points)).reshape(-1, 1))],  # Create vertex cells for each point
            point_data=point_data
        )
        
        # Write to VTK file
        mesh.write(vtkFileName)
        
        if gVerbose:
            print(f"Saved visualization data to: {vtkFileName}")

    def readSRHCFiles(self, case_name):
        """Read SRH-2D results from SRHC files (cell-centered data)
    
        Args:
            case_name (str): Base filename without _SIF.dat extension
        """
        # Get directory and base name
        directory = os.path.dirname(self.srhsif_obj.srhsif_filename)
        
        # Find all SRHC files
        pattern = os.path.join(directory, f"{case_name}_SRHC*.dat")
        srhc_files = sorted(glob.glob(pattern), key=self.natural_sort_key)

        #print("srhc_files = ", srhc_files)

        if gVerbose:
            print(f"Found {len(srhc_files)} SRHC files matching pattern: {pattern}")
            print(f"SRHC files: {srhc_files}")
        
        if not srhc_files:
            print(f"No SRHC files found matching pattern: {pattern}")
            return
            
        # Initialize data structures for cell-centered data
        self.xmdfTimeArray_Cell = []
        self.xmdfAllData_Cell = {}
        
        # Read first file to get variable names
        data = np.genfromtxt(srhc_files[0], delimiter=',', names=True, dtype=None, encoding='utf-8')
        # Remove the empty field name caused by trailing comma
        var_names = list(data.dtype.names)[:-1]  # Skip the last empty field
            
        # Initialize data arrays for each variable
        for var_name in var_names[3:]:  # 3: is for skipping Point_ID, X, and Y
            self.xmdfAllData_Cell[var_name] = []
            
        # Read each SRHC file (each represents a timestep)
        for timestep, srhc_file in enumerate(srhc_files):
            # Use genfromtxt to read the data
            data = np.genfromtxt(srhc_file, delimiter=',', skip_header=1, filling_values=0.0)
            # Remove the last column (empty due to trailing comma)
            data = data[:, :-1]
            
            # Store time (we'll use timestep number since actual time isn't in SRHC files)
            self.xmdfTimeArray_Cell.append(float(timestep))
            
            # Store data for each variable
            for i, var_name in enumerate(var_names[3:], 3):  # Skip Point_ID, X, and Y
                values = data[:, i]     

                #debug
                if gVerbose:
                    print("var_name = ", var_name)
                    print("values = ", values[0:10])

                # Fix water elevation -999 values
                if var_name == "Water_Elev_ft" or var_name == "Water_Elev_m":
                    values = np.where(values == -999, 
                                    self.srhgeom_obj.elementBedElevation,
                                    values)
                    
                self.xmdfAllData_Cell[var_name].append(values)
                
        # Convert lists to numpy arrays
        self.xmdfTimeArray_Cell = np.array(self.xmdfTimeArray_Cell)
        for var_name in self.xmdfAllData_Cell:
            self.xmdfAllData_Cell[var_name] = np.array(self.xmdfAllData_Cell[var_name])

    def outputSRHCDataToVTK(self, timeStep=-1, lastTimeStep=False, dir=''):
        """Output SRHC data to VTK format
        
        Args:
            timeStep (int): Specific timestep to output (-1 for all timesteps)
            lastTimeStep (bool): Only output last timestep if True
            dir (str): Output directory
        """
        # SRHC data is always at cell center
        bNodal = False    

        # Use existing outputXMDFDataToVTK function since data structure is same
        return self.outputXMDFDataToVTK(bNodal, timeStep, lastTimeStep, dir)


    def outputVTK(self, vtkFileName, resultVarNames, resultData, bNodal):
        """ Output result to VTK file

        The supplied resultVarNames and resultData should be compatible with the mesh. If resultVarNames is empty,
        it only outputs the mesh with no data.

        Parameters
        ----------
        vtkFileName : str
            name for the output vtk file
        resultVarNames : str
            result variable names
        resultData : `numpy.ndarray`
            2D array containing result data
        bNodal : bool
            whether the data is nodal (True) or at cell center (False)


        Returns
        -------

        """

        if gVerbose: print("Output to VTK ...")

        try:
            fid = open(vtkFileName, 'w')
        except IOError:
            print('vtk file open error')
            sys.exit()

        #get the units of the results
        if self.control_type == "SRHHydro":
            units = self.srhhydro_obj.srhhydro_content['OutputUnit']
        elif self.control_type == "SIF":
            units = self.srhsif_obj.srhsif_content['OutputUnit']
        else:
            raise NotImplementedError("TODO: implement the outputVTK for SIF file")


        if gVerbose: print("SRH-2D result units:", units)


        fid.write('# vtk DataFile Version 3.0\n')
        fid.write('Results from SRH-2D Modeling Run\n')
        fid.write('ASCII\n')
        fid.write('DATASET UNSTRUCTURED_GRID\n')
        fid.write('\n')

        # output points
        fid.write('POINTS %d double\n' % self.srhgeom_obj.nodeCoordinates.shape[0])

        point_id = 0  # point ID counter
        for k in range(self.srhgeom_obj.nodeCoordinates.shape[0]):
            point_id += 1
            fid.write(" ".join(map(str, self.srhgeom_obj.nodeCoordinates[k])))
            fid.write("\n")

        # output elements
        fid.write('CELLS %d %d \n' % (self.srhgeom_obj.elementNodesList.shape[0],
                                      self.srhgeom_obj.elementNodesList.shape[0] + np.sum(self.srhgeom_obj.elementNodesCount)))

        cell_id = 0  # cell ID counter
        for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
            cell_id += 1
            fid.write('%d ' % self.srhgeom_obj.elementNodesCount[k])
            fid.write(" ".join(map(str, self.srhgeom_obj.elementNodesList[k][:self.srhgeom_obj.elementNodesCount[k]] - 1)))
            fid.write("\n")

        # output element types
        fid.write('CELL_TYPES %d \n' % self.srhgeom_obj.elementNodesList.shape[0])

        for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
            fid.write('%d ' % self.srhgeom_obj.vtkCellTypeCode[k])
            if (((k + 1) % 20) == 0):
                fid.write('\n')

        fid.write('\n')

        # How many solution data entries to output depending on whether it is nodal or not
        entryCount = -1
        if bNodal:
            entryCount = self.srhgeom_obj.nodeCoordinates.shape[0]
        else:
            entryCount = self.srhgeom_obj.elementNodesList.shape[0]

        # output solution variables: only if there is solution variable
        if len(resultVarNames) != 0:
            if not bNodal:
                if gVerbose: print('Output variables are at cell centers. \n')
                fid.write('CELL_DATA %d\n' % self.srhgeom_obj.elementNodesList.shape[0])
            else:
                if gVerbose: print('Output variables are at vertices. \n')
                fid.write('POINT_DATA %d\n' % self.srhgeom_obj.nodeCoordinates.shape[0])

            # column numbers for Vel_X and Vel_Y for vector assemble
            nColVel_X = -1
            nColVel_Y = -1

            # First output all solution variables as scalars
            #print('The following solution variables are processed: \n')
            for k in range(len(resultVarNames)):
                #print('     %s\n' % resultVarNames[k])

                if resultVarNames[k].find('Vel_X') != -1:
                    nColVel_X = k
                elif resultVarNames[k].find('Vel_Y') != -1:
                    nColVel_Y = k

                fid.write('SCALARS %s double 1 \n' % resultVarNames[k])
                fid.write('LOOKUP_TABLE default\n')

                for entryI in range(entryCount):
                    fid.write('%f ' % resultData[entryI][k])

                    if (((entryI + 1) % 20) == 0):
                        fid.write('\n')

                fid.write('\n \n')

            # Then output Vel_X and Vel_Y as velocity vector (Vel_Z = 0.0)
            # print('nColVel_X, nColVel_Y = %d, %d' % (nColVel_X, nColVel_Y))
            if (nColVel_X != -1) and (nColVel_Y != -1):
                if units == "SI":
                    fid.write('VECTORS Velocity_m_p_s double \n')
                else:
                    fid.write('VECTORS Velocity_ft_p_s double \n')

                for entryI in range(entryCount):
                    fid.write('%f %f 0.0   ' % (resultData[entryI][nColVel_X],
                                                resultData[entryI][nColVel_Y]))

                    if (((entryI + 1) % 20) == 0):
                        fid.write('\n')

                fid.write('\n \n')

        fid.close()

    def output_2d_mesh_to_vtk(self, meshVTKFileName, bFlat=False, dir=''):
        """ output the flat mesh to vtk

        Parameters
        ----------
        meshVTKFileName : str
            file name for the mesh
        bFlat : bool
            whether to make the 2d mesh flat (node's z coordinate -> 0)
        dir : str, optional
            directory to write to

        Returns
        -------

        """
        #just call srhgeom_obj's function
        self.srhgeom_obj.output_2d_mesh_to_vtk(meshVTKFileName, bFlat, dir)

    def readOneSRHFile(self, srhFileName):
        """ Read one single SRH-2D result file in SRHC (cell center) or SRH (point) format.

        Note: SRH-2D outputs an extra "," to each line. As a result, Numpy's
        genfromtext(...) function adds a column of "nan" to the end. Need to take care of this.

        Parameters
        ----------
        srhFileName : str
            file name for the SRH result file

        Returns
        -------


        """

        if gVerbose: print("Reading the SRH/SRHC result file ...")

        data = np.genfromtxt(srhFileName, delimiter=',', names=True)

        return data.dtype.names[:-1], data

    def readTECFile(self, tecFileName):
        """Read SRH-2D results in Tecplot format

        Parameters
        ----------
        tecFileName : str
            Tecplot file name

        Returns
        -------

        """

        readerTEC = vtk.vtkTecplotReader()
        readerTEC.SetFileName(tecFileName)
        # 'update' the reader i.e. read the Tecplot data file
        readerTEC.Update()

        polydata = readerTEC.GetOutput()

        # print(polydata)  #this is a multibock dataset
        # print(polydata.GetBlock(0))  #we only need the first block

        # If there are no points in 'vtkPolyData' something went wrong
        if polydata.GetBlock(0).GetNumberOfPoints() == 0:
            raise ValueError(
                "No point data could be loaded from '" + tecFileName)
            return None

        #return of the first block (there should be only one block)
        return polydata.GetBlock(0)


    def cell_center_to_vertex(self, cell_data, vertex_data):
        """Interpolate result from cell center to vertex

        Not implemented yet.

        Returns
        -------

        """
        pass

    def vertex_to_cell_center(self, vertex_data, cell_data):
        """Interpolate result from vertex to cell center

        Not implemented yet.

        Returns
        -------

        """
        pass

    def natural_sort_key(self, s):
        return [int(text) if text.isdigit() else text.lower()
                for text in re.split(r'(\d+)', s)]
