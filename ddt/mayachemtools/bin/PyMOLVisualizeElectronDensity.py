#!/bin/env python
#
# File: PyMOLVisualizeElectronDensity.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
#
# The functionality available in this script is implemented using PyMOL, a
# molecular visualization system on an open source foundation originally
# developed by Warren DeLano.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

from __future__ import print_function

# Add local python path to the global path and import standard library modules...
import os
import sys;  sys.path.insert(0, os.path.join(os.path.dirname(sys.argv[0]), "..", "lib", "Python"))
import time
import re

# PyMOL imports...
try:
    import pymol
    # Finish launching PyMOL in  a command line mode for batch processing (-c)
    # along with the following options:  disable loading of pymolrc and plugins (-k);
    # suppress start up messages (-q)
    pymol.finish_launching(['pymol', '-ckq'])
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import PyMOL module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your PyMOL environment and try again.\n\n")
    sys.exit(1)

# MayaChemTools imports...
try:
    from docopt import docopt
    import MiscUtil
    import PyMOLUtil
except ImportError as ErrMsg:
    sys.stderr.write("\nFailed to import MayaChemTools module/package: %s\n" % ErrMsg)
    sys.stderr.write("Check/update your MayaChemTools environment and try again.\n\n")
    sys.exit(1)

ScriptName = os.path.basename(sys.argv[0])
Options = {}
OptionsInfo = {}

def main():
    """Start execution of the script"""
    
    MiscUtil.PrintInfo("\n%s (PyMOL v%s; %s) Starting...\n" % (ScriptName, pymol.cmd.get_version()[1], time.asctime()))
    
    (WallClockTime, ProcessorTime) = MiscUtil.GetWallClockAndProcessorTime()
    
    # Retrieve command line arguments and options...
    RetrieveOptions()
    
    # Process and validate command line arguments and options...
    ProcessOptions()

    # Perform actions required by the script...
    GenerateElectronDensityVisualization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateElectronDensityVisualization():
    """Generate electron density visualization."""

    Outfile = OptionsInfo["PMLOutfile"]
    OutFH = open(Outfile, "w")
    if OutFH is None:
        MiscUtil.PrintError("Failed to open output fie %s " % Outfile)
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % Outfile)

    # Setup PyMOL object names...
    PyMOLObjectNames = SetupPyMOLObjectNames()

    # Setup header and complex view...
    WritePMLHeader(OutFH, ScriptName)
    WritePyMOLParameters(OutFH)
    WriteComplexView(OutFH, PyMOLObjectNames)

    # Setup chain views...
    FirstChain = True
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        WriteChainView(OutFH, PyMOLObjectNames, ChainID)
        
        # Setup ligand views...
        FirstLigand = True
        for LigandID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]:
            WriteChainLigandView(OutFH, PyMOLObjectNames, ChainID, LigandID)
            
            # Set up ligand level group...
            Enable, Action = [False, "close"]
            if FirstLigand:
                FirstLigand = False
                Enable, Action = [True, "open"]
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroup"], PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"], Enable, Action)
            
        # Setup Chain level group...
        Enable, Action = [False, "close"]
        if FirstChain:
            FirstChain = False
            Enable, Action = [True, "open"]
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"], Enable, Action)
        
    # Set up complex level group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["PDBGroup"], PyMOLObjectNames["PDBGroupMembers"], True, "open")

    # Delete empty PyMOL objects...
    DeleteEmptyPyMOLObjects(OutFH, PyMOLObjectNames)
        
    OutFH.write("""\ncmd.orient("visible", animate = -1)\n""")
    
    OutFH.close()

    # Generate PSE file as needed...
    if OptionsInfo["PSEOut"]:
        GeneratePyMOLSessionFile()

def WritePMLHeader(OutFH, ScriptName):
    """Write out PML setting up complex view"""

    HeaderInfo = PyMOLUtil.SetupPMLHeaderInfo(ScriptName)
    OutFH.write("%s\n" % HeaderInfo)

def WritePyMOLParameters(OutFH):
    """Write out PyMOL global parameters. """

    PMLCmds = []
    PMLCmds.append("""cmd.set("mesh_width", %.2f)""" % (OptionsInfo["MeshWidth"]))
    PMLCmds.append("""cmd.set("transparency", %.2f, "", 0)""" % (OptionsInfo["SurfaceTransparency"]))
    PMLCmds.append("""cmd.set("label_font_id", %s)""" % (OptionsInfo["LabelFontID"]))
    PML = "\n".join(PMLCmds)
    
    OutFH.write("""\n""\n"Setting up PyMOL gobal parameters..."\n""\n""")
    OutFH.write("%s\n" % PML)
    
def WriteComplexView(OutFH, PyMOLObjectNames):
    """Write out PML for viewing polymer complex along with electron density."""

    # Setup complex...
    PML = PyMOLUtil.SetupPMLForPolymerComplexView(PyMOLObjectNames["Complex"], OptionsInfo["Infile"], True)
    OutFH.write("""\n""\n"Loading %s and setting up view for complex..."\n""\n""" % OptionsInfo["Infile"])
    OutFH.write("%s\n" % PML)

    # Setup electron density maps and meshes...
    WriteComplexCompositeElectronDensityMapView(OutFH, PyMOLObjectNames)
    if OptionsInfo["DifferenceMapAvailable"]:
        WriteComplexDifferenceElectronDensityMapView(OutFH, PyMOLObjectNames)

    # Setup complex group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexGroup"], PyMOLObjectNames["ComplexGroupMembers"], False, "close")

def WriteComplexCompositeElectronDensityMapView(OutFH, PyMOLObjectNames):
    """Write out PML for viewing composite maps."""

    # Load composite map (2Fo - Fc) and setup mesh views...
    MapFileName = OptionsInfo["EDMapFiles"]["CompositeMap"]
    Info = """\
""
"Loading composite map (2Fo - Fc) %s and setting up mesh view for complex..."
"" """ % MapFileName
    OutFH.write("\n%s\n" % Info)

    MapName = PyMOLObjectNames["ComplexCompositeEDMap"]
    ComplexName = PyMOLObjectNames["Complex"]
    
    ContourLevel = OptionsInfo["MeshLevelCompositeMap"]
    Color = OptionsInfo["MeshColorCompositeMap"]
    
    MeshName = PyMOLObjectNames["ComplexCompositeEDMesh"]
    SurfaceName = PyMOLObjectNames["ComplexCompositeEDSurface"]
    
    PML = SetupPMLForElectronDensityMap(MapFileName, MapName, True)
    OutFH.write("%s\n" % PML)

    EnableSurface = False if OptionsInfo["MeshComplex"] else True
    if OptionsInfo["MeshComplex"]:
        PML = SetupPMLForElectronDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = True, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    if OptionsInfo["SurfaceComplex"]:
        PML = SetupPMLForElectronDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurface, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexCompositeEDGroup"], PyMOLObjectNames["ComplexCompositeEDGroupMembers"], True, "close")
    
def WriteComplexDifferenceElectronDensityMapView(OutFH, PyMOLObjectNames):
    """Write out PML for viewing difference maps."""

    # Load difference map (Fo - Fc ) and setup mesh views...
    MapFileName = OptionsInfo["EDMapFiles"]["DifferenceMap"]
    Info = """\
""
"Loading difference map (Fo - Fc ) map %s and setting up mesh views..."
"" """ % MapFileName
    OutFH.write("\n%s\n" % Info)

    MapName = PyMOLObjectNames["ComplexDiffEDMap"]
    ComplexName = PyMOLObjectNames["Complex"]

    ContourLevel1 = OptionsInfo["Mesh1LevelDiffMap"]
    ContourLevel2 = OptionsInfo["Mesh2LevelDiffMap"]
    Color1 = OptionsInfo["Mesh1ColorDiffMap"]
    Color2 = OptionsInfo["Mesh2ColorDiffMap"]
    
    Mesh1Name = PyMOLObjectNames["ComplexDiffEDMesh1"]
    Surface1Name = PyMOLObjectNames["ComplexDiffEDSurface1"]
    Mesh2Name = PyMOLObjectNames["ComplexDiffEDMesh2"]
    Surface2Name = PyMOLObjectNames["ComplexDiffEDSurface2"]
    
    EnableSurface = False if OptionsInfo["MeshComplex"] else True
    
    PML = SetupPMLForElectronDensityMap(MapFileName, MapName, True)
    OutFH.write("%s\n" % PML)
    
    if OptionsInfo["MeshComplex"]:
        PML = SetupPMLForElectronDensityMesh(MapName, Mesh1Name, ContourLevel1, Color1, Enable = True, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)
        
    if OptionsInfo["SurfaceComplex"]:
        PML = SetupPMLForElectronDensitySurface(MapName, Surface1Name, ContourLevel1, Color1, Enable = EnableSurface, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)
    
    if OptionsInfo["MeshComplex"]:
        PML = SetupPMLForElectronDensityMesh(MapName, Mesh2Name, ContourLevel2, Color2, Enable = True, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)
    
    if OptionsInfo["SurfaceComplex"]:
        PML = SetupPMLForElectronDensitySurface(MapName, Surface2Name, ContourLevel2, Color2, Enable = EnableSurface, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)
    
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexDiffEDGroup"], PyMOLObjectNames["ComplexDiffEDGroupMembers"], True, "close")

def WriteChainView(OutFH, PyMOLObjectNames, ChainID):
    """Write out PML for viewing chain."""
    
    OutFH.write("""\n""\n"Setting up views for chain %s..."\n""\n""" % ChainID)
    
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    
    # Setup chain complex group view...
    WriteChainComplexAndMeshViews(OutFH, PyMOLObjectNames, ChainID)

    # Setup chain view...
    WriteChainAloneViews(OutFH, PyMOLObjectNames, ChainID)
    
    # Setup chain solvent view...
    PML = PyMOLUtil.SetupPMLForSolventView(PyMOLObjectNames["Chains"][ChainID]["Solvent"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

    # Setup chain inorganic view...
    PML = PyMOLUtil.SetupPMLForInorganicView(PyMOLObjectNames["Chains"][ChainID]["Inorganic"], ChainComplexName, False)
    OutFH.write("\n%s\n" % PML)

def WriteChainComplexAndMeshViews(OutFH, PyMOLObjectNames, ChainID):
    """Write chain complex and mesh views. """
    
    # Setup chain complex...
    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    PML = PyMOLUtil.SetupPMLForPolymerChainComplexView(ChainComplexName, PyMOLObjectNames["Complex"], ChainID, True)
    OutFH.write("%s\n" % PML)

    MeshChainComplex = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["MeshChainComplex"][ChainID]
    SurfaceChainComplex = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["SurfaceChainComplex"][ChainID]
    
    EnableSurface = False if MeshChainComplex else True
    
    if MeshChainComplex or SurfaceChainComplex:
        # Set up composite mesh and group...
        MapName = PyMOLObjectNames["ComplexCompositeEDMap"]
        ContourLevel = OptionsInfo["MeshLevelCompositeMap"]
        Color = OptionsInfo["MeshColorCompositeMap"]
        
        MeshName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDMesh"]
        SurfaceName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDSurface"]
        
        if MeshChainComplex:
            PML = SetupPMLForElectronDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = True, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        if SurfaceChainComplex:
            PML = SetupPMLForElectronDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurface, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroupMembers"], True, "close")
        
        if OptionsInfo["DifferenceMapAvailable"]:
            # Set up difference meshes and group...
            MapName = PyMOLObjectNames["ComplexDiffEDMap"]
            
            ContourLevel1 = OptionsInfo["Mesh1LevelDiffMap"]
            ContourLevel2 = OptionsInfo["Mesh2LevelDiffMap"]
            Color1 = OptionsInfo["Mesh1ColorDiffMap"]
            Color2 = OptionsInfo["Mesh2ColorDiffMap"]
        
            Mesh1Name = PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDMesh1"]
            Surface1Name = PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDSurface1"]
            Mesh2Name = PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDMesh2"]
            Surface2Name = PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDSurface2"]
        
            if MeshChainComplex:
                PML = SetupPMLForElectronDensityMesh(MapName, Mesh1Name, ContourLevel1, Color1, Enable = True, Selection = ChainComplexName)
                OutFH.write("\n%s\n" % PML)
            
            if SurfaceChainComplex:
                PML = SetupPMLForElectronDensitySurface(MapName, Surface1Name, ContourLevel1, Color1, Enable = EnableSurface, Selection = ChainComplexName)
                OutFH.write("\n%s\n" % PML)
            
            if MeshChainComplex:
                PML = SetupPMLForElectronDensityMesh(MapName, Mesh2Name, ContourLevel2, Color2, Enable = True, Selection = ChainComplexName)
                OutFH.write("\n%s\n" % PML)
            
            if SurfaceChainComplex:
                PML = SetupPMLForElectronDensitySurface(MapName, Surface2Name, ContourLevel2, Color2, Enable = EnableSurface, Selection = ChainComplexName)
                OutFH.write("\n%s\n" % PML)
        
            GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"], True, "close")
    
    # Setup chain complex group...
    EnableChainComplexGroup = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainComplexGroup"][ChainID]
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"], EnableChainComplexGroup, "close")
    
def WriteChainAloneViews(OutFH, PyMOLObjectNames, ChainID):
    """Write individual chain views. """

    ChainComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
    
    # Setup chain view...
    ChainName = PyMOLObjectNames["Chains"][ChainID]["ChainAlone"]
    PML = PyMOLUtil.SetupPMLForPolymerChainView(ChainName, ChainComplexName, Enable = True)
    OutFH.write("\n%s\n" % PML)

    # Setup chain putty by B-factor view...
    if OptionsInfo["BFactorChainCartoonPutty"]:
        BFactorPuttyName = PyMOLObjectNames["Chains"][ChainID]["ChainAloneBFactorPutty"]
        PML = PyMOLUtil.SetupPMLForBFactorPuttyView(BFactorPuttyName, ChainName, ColorPalette = OptionsInfo["BFactorColorPalette"], Enable = False)
        OutFH.write("\n%s\n" % PML)
        
    # Setup chain group...
    EnableChainAloneGroup = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainAloneGroup"][ChainID]
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"], EnableChainAloneGroup, "close")
    
def WriteChainLigandView(OutFH, PyMOLObjectNames, ChainID, LigandID):
    """Write out PML for viewing ligand in a chain."""
    
    for GroupID in ["Ligand", "Pocket", "PocketSolvent", "PocketInorganic"]:
        ComplexName = PyMOLObjectNames["Chains"][ChainID]["ChainComplex"]
        LigandName = PyMOLObjectNames["Ligands"][ChainID][LigandID]["Ligand"]
        
        # Setup main object...
        GroupTypeObjectID = "%s" % (GroupID)
        GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
        
        if re.match("^Ligand$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandView(GroupTypeObjectName, ComplexName, LigandID, True)
            OutFH.write("%s\n" % PML)
        elif re.match("^Pocket$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for pocket around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], True)
            OutFH.write("%s\n" % PML)
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], GroupTypeObjectName))
        elif re.match("^PocketSolvent$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for solvent in pockect around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketSolventView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], Enable = True)
            OutFH.write("%s\n" % PML)
        elif re.match("^PocketInorganic$", GroupID, re.I):
            OutFH.write("""\n""\n"Setting up views for inorganic in pockect around ligand %s in chain %s..."\n""\n""" % (LigandID, ChainID))
            PML = PyMOLUtil.SetupPMLForLigandPocketInorganicView(GroupTypeObjectName, ComplexName, LigandName, OptionsInfo["PocketDistanceCutoff"], Enable = True)
            OutFH.write("%s\n" % PML)
        
        # Set up composite mesh and group...
        CompositeMeshGroupID = "%sCompositeEDMeshGroup" % (GroupID)
        CompositeMeshGroupMembersID = "%sCompositeEDMeshGroupMembers" % (GroupID)
        CompositeMeshID = "%sCompositeEDMesh" % (GroupID)
        CompositeSurfaceID = "%sCompositeEDSurface" % (GroupID)

        CompositeMapName = PyMOLObjectNames["ComplexCompositeEDMap"]
        CompositeMeshName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshID]
        CompositeSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeSurfaceID]
        CompositeMeshGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupID]
        CompositeMeshGroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupMembersID]
        
        PML = SetupPMLForElectronDensityMesh(CompositeMapName, CompositeMeshName, OptionsInfo["MeshLevelCompositeMap"], OptionsInfo["MeshColorCompositeMap"], Enable = True, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        PML = SetupPMLForElectronDensitySurface(CompositeMapName, CompositeSurfaceName, OptionsInfo["MeshLevelCompositeMap"], OptionsInfo["MeshColorCompositeMap"], Enable = False, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, CompositeMeshGroupName, CompositeMeshGroupMembers, True, "close")
        
        if OptionsInfo["DifferenceMapAvailable"]:
            # Set up difference meshes and group...
            DiffMeshGroupID = "%sDiffEDMeshGroup" % (GroupID)
            DiffMeshGroupMembersID = "%sDiffEDMeshGroupMembers" % (GroupID)
            DiffMesh1ID = "%sDiffEDMesh1" % (GroupID)
            DiffSurface1ID = "%sDiffEDSurface1" % (GroupID)
            DiffMesh2ID = "%sDiffEDMesh2" % (GroupID)
            DiffSurface2ID = "%sDiffEDSurface2" % (GroupID)
        
            DiffMapName = PyMOLObjectNames["ComplexDiffEDMap"]
            DiffMesh1Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMesh1ID]
            DiffSurface1Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffSurface1ID]
            DiffMesh2Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMesh2ID]
            DiffSurface2Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffSurface2ID]
            DiffMeshGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMeshGroupID]
            DiffMeshGroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMeshGroupMembersID]
            
            PML = SetupPMLForElectronDensityMesh(DiffMapName, DiffMesh1Name, OptionsInfo["Mesh1LevelDiffMap"], OptionsInfo["Mesh1ColorDiffMap"],  Enable = True, Selection = GroupTypeObjectName)
            OutFH.write("\n%s\n" % PML)
            PML = SetupPMLForElectronDensitySurface(DiffMapName, DiffSurface1Name, OptionsInfo["Mesh1LevelDiffMap"], OptionsInfo["Mesh1ColorDiffMap"],  Enable = False, Selection = GroupTypeObjectName)
            OutFH.write("\n%s\n" % PML)
            
            PML = SetupPMLForElectronDensityMesh(DiffMapName, DiffMesh2Name, OptionsInfo["Mesh2LevelDiffMap"], OptionsInfo["Mesh2ColorDiffMap"],  Enable = True, Selection = GroupTypeObjectName)
            OutFH.write("\n%s\n" % PML)
            PML = SetupPMLForElectronDensitySurface(DiffMapName, DiffSurface2Name, OptionsInfo["Mesh2LevelDiffMap"], OptionsInfo["Mesh2ColorDiffMap"],  Enable = False, Selection = GroupTypeObjectName)
            OutFH.write("\n%s\n" % PML)
        
            GenerateAndWritePMLForGroup(OutFH, DiffMeshGroupName, DiffMeshGroupMembers, True, "close")
        
        # Set up polar contacts...
        if re.match("^(Pocket|PocketSolvent|PocketInorganic)$", GroupID, re.I):
            PolarContactsID = "%sPolarContacts" % (GroupID)
            PolarContactsName = PyMOLObjectNames["Ligands"][ChainID][LigandID][PolarContactsID]
            
            PolarContactsColor = OptionsInfo["PocketContactsLigandColor"]
            if re.match("^PocketSolvent$", GroupID, re.I):
                PolarContactsColor = OptionsInfo["PocketContactsSolventColor"]
            elif re.match("^PocketInorganic$", GroupID, re.I):
                PolarContactsColor = OptionsInfo["PocketContactsInorganicColor"]
            
            PML = PyMOLUtil.SetupPMLForPolarContactsView(PolarContactsName, LigandName, GroupTypeObjectName, Enable = False, Color = PolarContactsColor)
            OutFH.write("\n%s\n" % PML)
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (PolarContactsColor, PolarContactsName))
            
        # Set up hydrophobic surface...
        if re.match("^Pocket$", GroupID, re.I) and OptionsInfo["PocketSurface"]:
            HydrophobicSurfaceID = "%sHydrophobicSurface" % (GroupID)
            HydrophobicSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID]
            PML = PyMOLUtil.SetupPMLForHydrophobicSurfaceView(HydrophobicSurfaceName, GroupTypeObjectName, ColorPalette = "RedToWhite", Enable = False)
            OutFH.write("\n%s\n" % PML)
            
            OutFH.write("""cmd.set("label_color", "%s", "%s")\n""" % (OptionsInfo["PocketLabelColor"], HydrophobicSurfaceName))
        
        # Setup group....
        GroupNameID = "%sGroup" % (GroupID)
        GroupMembersID = "%sGroupMembers" % (GroupID)
        GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]
        GroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID]

        Action = "close"
        Enable = False
        if  re.match("^(Ligand|Pocket)$", GroupID, re.I):
            Action = "open"
            Enable = True
        GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable, Action)

def GenerateAndWritePMLForGroup(OutFH, GroupName, GroupMembers, Enable = False, Action = "close"):
    """Generate and write PML for group. """
    
    PML = PyMOLUtil.SetupPMLForGroup(GroupName, GroupMembers, Enable, Action)
    OutFH.write("""\n""\n"Setting up group %s..."\n""\n""" % GroupName)
    OutFH.write("%s\n" % PML)

def SetupPMLForElectronDensityMap(MapFileName, MapName, Enable = True):
    """Setup PML for loading and viewing electron density map. """

    PMLCmds = []
    PMLCmds.append("""cmd.load("%s", "%s")""" % (MapFileName, MapName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MapName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML
    
def SetupPMLForElectronDensityMesh(MapName, MeshName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for electron density mesh. """

    Carve = OptionsInfo["MeshCarveRadius"]
    
    PMLCmds = []
    if Selection is None:
        PMLCmds.append("""cmd.isomesh("%s", "%s", %.1f)""" % (MeshName, MapName, SigmaLevel))
    else:
        PMLCmds.append("""cmd.isomesh("%s", "%s", %.1f, "(%s)", carve = %.1f)""" % (MeshName, MapName, SigmaLevel, Selection, Carve))
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, MeshName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MeshName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForElectronDensitySurface(MapName, SurfaceName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for electron density surface. """

    Carve = OptionsInfo["MeshCarveRadius"]
    
    PMLCmds = []
    if Selection is None:
        PMLCmds.append("""cmd.isosurface("%s", "%s", %.1f)""" % (SurfaceName, MapName, SigmaLevel))
    else:
        PMLCmds.append("""cmd.isosurface("%s", "%s", %.1f, "(%s)", carve = %.1f)""" % (SurfaceName, MapName, SigmaLevel, Selection, Carve))
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, SurfaceName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(SurfaceName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def GeneratePyMOLSessionFile():
    """Generate PME file from PML file. """

    PSEOutfile = OptionsInfo["PSEOutfile"]
    PMLOutfile = OptionsInfo["PMLOutfile"]
    
    MiscUtil.PrintInfo("\nGenerating file %s..." % PSEOutfile)
    
    PyMOLUtil.ConvertPMLFileToPSEFile(PMLOutfile, PSEOutfile)
    
    if not os.path.exists(PSEOutfile):
        MiscUtil.PrintWarning("Failed to generate PSE file, %s..." % (PSEOutfile))
    
    if not OptionsInfo["PMLOut"]:
        MiscUtil.PrintInfo("Deleting file %s..." % PMLOutfile)
        os.remove(PMLOutfile)

def DeleteEmptyPyMOLObjects(OutFH, PyMOLObjectNames):
    """Delete empty PyMOL objects. """
    
    if OptionsInfo["AllowEmptyObjects"]:
        return
        
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        OutFH.write("""\n""\n"Checking and deleting empty objects for chain %s..."\n""\n""" % (ChainID))
        # Delete any chain level objects...
        WritePMLToCheckAndDeleteEmptyObject(OutFH, PyMOLObjectNames["Chains"][ChainID]["Solvent"])
        WritePMLToCheckAndDeleteEmptyObject(OutFH, PyMOLObjectNames["Chains"][ChainID]["Inorganic"])
        
        for LigandID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]:
            # Delete ligand level objects...
            for GroupID in ["Pocket", "PocketSolvent", "PocketInorganic"]:
                GroupNameID = "%sGroup" % (GroupID)
                GroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID]

                GroupTypeObjectID = "%s" % (GroupID)
                GroupTypeObjectName = PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID]
                
                WritePMLToCheckAndDeleteEmptyObject(OutFH, GroupTypeObjectName, GroupName)

def WritePMLToCheckAndDeleteEmptyObject(OutFH, ObjectName, ParentObjectName = None):
    """Write PML to check and delete empty PyMOL objects. """
    
    if ParentObjectName is None:
        PML = """CheckAndDeleteEmptyObject("%s")""" % (ObjectName)
    else:
        PML = """CheckAndDeleteEmptyObject("%s", "%s")""" % (ObjectName, ParentObjectName)
    
    OutFH.write("%s\n" % PML)

def SetupPyMOLObjectNames():
    """Setup hierarchy of PyMOL groups and objects for ligand centric views of
    electron density for chains and ligands present in input file.
    """

    PyMOLObjectNames = {}
    PyMOLObjectNames["Chains"] = {}
    PyMOLObjectNames["Ligands"] = {}

    # Setup groups and objects for complex...
    SetupPyMOLObjectNamesForComplex(PyMOLObjectNames)
    
    # Setup groups and objects for chain...
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        SetupPyMOLObjectNamesForChain(PyMOLObjectNames, ChainID)
        
        # Setup groups and objects for ligand...
        for LigandID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]:
            SetupPyMOLObjectNamesForLigand(PyMOLObjectNames, ChainID, LigandID)

    return PyMOLObjectNames

def SetupPyMOLObjectNamesForComplex(PyMOLObjectNames):
    """Stetup groups and objects for complex. """
    
    PDBFileRoot  = OptionsInfo["InfileRoot"]
    
    PDBGroupName = "%s" % PDBFileRoot
    PyMOLObjectNames["PDBGroup"] = PDBGroupName
    PyMOLObjectNames["PDBGroupMembers"] = []

    ComplexGroupName = "%s.Complex" % PyMOLObjectNames["PDBGroup"]
    PyMOLObjectNames["ComplexGroup"] = ComplexGroupName
    PyMOLObjectNames["PDBGroupMembers"].append(ComplexGroupName)
    
    PyMOLObjectNames["Complex"] = "%s.Complex" % ComplexGroupName

    CompositeMeshGroupName = "%s.2Fo-Fc" % (ComplexGroupName)
    CompositeMapName = "%s.Map" % (CompositeMeshGroupName)
    CompositeMeshName = "%s.Mesh" % (CompositeMeshGroupName)
    CompositeSurfaceName = "%s.Surface" % (CompositeMeshGroupName)
    
    PyMOLObjectNames["ComplexCompositeEDGroup"] = CompositeMeshGroupName
    PyMOLObjectNames["ComplexCompositeEDMap"] = CompositeMapName
    PyMOLObjectNames["ComplexCompositeEDMesh"] = CompositeMeshName
    PyMOLObjectNames["ComplexCompositeEDSurface"] = CompositeSurfaceName

    PyMOLObjectNames["ComplexCompositeEDGroupMembers"] = []
    PyMOLObjectNames["ComplexCompositeEDGroupMembers"].append(CompositeMapName)
    if OptionsInfo["MeshComplex"]:
        PyMOLObjectNames["ComplexCompositeEDGroupMembers"].append(CompositeMeshName)
    if OptionsInfo["SurfaceComplex"]:
        PyMOLObjectNames["ComplexCompositeEDGroupMembers"].append(CompositeSurfaceName)
    
    if OptionsInfo["DifferenceMapAvailable"]:
        DiffMeshGroupName = "%s.Fo-Fc" % ComplexGroupName
        DiffMapName = "%s.Map" % DiffMeshGroupName
        DiffMesh1Name = "%s.Mesh1" % DiffMeshGroupName
        DiffSurface1Name = "%s.Surface1" % DiffMeshGroupName
        DiffMesh2Name = "%s.Mesh2" % DiffMeshGroupName
        DiffSurface2Name = "%s.Surface2" % DiffMeshGroupName
        
        PyMOLObjectNames["ComplexDiffEDGroup"] = DiffMeshGroupName
        PyMOLObjectNames["ComplexDiffEDMap"] = DiffMapName
        PyMOLObjectNames["ComplexDiffEDMesh1"] = DiffMesh1Name
        PyMOLObjectNames["ComplexDiffEDSurface1"] = DiffSurface1Name
        PyMOLObjectNames["ComplexDiffEDMesh2"] = DiffMesh2Name
        PyMOLObjectNames["ComplexDiffEDSurface2"] = DiffSurface2Name
        
        PyMOLObjectNames["ComplexDiffEDGroupMembers"] = []
        PyMOLObjectNames["ComplexDiffEDGroupMembers"].append(DiffMapName)
        if OptionsInfo["MeshComplex"]:
            PyMOLObjectNames["ComplexDiffEDGroupMembers"].append(DiffMesh1Name)
        if OptionsInfo["SurfaceComplex"]:
            PyMOLObjectNames["ComplexDiffEDGroupMembers"].append(DiffSurface1Name)
        if OptionsInfo["MeshComplex"]:
            PyMOLObjectNames["ComplexDiffEDGroupMembers"].append(DiffMesh2Name)
        if OptionsInfo["SurfaceComplex"]:
            PyMOLObjectNames["ComplexDiffEDGroupMembers"].append(DiffSurface2Name)
    
    PyMOLObjectNames["ComplexGroupMembers"] = []
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["Complex"])
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["ComplexCompositeEDGroup"])
    if OptionsInfo["DifferenceMapAvailable"]:
        PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["ComplexDiffEDGroup"])
    
def SetupPyMOLObjectNamesForChain(PyMOLObjectNames, ChainID):
    """Setup groups and objects for chain."""
    
    PDBGroupName = PyMOLObjectNames["PDBGroup"]
    
    MeshChainComplex = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["MeshChainComplex"][ChainID]
    SurfaceChainComplex = OptionsInfo["SpecifiedChainsAndLigandsInfo"]["SurfaceChainComplex"][ChainID]
    
    PyMOLObjectNames["Chains"][ChainID] = {}
    PyMOLObjectNames["Ligands"][ChainID] = {}
    
    # Set up chain group and chain objects...
    ChainGroupName = "%s.Chain%s" % (PDBGroupName, ChainID)
    PyMOLObjectNames["Chains"][ChainID]["ChainGroup"] = ChainGroupName
    PyMOLObjectNames["PDBGroupMembers"].append(ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"] = []
    
    # Setup chain complex group and objects...
    ChainComplexGroupName = "%s.Complex" % (ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroup"] = ChainComplexGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainComplexGroupName)
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplex"] = "%s.Complex" % (ChainComplexGroupName)
    
    CompositeMeshGroupName = "%s.2Fo-Fc" % (ChainComplexGroupName)
    CompositeMeshName = "%s.Mesh" % (CompositeMeshGroupName)
    CompositeSurfaceName = "%s.Surface" % (CompositeMeshGroupName)
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroup"] = CompositeMeshGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDMesh"] = CompositeMeshName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDSurface"] = CompositeSurfaceName
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroupMembers"] = []
    if MeshChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroupMembers"].append(CompositeMeshName)
    if SurfaceChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCompositeEDGroupMembers"].append(CompositeSurfaceName)
    
    if OptionsInfo["DifferenceMapAvailable"]:
        DiffMeshGroupName = "%s.Fo-Fc" % (ChainComplexGroupName)
        DiffMesh1Name = "%s.Mesh1" % (DiffMeshGroupName)
        DiffSurface1Name = "%s.Surface1" % (DiffMeshGroupName)
        DiffMesh2Name = "%s.Mesh2" % (DiffMeshGroupName)
        DiffSurface2Name = "%s.Surface2" % (DiffMeshGroupName)
        
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroup"] = DiffMeshGroupName
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDMesh1"] = DiffMesh1Name
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDSurface1"] = DiffSurface1Name
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDMesh2"] = DiffMesh2Name
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDSurface2"] = DiffSurface2Name
        
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"] = []
        if MeshChainComplex:
            PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"].append(DiffMesh1Name)
        if SurfaceChainComplex:
            PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"].append(DiffSurface1Name)
        if MeshChainComplex:
            PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"].append(DiffMesh2Name)
        if SurfaceChainComplex:
            PyMOLObjectNames["Chains"][ChainID]["ChainComplexDiffEDGroupMembers"].append(DiffSurface2Name)
        
    NameIDs = ["ChainComplex"]
    if MeshChainComplex or SurfaceChainComplex :
        NameIDs.append("ChainComplexCompositeEDGroup")
        if OptionsInfo["DifferenceMapAvailable"]:
            NameIDs.append("ChainComplexDiffEDGroup")
        
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"] = []
    for NameID in NameIDs:
        Name = PyMOLObjectNames["Chains"][ChainID][NameID]
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexGroupMembers"].append(Name)

    # Setup up a group for individual chains...
    ChainAloneGroupName = "%s.Chain" % (ChainGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroup"] = ChainAloneGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainAloneGroupName)
        
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"] = []
        
    Name = "%s.Chain" % (ChainAloneGroupName)
    PyMOLObjectNames["Chains"][ChainID]["ChainAlone"] = Name
    PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(Name)
        
    if OptionsInfo["BFactorChainCartoonPutty"]:
        Name = "%s.BFactor" % (ChainAloneGroupName)
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneBFactorPutty"] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainAloneGroupMembers"].append(Name)
    
    # Setup solvent and inorganic objects for chain...
    for NameID in ["Solvent", "Inorganic"]:
        Name = "%s.%s" % (ChainGroupName, NameID)
        PyMOLObjectNames["Chains"][ChainID][NameID] = Name
        PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(Name)

def SetupPyMOLObjectNamesForLigand(PyMOLObjectNames, ChainID, LigandID):
    """Stetup groups and objects for ligand."""

    PyMOLObjectNames["Ligands"][ChainID][LigandID] = {}
    
    ChainGroupName = PyMOLObjectNames["Chains"][ChainID]["ChainGroup"]
    
    # Setup a chain level ligand group...
    ChainLigandGroupName = "%s.Ligand%s" % (ChainGroupName, LigandID)
    PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroup"] = ChainLigandGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainGroupMembers"].append(ChainLigandGroupName)
    
    PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"] = []

    # Set up groups and objects for a specific ligand group...
    for GroupType in ["Ligand", "Pocket", "Pocket_Solvent", "Pocket_Inorganic"]:
        GroupID = re.sub("_", "", GroupType)
        GroupName = "%s.%s" % (ChainLigandGroupName, GroupType)
                
        GroupNameID = "%sGroup" % (GroupID)
        GroupMembersID = "%sGroupMembers" % (GroupID)
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupNameID] = GroupName
        PyMOLObjectNames["Ligands"][ChainID][LigandID]["ChainLigandGroupMembers"].append(GroupName)
        
        GroupTypeObjectName = "%s.%s" % (GroupName, GroupType)
        GroupTypeObjectID = "%s" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupTypeObjectID] = GroupTypeObjectName
        
        CompositeMeshGroupName = "%s.2Fo-Fc" % (GroupName)
        CompositeMeshName = "%s.Mesh" % (CompositeMeshGroupName)
        CompositeSurfaceName = "%s.Surface" % (CompositeMeshGroupName)
                
        CompositeMeshGroupID = "%sCompositeEDMeshGroup" % (GroupID)
        CompositeMeshGroupMembersID = "%sCompositeEDMeshGroupMembers" % (GroupID)
        CompositeMeshID = "%sCompositeEDMesh" % (GroupID)
        CompositeSurfaceID = "%sCompositeEDSurface" % (GroupID)
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupID] = CompositeMeshGroupName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshID] = CompositeMeshName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeSurfaceID] = CompositeSurfaceName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupMembersID] = []
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupMembersID].append(CompositeMeshName)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CompositeMeshGroupMembersID].append(CompositeSurfaceName)
                
        if OptionsInfo["DifferenceMapAvailable"]:
            DiffMeshGroupName = "%s.Fo-Fc" % GroupName
            DiffMesh1Name = "%s.Mesh1" % DiffMeshGroupName
            DiffSurface1Name = "%s.Surface1" % DiffMeshGroupName
            DiffMesh2Name = "%s.Mesh2" % DiffMeshGroupName
            DiffSurface2Name = "%s.Surface2" % DiffMeshGroupName
                    
            DiffMeshGroupID = "%sDiffEDMeshGroup" % (GroupID)
            DiffMeshGroupMembersID = "%sDiffEDMeshGroupMembers" % (GroupID)
            DiffMesh1ID = "%sDiffEDMesh1" % (GroupID)
            DiffSurface1ID = "%sDiffEDSurface1" % (GroupID)
            DiffMesh2ID = "%sDiffEDMesh2" % (GroupID)
            DiffSurface2ID = "%sDiffEDSurface2" % (GroupID)
                    
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMeshGroupID] = DiffMeshGroupName
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMesh1ID] = DiffMesh1Name
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffSurface1ID] = DiffSurface1Name
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMesh2ID] = DiffMesh2Name
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffSurface2ID] = DiffSurface2Name
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMeshGroupMembersID] = []
            PyMOLObjectNames["Ligands"][ChainID][LigandID][DiffMeshGroupMembersID].extend([DiffMesh1Name, DiffSurface1Name, DiffMesh2Name, DiffSurface2Name])
                
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID] = []
        NameIDs = [GroupTypeObjectID, CompositeMeshGroupID]
        if OptionsInfo["DifferenceMapAvailable"]:
            NameIDs.append(DiffMeshGroupID)
        
        for NameID in NameIDs:
            Name = PyMOLObjectNames["Ligands"][ChainID][LigandID][NameID]
            PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(Name)
        
        if re.match("^Ligand$", GroupType, re.I):
            # No other object needed for Ligand group...
            continue
        
        PolarContactsName = "%s.Polar_Contacts" % (GroupName)
        PolarContactsID = "%sPolarContacts" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][PolarContactsID] = PolarContactsName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(PolarContactsName)
                
        if not re.match("^Pocket$", GroupType, re.I):
            # No other object needed for any other group besides Pocket...
            continue
        
        if not OptionsInfo["PocketSurface"]:
            continue
        
        HydrophobicSurfaceName = "%s.Surface" % (GroupName)
        HydrophobicSurfaceID = "%sHydrophobicSurface" % (GroupID)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][HydrophobicSurfaceID] = HydrophobicSurfaceName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID].append(HydrophobicSurfaceName)

def ProcessEDMapSuffixes():
    """Process suffixes for ED map files"""
    
    OptionsInfo["EDMapSuffixesMap"] = {"CompositeMap" : "", "DifferenceMap" : "_diff"}
    
    if re.match("^auto$", OptionsInfo["EDMapSuffixes"], re.I):
        SetupEDMapFileNames()
        return
    
    EDMapSuffixes = re.sub(" ", "", OptionsInfo["EDMapSuffixes"])
    if not EDMapSuffixes:
        MiscUtil.PrintError("No valid parameter name and value pairs specified using \"--EDMapSuffixes\" option.")
    
    EDMapSuffixesWords = EDMapSuffixes.split(",")
    if len(EDMapSuffixesWords) % 2:
        MiscUtil.PrintError("The number of comma delimited ED map types names and suffixes, %d, specified using \"--EDMapSuffixes\" option must be an even number." % (len(EDMapSuffixesWords)))

    for Index in range(0, len(EDMapSuffixesWords), 2):
        EDMapType = EDMapSuffixesWords[Index]
        EDMapSuffix = EDMapSuffixesWords[Index + 1]
        
        if re.match("^CompositeMap$", EDMapType, re.I):
            EDMapType = "CompositeMap"
        elif re.match("^DifferenceMap$", EDMapType, re.I):
            EDMapType = "DifferenceMap"
        else:
            MiscUtil.PrintError("The ED map type, %s, specified using \"--EDMapSuffixes\" option is not a valid ED map type. Supported ED map types: CompositeMap, DifferenceMap" % (EDMapType))
            
        if re.match(EDMapSuffix, "None", re.I):
            EDMapSuffix = ""
        
        OptionsInfo["EDMapSuffixesMap"][EDMapType] = EDMapSuffix
    
    SetupEDMapFileNames()

def SetupEDMapFileNames():
    """Set up ED map file names. """
    
    PDBFileRoot = OptionsInfo["InfileRoot"]
    
    OptionsInfo["EDMapFiles"] = {}
    OptionsInfo["DifferenceMapAvailable"] = True
    
    for EDMapType in sorted(OptionsInfo["EDMapSuffixesMap"]):
        EDMapSuffix = OptionsInfo["EDMapSuffixesMap"][EDMapType]
        EDMapFile = "%s%s.ccp4" % (PDBFileRoot, EDMapSuffix)
        
        EDMapFileExists = True if os.path.exists(EDMapFile) else False
        if not os.path.exists(EDMapFile):
            if re.match("^DifferenceMap$", EDMapType, re.I):
                OptionsInfo["DifferenceMapAvailable"] = False
                MiscUtil.PrintWarning("The ED map file, %s, corresponding to EDMapType, %s, for option \"--EDMapSuffix\" doesn't exist. PyMOL groups and objectes related to difference maps won't be created." % (EDMapFile, EDMapType))
            else:
                MiscUtil.PrintError("The ED map file, %s, corresponding to EDMapType, %s, for option \"--EDMapSuffix\" doesn't exist.\n" % (EDMapFile, EDMapType))
            
        OptionsInfo["EDMapFiles"][EDMapType] = EDMapFile

def ProcessChainAndLigandIDs():
    """Process chain and ligand IDs"""

    MolName = OptionsInfo["InfileRoot"]
    ChainsAndLigandsInfo = PyMOLUtil.GetChainsAndLigandsInfo(OptionsInfo["Infile"], MolName)
    OptionsInfo["ChainsAndLigandsInfo"] = ChainsAndLigandsInfo
    
    SpecifiedChainsAndLigandsInfo = PyMOLUtil.ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, "-c, --chainIDs", OptionsInfo["ChainIDs"], "-l, --ligandIDs", OptionsInfo["LigandIDs"])
    OptionsInfo["SpecifiedChainsAndLigandsInfo"] = SpecifiedChainsAndLigandsInfo
    
    CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo)
    
    ProcessChainMeshesAndSurfacesOptions()

def CheckPresenceOfValidLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo):
    """Check presence of valid ligand IDs."""

    MiscUtil.PrintInfo("\nSpecified chain IDs: %s" % (", ".join(SpecifiedChainsAndLigandsInfo["ChainIDs"])))
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        if len (SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID]):
            MiscUtil.PrintInfo("Chain ID: %s; Specified LigandIDs: %s" % (ChainID, ", ".join(SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID])))
        else:
            MiscUtil.PrintInfo("Chain IDs: %s; Specified LigandIDs: None" % (ChainID))
            MiscUtil.PrintWarning("No valid ligand IDs found for chain ID, %s. PyMOL groups and objects related to ligand and binding pockect won't be created." % (ChainID))

def ProcessChainMeshesAndSurfacesOptions():
    """Process options to create meshes and surfaces for chains."""

    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["MeshChainComplex"] = {}
    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["SurfaceChainComplex"] = {}
    
    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainComplexGroup"] = {}
    OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainAloneGroup"] = {}
    
    for ChainID in OptionsInfo["SpecifiedChainsAndLigandsInfo"]["ChainIDs"]:
        LigandsPresent = True if len(OptionsInfo["SpecifiedChainsAndLigandsInfo"]["LigandIDs"][ChainID]) else False

        if re.match("^auto$", OptionsInfo["MeshChainComplex"], re.I):
            MeshChainComplex = False if LigandsPresent else True
        else:
            MeshChainComplex = True if re.match("^Yes$", OptionsInfo["MeshChainComplex"], re.I) else False
        
        if re.match("^auto$", OptionsInfo["SurfaceChainComplex"], re.I):
            if re.match("^auto$", OptionsInfo["MeshChainComplex"], re.I):
                SurfaceChainComplex = False
            else:
                SurfaceChainComplex = False if (LigandsPresent or MeshChainComplex) else True
        else:
            SurfaceChainComplex = True if re.match("^Yes$", OptionsInfo["SurfaceChainComplex"], re.I) else False
        
        if LigandsPresent:
            EnableChainComplexGroup = False
            EnableChainAloneGroup = True
        else:
            EnableChainComplexGroup = True
            EnableChainAloneGroup = False

        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["MeshChainComplex"][ChainID] = MeshChainComplex
        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["SurfaceChainComplex"][ChainID] = SurfaceChainComplex
        
        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainComplexGroup"][ChainID] = EnableChainComplexGroup
        OptionsInfo["SpecifiedChainsAndLigandsInfo"]["EnableChainAloneGroup"][ChainID] = EnableChainAloneGroup

def ProcessOptions():
    """Process and validate command line arguments and options"""
    
    MiscUtil.PrintInfo("Processing options...")
    
    # Validate options...
    ValidateOptions()
    
    OptionsInfo["AllowEmptyObjects"] = True if re.match("^Yes$", Options["--allowEmptyObjects"], re.I) else False
    
    OptionsInfo["BFactorChainCartoonPutty"] = True if re.match("^Yes$", Options["--BFactorChainCartoonPutty"], re.I) else False
    OptionsInfo["BFactorColorPalette"] = Options["--BFactorColorPalette"]
    
    OptionsInfo["Infile"] = Options["--infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Infile"])
    OptionsInfo["InfileRoot"] = FileName

    OptionsInfo["Overwrite"] = Options["--overwrite"]
    OptionsInfo["PMLOut"] = True if re.match("^Yes$", Options["--PMLOut"], re.I) else False
    
    OptionsInfo["Outfile"] = Options["--outfile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(OptionsInfo["Outfile"])
    OptionsInfo["PSEOut"] = False 
    if re.match("^pml$", FileExt, re.I):
        OptionsInfo["PMLOutfile"] = OptionsInfo["Outfile"] 
        OptionsInfo["PMEOutfile"] = re.sub(".pml$", ".pme", OptionsInfo["Outfile"]) 
    elif re.match("^pse$", FileExt, re.I):
        OptionsInfo["PSEOut"] = True 
        OptionsInfo["PSEOutfile"] = OptionsInfo["Outfile"] 
        OptionsInfo["PMLOutfile"] = re.sub(".pse$", ".pml", OptionsInfo["Outfile"]) 
        if os.path.exists(OptionsInfo["PMLOutfile"]) and (not OptionsInfo["Overwrite"]):
            MiscUtil.PrintError("The intermediate output file to be generated, %s, already exist. Use option \"--ov\" or \"--overwrite\" and try again." % OptionsInfo["PMLOutfile"] )

    OptionsInfo["EDMapSuffixes"] = Options["--EDMapSuffixes"]
    ProcessEDMapSuffixes()
    
    OptionsInfo["LabelFontID"] = int(Options["--labelFontID"])
    
    # Process mesh parameters...
    OptionsInfo["MeshCarveRadius"] = float(Options["--meshCarveRadius"])
    OptionsInfo["MeshComplex"] = True if re.match("^Yes$", Options["--meshComplex"], re.I) else False
    OptionsInfo["MeshChainComplex"] = Options["--meshChainComplex"]
    OptionsInfo["MeshWidth"] = float(Options["--meshWidth"])
    
    OptionsInfo["MeshColorCompositeMap"] = Options["--meshColorCompositeMap"]
    OptionsInfo["MeshLevelCompositeMap"] = float(Options["--meshLevelCompositeMap"])
    OptionsInfo["Mesh1ColorDiffMap"] = Options["--mesh1ColorDiffMap"]
    OptionsInfo["Mesh1LevelDiffMap"] = float(Options["--mesh1LevelDiffMap"])
    OptionsInfo["Mesh2ColorDiffMap"] = Options["--mesh2ColorDiffMap"]
    OptionsInfo["Mesh2LevelDiffMap"] = float(Options["--mesh2LevelDiffMap"])
    
    OptionsInfo["SurfaceComplex"] = True if re.match("^Yes$", Options["--surfaceComplex"], re.I) else False
    OptionsInfo["SurfaceChainComplex"] = Options["--surfaceChainComplex"]
    OptionsInfo["SurfaceTransparency"] = float(Options["--surfaceTransparency"])
    
    OptionsInfo["PocketContactsLigandColor"] = Options["--pocketContactsLigandColor"]
    OptionsInfo["PocketContactsSolventColor"] = Options["--pocketContactsSolventColor"]
    OptionsInfo["PocketContactsInorganicColor"] = Options["--pocketContactsInorganicColor"]
    
    OptionsInfo["PocketDistanceCutoff"] = float(Options["--pocketDistanceCutoff"])
    OptionsInfo["PocketLabelColor"] = Options["--pocketLabelColor"]
    OptionsInfo["PocketSurface"] = True if re.match("^Yes$", Options["--pocketSurface"], re.I) else False
    
    OptionsInfo["ChainIDs"] = Options["--chainIDs"]
    OptionsInfo["LigandIDs"] = Options["--ligandIDs"]
    ProcessChainAndLigandIDs()

def RetrieveOptions(): 
    """Retrieve command line arguments and options"""
    
    # Get options...
    global Options
    Options = docopt(_docoptUsage_)
    
    # Set current working directory to the specified directory...
    WorkingDir = Options["--workingdir"]
    if WorkingDir:
        os.chdir(WorkingDir)
    
    # Handle examples option...
    if "--examples" in Options and Options["--examples"]:
        MiscUtil.PrintInfo(MiscUtil.GetExamplesTextFromDocOptText(_docoptUsage_))
        sys.exit(0)

def ValidateOptions():
    """Validate option values"""
    
    MiscUtil.ValidateOptionTextValue("--allowEmptyObjects", Options["--allowEmptyObjects"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--BFactorChainCartoonPutty", Options["--BFactorChainCartoonPutty"], "yes no")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb cif")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pml pse")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionTextValue("--PMLOut", Options["--PMLOut"], "yes no")
    MiscUtil.ValidateOptionIntegerValue("--labelFontID", Options["--labelFontID"], {})

    MiscUtil.ValidateOptionFloatValue("--meshCarveRadius", Options["--meshCarveRadius"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--meshComplex", Options["--meshComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--meshChainComplex", Options["--meshChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--meshWidth", Options["--meshWidth"], {">": 0.0})
    
    MiscUtil.ValidateOptionFloatValue("--meshLevelCompositeMap", Options["--meshLevelCompositeMap"], {})
    MiscUtil.ValidateOptionFloatValue("--mesh1LevelDiffMap", Options["--mesh1LevelDiffMap"], {})
    MiscUtil.ValidateOptionFloatValue("--mesh2LevelDiffMap", Options["--mesh2LevelDiffMap"], {})
    
    MiscUtil.ValidateOptionTextValue("--surfaceComplex", Options["--surfaceComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChainComplex", Options["--surfaceChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--surfaceTransparency", Options["--surfaceTransparency"], {">=": 0.0, "<=": 1.0})
    
    MiscUtil.ValidateOptionFloatValue("--pocketDistanceCutoff", Options["--pocketDistanceCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--pocketSurface", Options["--pocketSurface"], "yes no")
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLVisualizeElectronDensity.py - Visualize electron density

Usage:
    PyMOLVisualizeElectronDensity.py [--allowEmptyObjects <yes or no>]
                                     [--BFactorChainCartoonPutty <yes or no>] [--BFactorColorPalette <text> ]
                                     [--chainIDs <First, All or ID1,ID2...>] [--EDMapSuffixes <CompositeMap,None,...>]
                                     [--ligandIDs <Largest, All or ID1,ID2...>] [--labelFontID <number>]
                                     [--meshCarveRadius <number>] [--meshComplex <yes or no>]
                                     [--meshChainComplex <yes or no>] [--meshColorCompositeMap <text>]
                                     [--meshLevelCompositeMap <number>] [--meshWidth <number>]
                                     [--mesh1ColorDiffMap <text>] [--mesh1LevelDiffMap <number>]
                                     [--mesh2ColorDiffMap <text>] [--mesh2LevelDiffMap <number>]
                                     [--PMLOut <yes or no>] [--pocketContactsLigandColor <text>]
                                     [--pocketContactsSolventColor <text>] [--pocketContactsInorganicColor <text>]
                                     [--pocketDistanceCutoff <number>] [--pocketLabelColor <text>]
                                     [--pocketSurface <yes or no>] [--surfaceComplex <yes or no> ] [--surfaceChainComplex <yes or no>]
                                     [--surfaceTransparency <number>] [--overwrite] [-w <dir>] -i <infile> -o <outfile>
    PyMOLVisualizeElectronDensity.py -h | --help | -e | --examples

Description:
    Generate PyMOL visualization files for viewing X-ray electron density around
    chains, ligands, and ligand binding pockets in macromolecules including proteins
    and nucleic acids.

    The supported input file formats are: Macromolecule - PDB (.pdb) or CIF(.cif),
    Electron Density - Collaborative Computational Project Number 4 (CCP4) ( .ccp4)

    The supported output file formats are: PyMOL script file (.pml), PyMOL session
    file (.pse)

    Two types of CCP4 electron density map files may be used for visualizing electron
    density. These file types along with default file names are shown below:
    
        CompositeMap (2Fobs - Fcalc) - <InfileRoot>.ccp4 (required)
        DifferenceMap (Fobs - Fcalc) - <InfileRoot>_diff.ccp4 (optional)
    
    The compsite map file must be present. The difference map file is optional.
    The mesh and surface PyMOL objects are not generated for missing difference
    map file.

    The electron density present in composite map file is generated by adding two
    difference maps to a calculated map (Fcalc) as shown below:
    
        Fcalc + 2(Fobs - Fcalc) = 2Fobs - Fcalc
    
    The following types of meshes may be created by default for electron density
    present in composite and difference map files:
    
        CompositeMesh - ContourLevel: 1; Color: Blue
        DiffMesh1 - ContourLevel: 3; Color: Green
        DiffMesh2 - ContourLevel: -3; Color: Red

    The two meshes created for difference maps correspond to false negative and
    false positive in terms of electron density present in the model. The first mesh
    shown in  green color corresponds to observed electron density missing in the
    model. The second mesh in in red color indicates model electron density not
    observed in the experiment.
    
    A variety of PyMOL groups and objects may be  created for visualization of
    electron density present in map files. These groups and objects correspond to
    maps, meshes, surfaces,chains, ligands, inorganics, ligand binding pockets,
    pocket, polar interactions, and pocket hydrophobic surfaces. A complete
    hierarchy of all possible PyMOL groups and objects is shown below:
    
        <PDBFileRoot>
            .Complex
                .Complex
                .2Fo-Fc
                    .Map
                    .Mesh
                    .Surface
                .Fo-Fc
                    .Map
                    .Mesh1
                    .Surface1
                    .Mesh2
                    .Surface2
            .Chain<ID>
                .Complex
                    .Complex
                    .2Fo-Fc
                        .Mesh
                        .Surface
                    .Fo-Fc
                        .Mesh1
                        .Surface1
                        .Mesh2
                        .Surface2
                .Chain
                    .Chain
                    .BFactor
                .Solvent
                .Inorganic
                .Ligand<ID>
                    .Ligand
                        .Ligand
                        .2Fo-Fc
                            .Mesh
                            .Surface
                        .Fo-Fc
                            .Mesh1
                            .Surface1
                            .Mesh2
                            .Surface2
                    .Pocket
                        .Pocket
                        .2Fo-Fc
                            .Mesh
                            .Surface
                        .Fo-Fc
                            .Mesh1
                            .Surface1
                            .Mesh2
                            .Surface2
                        .Polar_Contacts
                        .Surface
                    .Pocket_Solvent
                        .Pocket_Solvent
                        .2Fo-Fc
                            .Mesh
                            .Surface
                        .Fo-Fc
                            .Mesh1
                            .Surface1
                            .Mesh2
                            .Surface2
                        .Polar_Contacts
                    .Pocket_Inorganic
                        .Pocket_Inorganic
                        .2Fo-Fc
                            .Mesh
                            .Surface
                        .Fo-Fc
                            .Mesh1
                            .Surface1
                            .Mesh2
                            .Surface2
                        .Polar_Contacts
                .Ligand<ID>
                    .Ligand
                        ... ... ...
                    .Pocket
                        ... ... ...
                    .Pocket_Solvent
                        ... ... ...
                    .Pocket_Inorganic
                        ... ... ...
            .Chain<ID>
                ... ... ...
                .Ligand<ID>
                    ... ... ...
                .Ligand<ID>
                    ... ... ...
            .Chain<ID>
                ... ... ...
    
    The meshes and surfaces  are not created for complete complex in input file
    by default. A word to the wise: The creation of these surface and mesh objects
    may slow down loading of PML file and generation of PSE file, based on the size
    of input complex and map files. The generation of PSE file may also fail.

Options:
    -a, --allowEmptyObjects <yes or no>  [default: no]
        Allow creation of empty PyMOL objects corresponding to solvent and
        inorganic atom selections across chains, ligands, and ligand binding pockets
        in input file.
    -c, --chainIDs <First, All or ID1,ID2...>  [default: First]
        List of chain IDs to use for visualizing electron density. Possible values:
        First, All, or a comma delimited list of chain IDs. The default is to use the
        chain ID for the first chain in input file.
    -b, --BFactorChainCartoonPutty <yes or no>  [default: yes]
        A cartoon putty around individual chains colored by B factors. The minimum
        and maximum values for B factors are automatically detected. These values
        indicate spread of electron density around atoms. The 'blue_white_red' color
        palette is deployed for coloring the cartoon putty.
    --BFactorColorPalette <text>  [default: blue_white_red]
        Color palette for coloring cartoon putty around chains generated using B
        factors. An valid PyMOL color palette name is allowed. No validation is
        performed. The complete list of valid color palette names is a available
        at: pymolwiki.org/index.php/Spectrum. Examples: blue_white_red,
        blue_white_magenta, blue_red, green_white_red, green_red.
    -e, --examples
        Print examples.
    --EDMapSuffixes <CompositeMap,None,...>  [default: auto]
        Electron density map file suffixes for generating names of map files from
        the root of input file. It is a pairwise comma delimited list of 'EDMapType'
        and file suffix.
        
        Supported values for 'EDMapType': 'CompositeMap, DifferenceMap'.
        Supported value for file suffix: Any valid string.
        
        Default value: 'CompositeMap,None,DifferenceMap,_diff'
        
        The default names of the map files, generated form a combination of
        'InfileRoot' and 'EDSMapType' are shown below:
            
            CompositeMap (2Fobs - Fcalc) - <InfileRoot>.ccp4
            DifferenceMap (Fobs - Fcalc) - <InfileRoot>_diff.ccp4
            
        The composite map file must be present. The difference map file is
        optional.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -l, --ligandIDs <Largest, All or ID1,ID2...>  [default: Largest]
        List of ligand IDs present in chains for visualizing electron density across
        ligands and ligand binding pockets. Possible values: Largest, All, or a comma
        delimited list of ligand IDs. The default is to use the largest ligand present
        in all or specified chains in input file.
        
        Ligands are identified using organic selection operator available in PyMOL.
        It'll also  identify buffer molecules as ligands. The largest ligand contains
        the highest number of heavy atoms.
    --labelFontID <number>  [default: 7]
        Font ID for drawing labels. Default: 7 (Sans Bold). Valid values: 5 to 16.
        The specified value must be a valid PyMOL font ID. No validation is
        performed. The complete lists of valid font IDs is available at:
        pymolwiki.org/index.php/Label_font_id. Examples: 5 - Sans;
        7 - Sans Bold; 9 - Serif; 10 - Serif Bold.
    --meshCarveRadius <number>  [default: 1.6]
        Radius in Angstroms around atoms for including electron density.
    --meshComplex <yes or no>  [default: no]
        Create meshes for complete complex in input file corresponding to composite
        and difference maps. A total of three meshes, one for composite map and
        two for difference map, are created for the complete complex.
        
        The composite and difference maps are always loaded for the complex.
    --meshChainComplex <yes, no, or auto>  [default: auto]
        Create meshes for individual chain complex in input file corresponding to
        composite and difference maps. A total of three meshes, one for composite map
        map and two for difference map, are created for each chain complex. By default,
        the meshes are automatically created for chain complexes without any ligands. 
    --meshColorCompositeMap <text>  [default: blue]
        Line color for mesh corresponding to composite map. The specified value
        must be valid color. No validation is performed.
    --meshLevelCompositeMap <number>  [default: 1.0]
        Contour level in sigma units for generating mesh corresponding to composite
        map.
    --meshWidth <number>  [default: 0.5]
        Line width for mesh lines corresponding to composite and difference maps.
    --mesh1ColorDiffMap <text>  [default: green]
        Line color for first mesh corresponding to difference map at contour level
        specified by '--mesh1LevelDiffMap'. The specified value must be valid color.
        No validation is performed.
    --mesh1LevelDiffMap <number>  [default: 3.0]
        Contour level in sigma units for generating first mesh corresponding to 
        to  difference map.
    --mesh2ColorDiffMap <text>  [default: red]
        Line color for second mesh corresponding to difference map at contour level
        specified by '--mesh2LevelDiffMap'. The specified value must be valid color.
        No validation is performed.
    --mesh2LevelDiffMap <number>  [default: -3.0]
        Contour level in sigma units for generating second mesh corresponding to
        difference map.
    -o, --outfile <outfile>
        Output file name.
    -p, --PMLOut <yes or no>  [default: yes]
        Save PML file during generation of PSE file.
    --pocketContactsLigandColor <text>  [default: orange]
        Color for drawing polar contacts between ligand and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketContactsSolventColor <text>  [default: marine]
        Color for drawing polar contacts between solvent and pocket residues..
        The specified value must be valid color. No validation is performed.
    --pocketContactsInorganicColor <text>  [default: deepsalmon]
        Color for drawing polar contacts between inorganic and pocket residues.
        The specified value must be valid color. No validation is performed.
    --pocketDistanceCutoff <number>  [default: 5.0]
        Distance in Angstroms for identifying pocket residues around ligands.
    --pocketLabelColor <text>  [default: magenta]
        Color for drawing residue or atom level labels for a pocket. The specified
        value must be valid color. No validation is performed.
    --pocketSurface <yes or no>  [default: yes]
        Hydrophobic surface around pocket. The pocket surface is colored by
        hydrophobicity. It is only valid for proteins. The color of amino acids is
        set using the Eisenberg hydrophobicity scale. The color varies from red
        to white, red being the most hydrophobic amino acid.
    --surfaceComplex <yes or no>  [default: no]
        Create surfaces for complete complex in input file corresponding to composite
        and difference maps. A total of three surfaces, one for composite map and
        two for difference map, are created for the complete complex.
        
        The composite and difference maps are always loaded for the complex.
    --surfaceChainComplex <yes, no or auto>  [default: auto]
        Create surfaces for individual chain complexes in input file corresponding to
        composite and difference maps. A total of three surfaces, one for composite
        map and two for difference map, are created for each chain complex. By default,
        the surfaces are automatically created for chain complexes without any ligands. 
    --surfaceTransparency <number>  [default: 0.25]
        Surface transparency for molecular and electron density surfaces.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To visualize electron density for the largest ligand in the first chain, and
    ligand binding pockets to highlight ligand interactions with pockect residues,
    solvents and inorganics, in a PDB and using corresponding map files, and
    generate a PML file, type:

        % PyMOLVisualizeElectronDensity.py -i Sample3.pdb -o Sample3.pml

    To visualize electron density for all ligands in all chains, and ligand binding
    pockets to highlight ligand interactions with pockect residues, solvents
    and inorganics, in a PDB and using corresponding map files, and generate
    a PML file, type:

        % PyMOLVisualizeElectronDensity.py -i Sample3.pdb -o Sample3.pml
          -c All -l All

    To visualize electron density for all chains and ligands, along with displaying
    meshes and surfaces for complete complex and individual chains, in a PDB
    file and using  corresponding map files, and generate a PML file, type:

        % PyMOLVisualizeElectronDensity.py -i Sample3.pdb -o Sample3.pml
          --chainIDs All --ligandIDs All --meshComplex yes --surfaceComplex yes
          --meshChainComplex yes --surfaceChainComplex yes

    To visualize electron density for ligand ADP in chain E along with ligand binding
    pocket, in a PDB and using  corresponding map files, and generate a PSE file, type:

        % PyMOLVisualizeElectronDensity.py -i Sample3.pdb -o Sample3.pse
          --chainIDs E --ligandIDs ADP

    To visualize electron density for all igands in all chains along with their binding
    pockets in a PDB file and using explicit file name suffixes for map files, and
    generate a PML file, type:

        % PyMOLVisualizeElectronDensity.py -i Sample3.pdb -o Sample3.pml
          --chainIDs All --ligandIDs All --EDMapSuffixes "CompositeMap,None,
          DifferenceMap,_diff"

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl, PyMOLVisualizeCryoEMDensity.py,
    PyMOLVisualizeMacromolecules.py

Copyright:
    Copyright (C) 2018 Manish Sud. All rights reserved.

    The functionality available in this script is implemented using PyMOL, a
    molecular visualization system on an open source foundation originally
    developed by Warren DeLano.

    This file is part of MayaChemTools.

    MayaChemTools is free software; you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option) any
    later version.

"""

if __name__ == "__main__":
    main()
