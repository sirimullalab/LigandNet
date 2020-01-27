#!/bin/env python
#
# File: PyMOLVisualizeCryoEMDensity.py
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
import xml.etree.ElementTree as ElementTree

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
    GenerateCryoEMDensityVisualization()
    
    MiscUtil.PrintInfo("\n%s: Done...\n" % ScriptName)
    MiscUtil.PrintInfo("Total time: %s" % MiscUtil.GetFormattedElapsedTime(WallClockTime, ProcessorTime))

def GenerateCryoEMDensityVisualization():
    """Generate cryoEM density visualization."""
    
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
    """Write out PML for viewing polymer complex along with cryo-EM density."""

    # Setup complex...
    PML = PyMOLUtil.SetupPMLForPolymerComplexView(PyMOLObjectNames["Complex"], OptionsInfo["Infile"], True)
    OutFH.write("""\n""\n"Loading %s and setting up view for complex..."\n""\n""" % OptionsInfo["Infile"])
    OutFH.write("%s\n" % PML)

    # Setup cryo-EM density maps and meshes...
    WriteComplexCryoEMDensityMapView(OutFH, PyMOLObjectNames)

    # Setup complex group...
    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexGroup"], PyMOLObjectNames["ComplexGroupMembers"], False, "close")

def WriteComplexCryoEMDensityMapView(OutFH, PyMOLObjectNames):
    """Write out PML for viewing cryoEM density map."""

    # Load cryo-EM density map and setup mesh views...
    MapFileName = OptionsInfo["DensityMapFileName"]
    Info = """\
""
"Loading cryo-EM density map %s and setting up mesh view for complex..."
"" """ % MapFileName
    OutFH.write("\n%s\n" % Info)

    MapName = PyMOLObjectNames["ComplexCryoEMMap"]
    ComplexName = PyMOLObjectNames["Complex"]
    
    ContourLevel = OptionsInfo["MeshLevel"]
    Color = OptionsInfo["MeshColor"]
    
    MeshName = PyMOLObjectNames["ComplexCryoEMMesh"]
    SurfaceName = PyMOLObjectNames["ComplexCryoEMSurface"]
    
    PML = SetupPMLForCryoEMDensityMap(MapFileName, MapName, True)
    OutFH.write("%s\n" % PML)

    EnableSurface = False if OptionsInfo["MeshComplex"] else True
    if OptionsInfo["MeshComplex"]:
        PML = SetupPMLForCryoEMDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = True, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    if OptionsInfo["SurfaceComplex"]:
        PML = SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurface, Selection = ComplexName)
        OutFH.write("\n%s\n" % PML)

    GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["ComplexCryoEMGroup"], PyMOLObjectNames["ComplexCryoEMGroupMembers"], True, "close")
    
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
        # Set up cryoEM mesh and group...
        MapName = PyMOLObjectNames["ComplexCryoEMMap"]
        ContourLevel = OptionsInfo["MeshLevel"]
        Color = OptionsInfo["MeshColor"]
        
        MeshName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMMesh"]
        SurfaceName = PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMSurface"]
        
        if MeshChainComplex:
            PML = SetupPMLForCryoEMDensityMesh(MapName, MeshName, ContourLevel, Color, Enable = True, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        if SurfaceChainComplex:
            PML = SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, ContourLevel, Color, Enable = EnableSurface, Selection = ChainComplexName)
            OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroup"], PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"], True, "close")
        
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
        
        # Set up cryoEM mesh and group...
        CryoEMMeshGroupID = "%sCryoEMMeshGroup" % (GroupID)
        CryoEMMeshGroupMembersID = "%sCryoEMMeshGroupMembers" % (GroupID)
        CryoEMMeshID = "%sCryoEMMesh" % (GroupID)
        CryoEMSurfaceID = "%sCryoEMSurface" % (GroupID)

        CryoEMMapName = PyMOLObjectNames["ComplexCryoEMMap"]
        CryoEMMeshName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshID]
        CryoEMSurfaceName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMSurfaceID]
        CryoEMMeshGroupName = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupID]
        CryoEMMeshGroupMembers = PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID]
        
        PML = SetupPMLForCryoEMDensityMesh(CryoEMMapName, CryoEMMeshName, OptionsInfo["MeshLevel"], OptionsInfo["MeshColor"], Enable = True, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        PML = SetupPMLForCryoEMDensitySurface(CryoEMMapName, CryoEMSurfaceName, OptionsInfo["MeshLevel"], OptionsInfo["MeshColor"], Enable = False, Selection = GroupTypeObjectName)
        OutFH.write("\n%s\n" % PML)
        
        GenerateAndWritePMLForGroup(OutFH, CryoEMMeshGroupName, CryoEMMeshGroupMembers, True, "close")
        
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

def SetupPMLForCryoEMDensityMap(MapFileName, MapName, Enable = True):
    """Setup PML for loading and viewing cryo-EM density map. """

    PMLCmds = []
    PMLCmds.append("""cmd.load("%s", "%s")""" % (MapFileName, MapName))
    PMLCmds.append(PyMOLUtil.SetupPMLForEnableDisable(MapName, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML
    
def SetupPMLForCryoEMDensityMesh(MapName, MeshName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for cryo-EM density mesh. """

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

def SetupPMLForCryoEMDensitySurface(MapName, SurfaceName, SigmaLevel, Color, Enable = True, Selection = None):
    """Setup PML for cryo-EM density surface. """

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
    cryo-EM density for chains and ligands present in input file.
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

    CryoEMMeshGroupName = "%s.CryoEM" % (ComplexGroupName)
    CryoEMMapName = "%s.Map" % (CryoEMMeshGroupName)
    CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
    CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
    
    PyMOLObjectNames["ComplexCryoEMGroup"] = CryoEMMeshGroupName
    PyMOLObjectNames["ComplexCryoEMMap"] = CryoEMMapName
    PyMOLObjectNames["ComplexCryoEMMesh"] = CryoEMMeshName
    PyMOLObjectNames["ComplexCryoEMSurface"] = CryoEMSurfaceName

    PyMOLObjectNames["ComplexCryoEMGroupMembers"] = []
    PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMMapName)
    if OptionsInfo["MeshComplex"]:
        PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMMeshName)
    if OptionsInfo["SurfaceComplex"]:
        PyMOLObjectNames["ComplexCryoEMGroupMembers"].append(CryoEMSurfaceName)
    
    PyMOLObjectNames["ComplexGroupMembers"] = []
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["Complex"])
    PyMOLObjectNames["ComplexGroupMembers"].append(PyMOLObjectNames["ComplexCryoEMGroup"])
    
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
    
    CryoEMMeshGroupName = "%s.CryoEM" % (ChainComplexGroupName)
    CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
    CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroup"] = CryoEMMeshGroupName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMMesh"] = CryoEMMeshName
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMSurface"] = CryoEMSurfaceName
    
    PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"] = []
    if MeshChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"].append(CryoEMMeshName)
    if SurfaceChainComplex:
        PyMOLObjectNames["Chains"][ChainID]["ChainComplexCryoEMGroupMembers"].append(CryoEMSurfaceName)
    
    NameIDs = ["ChainComplex"]
    if MeshChainComplex or SurfaceChainComplex :
        NameIDs.append("ChainComplexCryoEMGroup")
    
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
        
        CryoEMMeshGroupName = "%s.CryoEM" % (GroupName)
        CryoEMMeshName = "%s.Mesh" % (CryoEMMeshGroupName)
        CryoEMSurfaceName = "%s.Surface" % (CryoEMMeshGroupName)
                
        CryoEMMeshGroupID = "%sCryoEMMeshGroup" % (GroupID)
        CryoEMMeshGroupMembersID = "%sCryoEMMeshGroupMembers" % (GroupID)
        CryoEMMeshID = "%sCryoEMMesh" % (GroupID)
        CryoEMSurfaceID = "%sCryoEMSurface" % (GroupID)
        
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupID] = CryoEMMeshGroupName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshID] = CryoEMMeshName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMSurfaceID] = CryoEMSurfaceName
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID] = []
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID].append(CryoEMMeshName)
        PyMOLObjectNames["Ligands"][ChainID][LigandID][CryoEMMeshGroupMembersID].append(CryoEMSurfaceName)
                
        PyMOLObjectNames["Ligands"][ChainID][LigandID][GroupMembersID] = []
        NameIDs = [GroupTypeObjectID, CryoEMMeshGroupID]
        
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

def RetrieveDensityMapFileName():
    """Retrieve density map file name."""
    
    EMDBID = RetrieveEMDBID()
    if EMDBID is None:
        MiscUtil.PrintError("Failed to retrieve EMDB ID from input file %s to automatically set density map file name. Use option \"-d, --densityMapFile \" to specify density map file name and try again." % OptionsInfo["Infile"])

    DensityMapFile = None
    MapFileRoot = "emd_%s" % EMDBID
    MapFile1 = "%s.map.gz" % MapFileRoot
    MapFile2 = "%s.map" % MapFileRoot
    if os.path.exists(MapFile1):
        DensityMapFile = MapFile1
    elif os.path.exists(MapFile2):
        DensityMapFile = MapFile2
    else:
        MiscUtil.PrintError("Density map files %s or %s don't exist. Use option \"-d, --densityMapFile \" to specify density map file name and try again" % (MapFile1, MapFile2))
    
    MiscUtil.PrintInfo("Setting density map file name as  %s..." % DensityMapFile)
    
    return DensityMapFile

def RetrieveMeshLevel():
    """Retrieve recommened mesh contour level."""

    MeshLevel = None
    EMDBID = RetrieveEMDBID()
    if EMDBID is None:
        MiscUtil.PrintWarning("Failed to retrieve EMDB ID from input file %s to detect local header file already downloaded from EMDB server..." % OptionsInfo["Infile"])
        return MeshLevel

    MetadataHeaderFile = "emd-%s.xml" % (EMDBID)
    if not os.path.exists(MetadataHeaderFile):
        MiscUtil.PrintWarning("Failed to find a local header file, %s, for EMDB ID %s..." % (MetadataHeaderFile, EMDBID))
        return MeshLevel

    MiscUtil.PrintInfo("Retrieving recommeded mesh contour level from header file %s..." % MetadataHeaderFile)

    ContourLevel = None
    Source = None
    XMLTree = ElementTree.parse(MetadataHeaderFile)
    XMLRoot = XMLTree.getroot()

    MapElement = XMLTree.find("map")
    if MapElement is not None:
        ContourLevelElement = MapElement.find("contourLevel")
        if ContourLevelElement is not None:
            ContourLevel = ContourLevelElement.text
            Source = ContourLevelElement.get("source")

    if ContourLevel is not None:
        if Source is None:
            Source = "NA"
        MiscUtil.PrintInfo("Setting mesh level to recommended (Source: %s) mesh contour level value of %s..." % (Source, ContourLevel))
        MeshLevel = ContourLevel
    
    return MeshLevel

def RetrieveEMDBID():
    """Retrieve EMDB ID from input file. """

    if "EMDBID" in OptionsInfo:
        return OptionsInfo["EMDBID"]
    
    EMDBID = None
    
    Infile = OptionsInfo["Infile"]
    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)

    if re.match("^pdb$", FileExt, re.I):
        EMDBID = RetriveEMDBIDFromPDBFile(Infile)
    elif re.match("^cif$", FileExt, re.I):
        EMDBID = RetriveEMDBIDFromCIFFile(Infile)
    else:
        EMDBID = None

    OptionsInfo["EMDBID"] = EMDBID
    
    return EMDBID

def RetriveEMDBIDFromPDBFile(Infile):
    """Retrieve EMDB ID from PDB file. """

    EMDBID = None
    InfileFH = open(Infile, "r")
    if InfileFH is None:
        MiscUtil.PrintError("Couldn't open input file: %s.\n" % (Infile))

    MiscUtil.PrintInfo("\nRetrieving EMDB ID from input file %s..." % Infile)
    
    EMDBID = None
    for Line in InfileFH:
        Line = Line.rstrip()
        if re.match("^REMARK", Line, re.I):
            if re.search("DB: EMDB", Line, re.I):
                for Word in Line.split(" "):
                    # Retrieve string with EMD-
                    if re.search("EMD-", Word, re.I):
                        Word = Word.strip()
                        EMDBID = re.sub("EMD-", "", Word)
                        break
                break
    InfileFH.close()
        
    return EMDBID

def RetriveEMDBIDFromCIFFile(Infile):
    """Retrieve EMDB ID from CIF file. """

    InfileFH = open(Infile, "r")
    if InfileFH is None:
        MiscUtil.PrintError("Couldn't open input file: %s.\n" % (Infile))

    MiscUtil.PrintInfo("\nRetrieving EMDB ID from input file %s..." % Infile)
    
    EMDBID = None
    for Line in InfileFH:
        Line = Line.rstrip()
        if re.match("^EMDB  EMD", Line, re.I):
            for Word in Line.split(" "):
                # Retrieve string with EMD-
                if re.search("EMD-", Word, re.I):
                    Word = Word.strip()
                    EMDBID = re.sub("EMD-", "", Word)
                    break
            break
    InfileFH.close()
    
    return EMDBID

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

    OptionsInfo["LabelFontID"] = int(Options["--labelFontID"])
    
    # Process mesh parameters...
    OptionsInfo["MeshCarveRadius"] = float(Options["--meshCarveRadius"])
    OptionsInfo["MeshComplex"] = True if re.match("^Yes$", Options["--meshComplex"], re.I) else False
    OptionsInfo["MeshChainComplex"] = Options["--meshChainComplex"]
    
    OptionsInfo["MeshWidth"] = float(Options["--meshWidth"])
    OptionsInfo["MeshColor"] = Options["--meshColor"]
    
    OptionsInfo["SurfaceComplex"] = True if re.match("^Yes$", Options["--surfaceComplex"], re.I) else False
    OptionsInfo["SurfaceChainComplex"] = Options["--surfaceChainComplex"]
    OptionsInfo["SurfaceTransparency"] = float(Options["--surfaceTransparency"])
    
    OptionsInfo["PocketContactsLigandColor"] = Options["--pocketContactsLigandColor"]
    OptionsInfo["PocketContactsSolventColor"] = Options["--pocketContactsSolventColor"]
    OptionsInfo["PocketContactsInorganicColor"] = Options["--pocketContactsInorganicColor"]
    
    OptionsInfo["PocketDistanceCutoff"] = float(Options["--pocketDistanceCutoff"])
    OptionsInfo["PocketLabelColor"] = Options["--pocketLabelColor"]
    OptionsInfo["PocketSurface"] = True if re.match("^Yes$", Options["--pocketSurface"], re.I) else False

    DensityMapFile = Options["--densityMapFile"]
    if re.match("^auto$", DensityMapFile, re.I):
        DensityMapFile = RetrieveDensityMapFileName()
    OptionsInfo["DensityMapFile"] = Options["--densityMapFile"]
    OptionsInfo["DensityMapFileName"] = DensityMapFile
    
    MeshLevel = Options["--meshLevel"]
    if re.match("^auto$", MeshLevel, re.I):
        MeshLevel = RetrieveMeshLevel()
        if MeshLevel is None:
            MiscUtil.PrintWarning("Failed to retrieve recommended mesh contour level from header file. It's being set to 1.0. Use \"--meshLevel\" option to specify a different contour mesh level.")
            MeshLevel = 1.0
    OptionsInfo["MeshLevel"] = float(MeshLevel)

    # Process specified chains and ligands...
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

    if not re.match("^auto$", Options["--densityMapFile"], re.I):
        MiscUtil.ValidateOptionFilePath("-d, --densityMapFile", Options["--densityMapFile"])
        MiscUtil.ValidateOptionFileExt("-d, --densityMapFile", Options["--densityMapFile"], "map map.gz")
    
    MiscUtil.ValidateOptionFilePath("-i, --infile", Options["--infile"])
    MiscUtil.ValidateOptionFileExt("-i, --infile", Options["--infile"], "pdb cif")
    
    MiscUtil.ValidateOptionFileExt("-o, --outfile", Options["--outfile"], "pml pse")
    MiscUtil.ValidateOptionsOutputFileOverwrite("-o, --outfile", Options["--outfile"], "--overwrite", Options["--overwrite"])
    MiscUtil.ValidateOptionsDistinctFileNames("-i, --infile", Options["--infile"], "-o, --outfile", Options["--outfile"])
    
    MiscUtil.ValidateOptionIntegerValue("--labelFontID", Options["--labelFontID"], {})

    MiscUtil.ValidateOptionFloatValue("--meshCarveRadius", Options["--meshCarveRadius"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--meshComplex", Options["--meshComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--meshChainComplex", Options["--meshChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--meshWidth", Options["--meshWidth"], {">": 0.0})
    
    if not re.match("^auto$", Options["--meshLevel"], re.I):
        MiscUtil.ValidateOptionFloatValue("--meshLevel", Options["--meshLevel"], {})
    
    MiscUtil.ValidateOptionTextValue("--PMLOut", Options["--PMLOut"], "yes no")
    
    MiscUtil.ValidateOptionFloatValue("--pocketDistanceCutoff", Options["--pocketDistanceCutoff"], {">": 0.0})
    MiscUtil.ValidateOptionTextValue("--pocketSurface", Options["--pocketSurface"], "yes no")
    
    MiscUtil.ValidateOptionTextValue("--surfaceComplex", Options["--surfaceComplex"], "yes no")
    MiscUtil.ValidateOptionTextValue("--surfaceChainComplex", Options["--surfaceChainComplex"], "yes no auto")
    MiscUtil.ValidateOptionFloatValue("--surfaceTransparency", Options["--surfaceTransparency"], {">=": 0.0, "<=": 1.0})
    
# Setup a usage string for docopt...
_docoptUsage_ = """
PyMOLVisualizeCryoEMDensity.py - Visualize cryo-EM density

Usage:
    PyMOLVisualizeCryoEMDensity.py [--allowEmptyObjects <yes or no>]
                                   [--BFactorChainCartoonPutty <yes or no>] [--BFactorColorPalette <text> ]
                                   [--chainIDs <First, All or ID1,ID2...>] [--densityMapFile <text>]
                                   [--ligandIDs <Largest, All or ID1,ID2...>] [--labelFontID <number>]
                                   [--meshCarveRadius <number>] [--meshComplex <yes or no>]
                                   [--meshChainComplex <yes, no, or auto>] [--meshColor <text>]
                                   [--meshLevel <number>] [--meshWidth <number>] [--PMLOut <yes or no>]
                                   [--pocketContactsLigandColor <text>] [--pocketContactsSolventColor <text>]
                                   [--pocketContactsInorganicColor <text>] [--pocketDistanceCutoff <number>]
                                   [--pocketLabelColor <text>] [--pocketSurface <yes or no>]
                                   [--surfaceComplex <yes or no>] [--surfaceChainComplex <yes, no or auto>]
                                   [--surfaceTransparency <number>] [--overwrite] [-w <dir>] -i <infile> -o <outfile>
    PyMOLVisualizeCryoEMDensity.py -h | --help | -e | --examples

Description:
    Generate PyMOL visualization files for viewing electron microscopy (EM) or
    cryo-EM density around chains, ligands, and ligand binding pockets in
    macromolecules including proteins and nucleic acids.

    The supported input file formats are: Macromolecule - PDB (.pdb) or CIF(.cif),
    Cryo-EM Density - Collaborative Computational Project Number 4 (CCP4) ( .map)

    The supported output file formats are: PyMOL script file (.pml), PyMOL session
    file (.pse)

    The cryo-EM density and header files along with PDB files may be downloaded
    from appropriate servers using DownloadPDBFiles.pl script.

    A variety of PyMOL groups and objects may be  created for visualization of
    cryo-EM density present in map files. These groups and objects correspond to
    maps, meshes, surfaces,chains, ligands, inorganics, ligand binding pockets,
    pocket, polar interactions, and pocket hydrophobic surfaces. A complete
    hierarchy of all possible PyMOL groups and objects is shown below:
    
        <PDBFileRoot>
            .Complex
                .Complex
                .CryoEM
                    .Map
                    .Mesh
                    .Surface
            .Chain<ID>
                .Complex
                    .Complex
                    .CryoEM
                        .Mesh
                        .Surface
                .Chain
                    .Chain
                    .BFactor
                .Solvent
                .Inorganic
                .Ligand<ID>
                    .Ligand
                        .Ligand
                        .CryoEM
                            .Mesh
                            .Surface
                    .Pocket
                        .Pocket
                        .CryoEM
                            .Mesh
                            .Surface
                        .Polar_Contacts
                        .Surface
                    .Pocket_Solvent
                        .Pocket_Solvent
                        .CryoEM
                            .Mesh
                            .Surface
                        .Polar_Contacts
                    .Pocket_Inorganic
                        .Pocket_Inorganic
                        .CryoEM
                            .Mesh
                            .Surface
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
    of input complex and map files. The generation of PSE file may also fail. In 
    addition, you may want to interactively manipulate the contour level for meshes
    and surfaces. The recommended value for contour level is automatically retrieved
    from the header file available from EM density server. The recommended value
    may not always work.

Options:
    -a, --allowEmptyObjects <yes or no>  [default: no]
        Allow creation of empty PyMOL objects corresponding to solvent and
        inorganic atom selections across chains, ligands, and ligand binding pockets
        in input file.
    -b, --BFactorChainCartoonPutty <yes or no>  [default: yes]
        A cartoon putty around individual chains colored by B factors. The minimum
        and maximum values for B factors are automatically detected. These values
        indicate spread of cryo-EM density around atoms. The 'blue_white_red' color
        palette is deployed for coloring the cartoon putty.
    --BFactorColorPalette <text>  [default: blue_white_red]
        Color palette for coloring cartoon putty around chains generated using B
        factors. An valid PyMOL color palette name is allowed. No validation is
        performed. The complete list of valid color palette names is a available
        at: pymolwiki.org/index.php/Spectrum. Examples: blue_white_red,
        blue_white_magenta, blue_red, green_white_red, green_red.
    -c, --chainIDs <First, All or ID1,ID2...>  [default: First]
        List of chain IDs to use for visualizing cryo-EM density. Possible values:
        First, All, or a comma delimited list of chain IDs. The default is to use the
        chain ID for the first chain in input file.
    -d, --densityMapFile <text>  [default: auto]
        CryoEM density map file name. The EMDB ID is retrieved from PDB and CIF
        file to set the cryo-EM density file name during automatic detection of
        density file. The format of the file name is as follows:
        
            emd_<EMDBID>.map.gz or emd_<EMDBID>.map
         
        The density file must be present in the working directory. 
    -e, --examples
        Print examples.
    -h, --help
        Print this help message.
    -i, --infile <infile>
        Input file name.
    -l, --ligandIDs <Largest, All or ID1,ID2...>  [default: Largest]
        List of ligand IDs present in chains for visualizing cryo-EM density across
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
        Radius in Angstroms around atoms for including cryo-EM density.
    --meshComplex <yes or no>  [default: no]
        Create meshes for complete complex in input file corresponding to density
        map.
    --meshChainComplex <yes, no, or auto>  [default: auto]
        Create meshes for individual chain complex in input file corresponding to
        density map. By default, the meshes are automatically created for chain
        complexes without any ligands. 
    --meshColor <text>  [default: blue]
        Line color for mesh corresponding to density map. The specified value
        must be valid color. No validation is performed.
    --meshLevel <number>  [default: auto]
        Contour level in sigma units for generating mesh corresponding to density
        map. The default is to automatically retrieve the recommended contour
        level. The header file emd-<EMDBID>.xml must be present in the working
        directory  to automatically retrieve recommended value for mesh contour
        level. Otherwise, the default contour level is set to 1.
        
        You may want to interactively manipulate the contour level for meshes and
        surfaces. The default recommended value may not always work.
    --meshWidth <number>  [default: 0.5]
        Line width for mesh lines corresponding to density map.
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
        Create surfaces for complete complex in input file corresponding to density
        map.
    --surfaceChainComplex <yes, no or auto>  [default: auto]
        Create surfaces for individual chain complexes in input file corresponding to
        density map. By default, the surfaces are automatically created for chain complexes
        without any ligands.
    --surfaceTransparency <number>  [default: 0.25]
        Surface transparency for molecular and cryo-EM density surfaces.
    --overwrite
        Overwrite existing files.
    -w, --workingdir <dir>
        Location of working directory which defaults to the current directory.

Examples:
    To visualize cryo-EM density at recommended contour level for the first
    chain complex in a PDB file using corresponding density map and header
    file, and generate a PML file type:

        % PyMOLVisualizeCryoEMDensity.py -i 5K12.pdb -o 5K12.pml

    To visualize electron density for the largest ligand in  chain K, and ligand
    binding pocket to highlight ligand interactions with pockect residues,
    solvents and inorganics, in a PDB and using corresponding map files, and
    generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -c K -i 5UMD.cif -o 5UMD.pml

    To visualize electron density for all  chains along with any solvents in a
    PDB file and using corresponding map files, and generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -c all -i 5K12.pdb -o 5K12.pml

    To visualize cryo-EM density at a specific contour level for the first chain
    complex along with a mesh surface in a PDB file using corresponding to a
    specific density map file, and generate a PML file, type:

        % PyMOLVisualizeCryoEMDensity.py -d emd_8194.map.gz --meshLevel 1.0
          --surfaceChainComplex yes -i 5K12.pdb -o 5K12.pml

Author:
    Manish Sud(msud@san.rr.com)

See also:
    DownloadPDBFiles.pl,  PyMOLVisualizeElectronDensity.py,
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
