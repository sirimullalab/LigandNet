#
# File: PyMOLUtil.py
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
#
# The functionality available in this file is implemented using PyMOL, a
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

import os
import sys
import re

from pymol import cmd, stored, CmdException

import MiscUtil

__all__ = ["CalculateCenterOfMass", "ConvertFileFormat", "ConvertPMLFileToPSEFile", "GetChains", "GetChainsAndLigandsInfo", "GetInorganicResiduesInfo", "GetLargestLigand", "GetLigandResiduesInfo", "GetLigands", "GetMolecules", "GetPocketInorganicResiduesInfo", "GetPocketPolymerResiduesInfo", "GetPocketSolventResiduesInfo", "GetPolymerResiduesInfo", "GetSolventResiduesInfo", "ProcessChainsAndLigandsOptionsInfo", "SetupPMLForAlignment", "SetupPMLForBFactorPuttyView", "SetupPMLForBallAndStickView", "SetupPMLForEnableDisable", "SetupPMLForGroup", "SetupPMLForHydrophobicSurfaceView", "SetupPMLForInorganicView", "SetupPMLForLigandPocketInorganicView", "SetupPMLForLigandPocketSolventView", "SetupPMLForLigandPocketView", "SetupPMLForLigandView", "SetupPMLForPolarContactsView", "SetupPMLForPolymerChainComplexView", "SetupPMLForPolymerChainView", "SetupPMLForPolymerComplexView", "SetupPMLForSolventView", "SetupPMLForSurfaceView", "SetupPMLHeaderInfo"]

def GetMolecules(Selection = "all"):
    """Get names of molecule objects in a selection or all molecule objects.

    Arguments:
        Selection: (str): A PyMOL selection.

    Returns:
        list: Names of molecule objects.

    """
    Names = cmd.get_object_list('(' + Selection + ')')

    return Names

def GetChains(MoleculeName, RemoveEmpty = True):
    """Get chain identifiers present in a molecule.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        RemoveEmpty (bool): Remove empty chain ID from the list of chain IDs
            returned by PyMOL.

    Returns:
        list: Names of chains present in a molecule, sorted alphabetically in a
            ascending order.

    """
    if not len(MoleculeName):
        return None

    ChainIDs = []
    try:
        ChainIDs = cmd.get_chains('model %s' % MoleculeName)
    except CmdException as ErrMsg:
        MiscUtil.PrintWarning("PyMOLUtil.GetChains: Invalid molecule name: %s" % MoleculeName)

    if not len(ChainIDs):
        return None

    # Remove empty Chain IDs from the list...
    if RemoveEmpty:
        NonEmptyChainIDs = []
        for ChainID in ChainIDs:
            if len(ChainID):
                NonEmptyChainIDs.append(ChainID)
        if len(NonEmptyChainIDs) != len(ChainIDs):
            MiscUtil.PrintInfo("PyMOLUtil.GetChains: Removing non-empty chain IDs from the list of chain IDs...")
            
        ChainIDs = NonEmptyChainIDs
        
    return ChainIDs

def GetChainsAndLigandsInfo(Infile, MolName, Quite = False, LigandSortBy = "Size", LigandSortOrder = "Auto", LigandIgnoreHydrogens = "Yes"):
    """Get chain identifiers present in a molecule along with names of the
    ligands present in chains. Ligands are identified using PyMOL 'organic'
    selection.

    Arguments:
        Infile (str) : Name of a file.
        MolName (str) : Name to use for PyMOL molecule object.
        Quite (bool) : Flag 
        LigandSortBy (str): Sort ligand names alphabetically or by size. Possible
            values: Alphabetical or Size
        LigandSortOrder (str): Sort order for sorting ligands. Possible values:
            Ascending, Descending, Auto. The 'Auto' value implies automatic
            determination of sort order based on the value of 'SortBy'.
            Automatic defaults: Descending for SortBy value of Size; Ascending
            for SortBy value of Alphabetical.
        LigandIgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        dict: A dictionary containing list of chain identifiers and dictionaries
            of chains containing lists of ligand names for each chain. Names of
            ligands present in chain for a molecule sorted by size or
            alphabetically.

    Examples:

        ChainsAndLigandsInfo = GetChainsAndLigandsInfo(Infile, MolName)
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            for LigandID in ChainsAndLigandsInfo["LigandIDs"][ChainID]:
                MiscUtil.PrintInfo("ChainID: %s; LigandID: %s" % (ChainID,
                    LigandID))

    """
    if not Quite:
        MiscUtil.PrintInfo("\nRetrieving chain and ligand information from input file %s..." % Infile)
    
    # Collect chains and ligands information with ligands sorted by size to be used for
    # identification of largest ligand at the top...
    cmd.load(Infile, MolName)
    ChainsAndLigandsInfo = _GetChainsAndLigands(MolName, LigandSortBy = "size", LigandSortOrder = "descending")
    cmd.delete(MolName)

    # Print out chain and ligand IDs...
    if not Quite:
        ChainIDs = ", ".join(ChainsAndLigandsInfo["ChainIDs"]) if len(ChainsAndLigandsInfo["ChainIDs"]) else "None"
        MiscUtil.PrintInfo("Chain IDs: %s" % ChainIDs)
        
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            LigandIDs = ", ".join(ChainsAndLigandsInfo["LigandIDs"][ChainID]) if len(ChainsAndLigandsInfo["LigandIDs"][ChainID]) else "None"
            MiscUtil.PrintInfo("Chain ID: %s; LigandIDs: %s" % (ChainID, LigandIDs))
    
    return ChainsAndLigandsInfo

def GetLigands(MoleculeName, ChainName, SortBy = "Size", SortOrder = "Auto", IgnoreHydrogens = "Yes"):
    """Get names of ligands present in a chain of a  molecule. Ligands are
    identified using PyMOL 'organic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        SortBy (str): Sort ligand names alphabetically or by size. Possible
            values: Alphabetical or Size
        SortOrder (str): Sort order for sorting ligands. Possible values:
            Ascending, Descending, Auto. The 'Auto' value implies automatic
            determination of sort order based on the value of 'SortBy'.
            Automatic defaults: Descending for SortBy value of Size; Ascending
            for SortBy value of Alphabetical.
        IgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        list: Names of ligands present in chain for a molecule sorted by size
            or alphabetically.

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    LigandsInfoMap = _GetLigandsInfo(MoleculeName, ChainName, SortBy, SortOrder, IgnoreHydrogens)
    
    LigandIDs = LigandsInfoMap["LigandResNames"]
    if not len(LigandIDs):
        LigandIDs = None
    
    return LigandIDs

def GetLargestLigand(MoleculeName, ChainName, IgnoreHydrogens = 'Yes'):
    """Get name of the largest ligand for a chain present in a molecule. Ligands
    are identified using PyMOL 'organic' selection.

    Arguments:
        IgnoreHydrogens (str): Ignore hydrogens during determination of ligand
            size.

    Returns:
        str: Name of the largest ligand present in a chain. 

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SortBy = "Size"
    SortOrder = "Descending"
    LigandsInfoMap = _GetLigandsInfo(MoleculeName, ChainName, SortBy, SortOrder, IgnoreHydrogens)
    LigandIDs = LigandsInfoMap["LigandResNames"]
    
    if len(LigandIDs):
        LigandID = LigandIDs[0]
    else:
        LigandID = None

    return LigandID

def _GetChainsAndLigands(MoleculeName, LigandSortBy = "Size", LigandSortOrder = "Auto", LigandIgnoreHydrogens = "Yes"):
    """Get chain identifiers in a molecule along with names of the ligands
    present in chains. Ligands are identified using PyMOL 'organic' selection.
    """
    if not len(MoleculeName):
        return None

    ChainIDs = GetChains(MoleculeName)
    if ChainIDs is None:
        return None
        
    ChainsAndLigandsMap = {}
    ChainsAndLigandsMap["ChainIDs"] = []
    ChainsAndLigandsMap["LigandIDs"] = {}
    
    for ChainID in ChainIDs:
        ChainsAndLigandsMap["ChainIDs"].append(ChainID)
        ChainsAndLigandsMap["LigandIDs"][ChainID] = []
        
        LigandIDs = GetLigands(MoleculeName, ChainID, SortBy = LigandSortBy, SortOrder = LigandSortOrder, IgnoreHydrogens = LigandIgnoreHydrogens)
        if LigandIDs is not None:
            ChainsAndLigandsMap["LigandIDs"][ChainID] = LigandIDs
        
    return ChainsAndLigandsMap

def _GetLigandsInfo(MoleculeName, ChainName, SortBy = "Size", SortOrder = "Auto", IgnoreHydrogens = "Yes"):
    """Retrieve information about ligands present in a chain of a molecule."""

    if not MiscUtil.CheckTextValue(SortBy, "Size Alphabetical"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter SortBy is not valid. SupportedValues: Size Alphabetical" % SortBy)
        
    if not MiscUtil.CheckTextValue(SortOrder, "Ascending Descending Auto"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter SortOrder is not valid. SupportedValues: Ascending Descending Auto" % SortOrder)
        
    if not MiscUtil.CheckTextValue(IgnoreHydrogens, "Yes No"):
        MiscUtil.PrintError("PyMOLUtil._GetLigandsInfo: The value specified, %s, for parameter IgnoreHydrogens is not valid. SupportedValues: Yes No" % IgnoreHydrogens)
        
    SortBySize = True if re.match("^Size$", SortBy, re.I) else False
    if re.match("^Auto$", SortOrder, re.I):
        SortOrderDescending = True if re.match("^Size$", SortBy, re.I) else False
    else:
        SortOrderDescending = True if re.match("^Descending$", SortOrder, re.I) else False
    IgnoreHydrogenAtoms = True if re.match("^Yes$", IgnoreHydrogens, re.I) else False

    # Set up a command to retrieve all appropriate ligand atoms in organic ligands...
    SelectionCmd = "%s and chain %s and organic" % (MoleculeName, ChainName)
    if IgnoreHydrogenAtoms:
        SelectionCmd = "%s and not hydro" % (SelectionCmd)

    # Retrieve atoms...
    stored.LigandsInfo = []
    cmd.iterate(SelectionCmd, "stored.LigandsInfo.append([resi, resn])")

    # Retrieve ligands...
    LigandsInfoMap = {}
    LigandsInfoMap["LigandResNames"] = []
    LigandsInfoMap["LigandAtomCount"] = {}
    LigandsInfoMap["LigandResNumber"] = {}
    
    for LigandResNum, LigandResName in stored.LigandsInfo:
        if LigandResName in LigandsInfoMap["LigandResNames"]:
            LigandsInfoMap["LigandAtomCount"][LigandResName] += 1
        else:
            LigandsInfoMap["LigandResNames"].append(LigandResName)
            LigandsInfoMap["LigandAtomCount"][LigandResName] = 1
            LigandsInfoMap["LigandResNumber"][LigandResName] = LigandResNum

    if not len(LigandsInfoMap["LigandResNames"]):
        return LigandsInfoMap
        
    # Sort ligand names...
    ReverseOrder = True if SortOrderDescending else False
    if SortBySize:
        SortedLigandResNames = sorted(LigandsInfoMap["LigandResNames"], key = lambda LigandResName: LigandsInfoMap["LigandAtomCount"][LigandResName], reverse = ReverseOrder)
    else:
        # Sort alphabetically...
        SortedLigandResNames = sorted(LigandsInfoMap["LigandResNames"], reverse = ReverseOrder)

    LigandsInfoMap["LigandResNames"] = SortedLigandResNames
    
    return LigandsInfoMap

def GetPolymerResiduesInfo(MoleculeName, ChainName):
    """Get information for residues present in a chain of a  molecule.
    Chains are identified using PyMOL 'polymer' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and polymer" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetSolventResiduesInfo(MoleculeName, ChainName):
    """Get information for solvent residues present in a chain of a  molecule.
    Solvents are identified using PyMOL 'solvent' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetSolventResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and solvent" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetInorganicResiduesInfo(MoleculeName, ChainName):
    """Get information for inorganic residues present in a chain of a  molecule.
    Inorganic residues are identified using PyMOL 'inorganic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetInorganicResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and inorganic" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetLigandResiduesInfo(MoleculeName, ChainName):
    """Get information for ligand residues present in a chain of a  molecule.
    Ligands are identified using PyMOL 'organic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetLigandResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName)):
        return None

    SelectionCmd = "%s and chain %s and organic" % (MoleculeName, ChainName)
    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketPolymerResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for chain residues present in a pocket around a ligand
    in a molecule. Polymer residues are identified using negation of PyMOL
    selection operators 'organic', 'solvent', and 'inorganic'.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around ligand to identify pocket
            residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketPolymerResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and (not solvent) and (not inorganic) and (not organic))" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketSolventResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for solvent residues present in a pocket around a ligand
    in a molecule. Solvent residues are identified using PyMOL 'solvent'
    selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around ligand to identify pocket
            residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketSolventResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and solvent)" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def GetPocketInorganicResiduesInfo(MoleculeName, ChainName, LigandResName, LigandResNum, PocketDistanceCutoff):
    """Get information for inorganic residues present in a pocket around a
    ligand in a molecule. Inorganic residues are identified using PyMOL
    'inorganic' selection.

    Arguments:
        MoleculeName (str): Name of a PyMOL molecule object.
        ChainName (str): Name of a chain in a molecule.
        LigandResName (str): Residue name of a ligand in a chain.
        LigandResNum (str): Residue number of a ligand in a chain.
        PocketDistanceCutoff (float): Distance around a ligand to identify
            pocket residues.

    Returns:
        dict: A dictionary containing list of residue names and dictionaries of
            residue numbers and residue count for each residue. Names of 
            residues in the dictionary are not sorted.

    Examples:

        ResiduesInfo = GetPocketInorganicResiduesInfo(MolName, ChainName)
        for ResName in ResiduesInfo["ResNames"]:
            ResCount = ResiduesInfo["ResCount"][ResName]
            ResNums = ResiduesInfo["ResNum"][ResName]
            MiscUtil.PrintInfo("ResName: %s; ResCount: %s; ResNums: %s" %
                (ResName, ResCount, ResNums))

    """
    if not (len(MoleculeName) and len(ChainName) and len(LigandResName) and len(LigandResNum)):
        return None

    LigandSelection = "%s and chain %s and organic and resn %s and resi %s" % (MoleculeName, ChainName, LigandResName, LigandResNum)
    MoleculeSelection = "%s and chain %s" % (MoleculeName, ChainName)
    SelectionCmd = "((byresidue (%s) within %.1f of (%s))  and inorganic)" % (MoleculeSelection, PocketDistanceCutoff, LigandSelection)

    ResiduesInfoMap = _GetSelectionResiduesInfo(SelectionCmd)

    return ResiduesInfoMap

def _GetSelectionResiduesInfo(SelectionCmd):
    """Get residue names and count information for a selection. """
    
    # Retrieve atoms...
    stored.SelectionInfo = []
    cmd.iterate(SelectionCmd, "stored.SelectionInfo.append([resi, resn])")

    # Retrieve residues...
    SelectionInfoMap = {}
    SelectionInfoMap["ResNames"] = []
    SelectionInfoMap["ResNum"] = {}
    SelectionInfoMap["ResCount"] = {}

    PreviousResNum = 0
    for ResNum, ResName in stored.SelectionInfo:
        if ResName in SelectionInfoMap["ResNames"]:
            if ResNum != PreviousResNum:
                # Same residue name but different residue number
                SelectionInfoMap["ResNum"][ResName].append(ResNum)
                PreviousResNum = ResNum
                
                SelectionInfoMap["ResCount"][ResName] += 1
        else:
            SelectionInfoMap["ResNames"].append(ResName)
            
            SelectionInfoMap["ResNum"][ResName] = []
            SelectionInfoMap["ResNum"][ResName].append(ResNum)
            PreviousResNum = ResNum
            
            SelectionInfoMap["ResCount"][ResName] = 1
                
    return SelectionInfoMap

def ProcessChainsAndLigandsOptionsInfo(ChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue, LigandsOptionName = None, LigandsOptionValue = None):
    """Process specified chain and ligand IDs using command line options.

    Arguments:
        ChainsAndLigandsInfo (dict): A dictionary containing information
            existing chains and ligands. 
        ChainsOptionName (str): Name of command line chains option.
        ChainsOptionValue (str): Value for command line chains option.
        LigandsOptionName (str): Name of command line ligands option.
        LigandsOptionValue (str): Value for command line ligands option.

    Returns:
        dict: A dictionary containing list of chain identifiers and dictionaries
            of chains containing lists of ligand names for each chain.

    Examples:

        ChainsAndLigandsInfo = ProcessChainsAndLigandsOptionsInfo(Infile,
            MolName)
        for ChainID in ChainsAndLigandsInfo["ChainIDs"]:
            for LigandID in ChainsAndLigandsInfo["LigandIDs"][ChainID]:
                MiscUtil.PrintInfo("ChainID: %s; LigandID: %s" % (ChainID,
                    LigandID))

    """
    SpecifiedChainsAndLigandsInfo = {}
    SpecifiedChainsAndLigandsInfo["ChainIDs"] = []
    SpecifiedChainsAndLigandsInfo["LigandIDs"] = {}

    if ChainsOptionValue is None:
        return SpecifiedChainsAndLigandsInfo
        
    _ProcessChainIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue)

    if LigandsOptionValue is None:
        return SpecifiedChainsAndLigandsInfo
    
    _ProcessLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, LigandsOptionName, LigandsOptionValue)
    
    return SpecifiedChainsAndLigandsInfo

def _ProcessChainIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, ChainsOptionName, ChainsOptionValue):
    """Process chain IDs"""

    MiscUtil.PrintInfo("Processing chain IDs...")
    
    if re.match("^All$", ChainsOptionValue, re.I):
        SpecifiedChainsAndLigandsInfo["ChainIDs"] = ChainsAndLigandsInfo["ChainIDs"]
        return
    elif re.match("^(First|Auto)$", ChainsOptionValue, re.I):
        FirstChainID = ChainsAndLigandsInfo["ChainIDs"][0] if (len(ChainsAndLigandsInfo["ChainIDs"])) else None
        if FirstChainID is not None:
            SpecifiedChainsAndLigandsInfo["ChainIDs"].append(FirstChainID)
        return
    
    ChainIDs = re.sub(" ", "", ChainsOptionValue)
    if not ChainIDs:
        MiscUtil.PrintError("No valid value specified using \"%s\" option." % ChainsOptionName)

    ChainIDsList = ChainsAndLigandsInfo["ChainIDs"]
    SpecifiedChainIDsList = []
    
    ChainIDsWords = ChainIDs.split(",")
    for ChainID in ChainIDsWords:
        if not ChainID in ChainIDsList:
            MiscUtil.PrintWarning("The chain ID, %s, specified using \"%s\" option is not valid. It'll be ignored. Valid chain IDs: %s" % (ChainID, ChainsOptionName, ", ".join(ChainIDsList)))
            continue
        if ChainID in SpecifiedChainIDsList:
            MiscUtil.PrintWarning("The chain ID, %s, has already been specified using \"%s\" option. It'll be ignored." % (ChainID, ChainsOptionName))
            continue
        SpecifiedChainIDsList.append(ChainID)
    
    if not len(SpecifiedChainIDsList):
        MiscUtil.PrintError("No valid chain IDs \"%s\"  specified using \"%s\" option." % (ChainsOptionValue, ChainsOptionName))
    
    SpecifiedChainsAndLigandsInfo["ChainIDs"] = SpecifiedChainIDsList
    
def _ProcessLigandIDs(ChainsAndLigandsInfo, SpecifiedChainsAndLigandsInfo, LigandsOptionName, LigandsOptionValue):
    """Process ligand IDs"""

    MiscUtil.PrintInfo("Processing ligand IDs...")
    
    # Intialize ligand IDs...
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"] :
        SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = []
        
    if re.match("^All$", LigandsOptionValue, re.I):
        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"] :
            SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = ChainsAndLigandsInfo["LigandIDs"][ChainID]
        return
    elif re.match("^(Largest|Auto)$", LigandsOptionValue, re.I):
        # Setup largest ligand ID for each chain...
        for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"] :
            LargestLigandID = ChainsAndLigandsInfo["LigandIDs"][ChainID][0] if (len(ChainsAndLigandsInfo["LigandIDs"][ChainID])) else None
            if LargestLigandID is not None:
                SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID].append(LargestLigandID)
        return
    
    LigandIDs = re.sub(" ", "", LigandsOptionValue)
    if not LigandIDs:
        MiscUtil.PrintError("No valid value specified using \"%s\" option." % LigandsOptionName)
    
    LigandIDsWords = LigandIDs.split(",")
    
    for ChainID in SpecifiedChainsAndLigandsInfo["ChainIDs"]:
        LigandIDsList = ChainsAndLigandsInfo["LigandIDs"][ChainID]
        SpecifiedLigandIDsList = []

        for LigandID in LigandIDsWords:
            if not LigandID in LigandIDsList:
                LigandIDsListNames = ",".join(LigandIDsList) if len(LigandIDsList) else "None"
                MiscUtil.PrintWarning("The ligand ID, %s, specified using \"%s\" option is not valid for chain, %s. It'll be ignored. Valid ligand IDs: %s" % (LigandID, LigandsOptionName, ChainID, LigandIDsListNames))
                continue
            if LigandID in SpecifiedLigandIDsList:
                MiscUtil.PrintWarning("The ligand ID, %s, has already been specified using \"%s\" option. It'll be ignored." % (LigandID, LigandsOptionName))
                continue
            SpecifiedLigandIDsList.append(LigandID)
            
        if not len(SpecifiedLigandIDsList):
            MiscUtil.PrintWarning("No valid ligand IDs \"%s\" specified using \"%s\" option for chain ID, %s." % (LigandsOptionValue, LigandsOptionName, ChainID))
        
        SpecifiedChainsAndLigandsInfo["LigandIDs"][ChainID] = SpecifiedLigandIDsList
    
def CalculateCenterOfMass(Selection = "all", Quiet = 0):
    """Calculate center of mass for a selection.

    Arguments:
        Selection (str): A PyMOL selection.
        Quiet (int): Print information.

    Returns:
        list: X, Y, Z coordinates for center of mass.

    """
    MassTotal = 0.0
    X, Y, Z = [0.0, 0.0, 0.0]
    
    Atoms = cmd.get_model(Selection)
    for Atom in Atoms.atom:
        Mass = Atom.get_mass()
        MassTotal += Mass

        X += Atom.coord[0] * Mass
        Y += Atom.coord[1] * Mass
        Z += Atom.coord[2] * Mass

    XCOM = X/MassTotal
    YCOM = Y/MassTotal
    ZCOM = Z/MassTotal

    if not Quiet:
        MiscUtil.PrintInfo("PyMOLUtil.CalculateCenterOfMass: %f, %f, %f" % (XCOM, YCOM, ZCOM))
    
    return [XCOM, YCOM, ZCOM]

def ConvertFileFormat(Infile, Outfile, Reinitialize = True, OutputFeedback = True):
    """Convert infile to outfile by automatically detecting their formats
    from the file extensions.

    The formats of both input and output files must be a valid format supported
    by PyMOL.

    Arguments:
        Infile (str): Name of input file.
        Outfile (str): Name of outfile file.
        Reinitialize (bool): Reinitialize PyMOL before loading input file.
        OutputFeedback (bool): Control output feedback.

    """
    
    if not os.path.exists(Infile):
        MiscUtil.PrintWarning("The input file, %s, doesn't exists.%s..." % (Infile))

    if Reinitialize:
        cmd.reinitialize()

    if not OutputFeedback:
        # Turn off output feedback...
        MiscUtil.PrintInfo("Disabling output feedback for PyMOL...")
        cmd.feedback("disable", "all", "output")

    FileDir, FileName, FileExt = MiscUtil.ParseFileName(Infile)
    MolName = FileName
    
    cmd.load(Infile, MolName)
    cmd.save(Outfile, MolName)
    cmd.delete(MolName)
    
    if not OutputFeedback:
        # Turn it back on...
        MiscUtil.PrintInfo("\nEnabling output feedback for PyMOL...")
        cmd.feedback("enable", "all", "output")

def ConvertPMLFileToPSEFile(PMLFile, PSEFile, Reinitialize = True, OutputFeedback = True):
    """Convert PML file to PME file.

    Arguments:
        PMLFile (str): Name of PML file.
        PSEFile (str): Name of PSE file.
        Reinitialize (bool): Reinitialize PyMOL before loading PML file.
        OutputFeedback (bool): Control output feedback.

    """
    
    if not os.path.exists(PMLFile):
        MiscUtil.PrintWarning("The PML file, %s, doesn't exists.%s..." % (PMLFile))

    if Reinitialize:
        cmd.reinitialize()

    if not OutputFeedback:
        # Turn off output feedback...
        MiscUtil.PrintInfo("Disabling output feedback for PyMOL...")
        cmd.feedback("disable", "all", "output")

    cmd.do("@%s" % PMLFile)
    cmd.save(PSEFile)
    
    if not OutputFeedback:
        # Turn it back on...
        MiscUtil.PrintInfo("\nEnabling output feedback for PyMOL...")
        cmd.feedback("enable", "all", "output")

def SetupPMLHeaderInfo(ScriptName = None, IncludeLocalPython = True):
    """Setup header information for generating PML files. The local Python
    functions are optionally embedded in the header information for their
    usage in PML files.

    Arguments:
        ScriptName (str): Name of script calling the function.
        IncludeLocalPython (bool): Include local Python functions.

    Returns:
        str: Text containing header information for generating PML files.

    """
    if ScriptName is None:
        HeaderInfo = """\
#
# This file is automatically generated by a script available in MayaChemTools.
#
cmd.reinitialize()"""
    else:
        HeaderInfo = """\
#
# This file is automatically generated by the following PyMOL script available in
# MayaChemTools: %s
#
cmd.reinitialize() """ % (ScriptName)

    if IncludeLocalPython:
        PMLForLocalPython = _SetupPMLForLocalPython()
        HeaderInfo = "%s\n\n%s" % (HeaderInfo, PMLForLocalPython)
    
    return HeaderInfo

def _SetupPMLForLocalPython():
    """Setup local Python functions for PML file.
    """
    
    PMLForPython = """\
""
"Setting up local Python functions  for PML script..."
""
python

from __future__ import print_function
import re

def ColorByHydrophobicity(Selection, ColorPalette = "RedToWhite"):
        \"""Color by hydrophobicity using hydrophobic values for amino acids
        corresponding to the Eisenberg hydrophobicity scale.

        Possible values for ColorPalette: RedToWhite or WhiteToGreen from most
        hydrophobic amino acid to least hydrophobic.

        The colors values for amino acids are taken from color_h script avaiable
        as part of the Script Library at PyMOL Wiki. 
        
        \"""
        
        if not re.match("^(RedToWhite|WhiteToGreen)$", ColorPalette, re.I):
                print("Invalid ColorPalette value: %s. Valid values: RedToWhite, WhiteToGreen" % ColorPalette)

        ResColors = {}
        ColorType = ""
        ColorsAvailable = True
        
        if re.match("^WhiteToGreen$", ColorPalette, re.I):
                ColorType = "H2"
                ResColors = {"ile" : [0.938,1,0.938], "phe" : [0.891,1,0.891], "val" : [0.844,1,0.844], "leu" : [0.793,1,0.793], "trp" : [0.746,1,0.746], "met" : [0.699,1,0.699], "ala" : [0.652,1,0.652], "gly" : [0.606,1,0.606], "cys" : [0.555,1,0.555], "tyr" : [0.508,1,0.508], "pro" : [0.461,1,0.461], "thr" : [0.414,1,0.414], "ser" : [0.363,1,0.363], "his" : [0.316,1,0.316], "glu" : [0.27,1,0.27], "asn" : [0.223,1,0.223], "gln" : [0.176,1,0.176], "asp" : [0.125,1,0.125], "lys" : [0.078,1,0.078], "arg" : [0.031,1,0.031]}
        else:
                ColorType = "H1"
                ResColors = {"ile" : [0.996,0.062,0.062], "phe" : [0.996,0.109,0.109], "val" : [0.992,0.156,0.156], "leu" : [0.992,0.207,0.207], "trp" : [0.992,0.254,0.254], "met" : [0.988,0.301,0.301], "ala" : [0.988,0.348,0.348], "gly" : [0.984,0.394,0.394], "cys" : [0.984,0.445,0.445], "tyr" : [0.984,0.492,0.492], "pro" : [0.980,0.539,0.539], "thr" : [0.980,0.586,0.586], "ser" : [0.980,0.637,0.637], "his" : [0.977,0.684,0.684], "glu" : [0.977,0.730,0.730], "asn" : [0.973,0.777,0.777], "gln" : [0.973,0.824,0.824], "asp" : [0.973,0.875,0.875], "lys" : [0.899,0.922,0.922], "arg" : [0.899,0.969,0.969]}
        
        # Set up colors...
        for ResName in ResColors:
                ColorName = "color_%s_%s" % (ResName, ColorType)
                cmd.set_color(ColorName, ResColors[ResName])

                ResSelection = "(%s and resn %s*)" % (Selection, ResName)
                cmd.color(ColorName, ResSelection)
        
cmd.extend("ColorByHydrophobicity", ColorByHydrophobicity)

def CheckAndDeleteEmptyObject(ObjectName, ParentObjectName = None):
    \"""Delete an empty object along with optionally deleting its parent.
        
    \"""
    
    if  cmd.count_atoms("(%s)" % ObjectName):
        return

    cmd.delete("%s" % ObjectName)
    
    if ParentObjectName is not None:
        cmd.delete("%s" % ParentObjectName)
    
cmd.extend("CheckAndDeleteEmptyObject", CheckAndDeleteEmptyObject)

python end"""
    
    return PMLForPython

def SetupPMLForEnableDisable(Name, Enable = True):
    """Setup PML command for enabling or disabling display of a PyMOL object.

    Arguments:
        Name (str): Name of a PyMOL object.
        Enable (bool): Display status.

    Returns:
        str: PML command for enabling or disabling display of an object.

    """
    
    if Enable:
        PML = """cmd.enable("%s")""" % Name
    else:
        PML = """cmd.disable("%s")""" % Name
        
    return PML

def SetupPMLForGroup(GroupName, GroupMembersList, Enable = None, Action = None):
    """Setup PML commands for creating a group from a list of group members. The
    display and open status of the group may be optionally set. The 'None' values
    for Enable and Action imply usage of PyMOL defaults for the creation of group.

    Arguments:
        GroupName (str): Name of a PyMOL group.
        GroupMembersList (list): List of group member names.
        Enable (bool): Display status of group.
        Action (str): Open or close status of group object.

    Returns:
        str: PML commands for creating a group object.

    """

    PMLCmds = []
    
    GroupMembers = " ".join(GroupMembersList)
    PMLCmds.append("""cmd.group("%s", "%s")""" % (GroupName, GroupMembers))
    
    if Enable is not None:
        if Enable:
            PMLCmds.append("""cmd.enable("%s")""" % GroupName)
        else:
            PMLCmds.append("""cmd.disable("%s")""" % GroupName)
    
    if Action is not None:
        PMLCmds.append("""cmd.group("%s", action="%s")""" % (GroupName, Action))

    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandView(Name, Selection, LigandResName, Enable = True):
    """Setup PML commands for creating a ligand view corresponding to a ligand 
    present in a selection. The ligand is identified using organic selection
    operator available in PyMOL in conjunction with the specified ligand ID.
    The ligand is colored by atom types and displayed as 'sticks'.

    Arguments:
        Name (str): Name of a new PyMOL ligand object.
        Selection (str): PyMOL selection containing ligand.
        LigandResName (str): Ligand ID.
        Enable (bool): Display status of ligand object.

    Returns:
        str: PML commands for a ligand view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and organic and (resn %s))")""" % (Name, Selection, LigandResName))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""util.cbag("%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding all residues present in a selection within a specified
    distance from a ligand selection. The solvent and inorganic portions of
    the selection are not included in the binding pocket. The pocket residues
    are shown as 'lines'. The hydrogen atoms are not displayed.

    Arguments:
        Name (str): Name of a new PyMOL binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pockect residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and (not solvent) and (not inorganic) and (not organic))")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "(%s)")""" % (Name))
    PMLCmds.append("""cmd.hide("(%s and hydro)")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketSolventView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding to only solvent residues present in a selection within a
    specified distance from a ligand selection. The solvent pocket residues
    are shown as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pocket solvent residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view only showing solvent
            residues.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and solvent)")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForLigandPocketInorganicView(Name, Selection, LigandSelection, DistanceCutoff, Enable = True):
    """Setup PML commands for creating a ligand binding pocket view
    corresponding to only inorganic residues present in a selection within a
    specified distance from a ligand selection. The inorganic pocket residues
    are shown as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent binding pocket object.
        Selection (str): PyMOL selection containing binding pocket residues.
        LigandSelection (str): PyMOL selection containing ligand.
        DistanceCutoff (float): Distance cutoff from ligand for selecting
            binding pocket inorganic residues.
        Enable (bool): Display status of binding pocket object.

    Returns:
        str: PML commands for a ligand binding pocket view only showing inorganic
            residues.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((byresidue (%s) within %.1f of (%s)) and inorganic)")""" % (Name, Selection, DistanceCutoff, LigandSelection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolarContactsView(Name, Selection1, Selection2, Enable = True, Color = "yellow"):
    """Setup PML commands for creating polar contacts view between a pair of
    selections. The polar contact view is generated using 'util.dist' command. The
    distance labels are shown by default.

    Arguments:
        Name (str): Name of a new PyMOL polar contacts object.
        Selection1 (str): First PyMOL selection.
        Selection2 (str): Second PyMOL selection.
        Enable (bool): Display status of polar contacts object.
        Colot (str): Color for polar contact lines and labels.

    Returns:
        str: PML commands for polar contacts view between a pair of selections.

    """
    
    PMLCmds = []
    PMLCmds.append("""cmd.dist("%s","(%s)","(%s)",quiet = 1, mode = 2, label = 1, reset = 1)""" % (Name, Selection1, Selection2))
    PMLCmds.append("""util.color_deep("%s", "%s")""" % (Color, Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForAlignment(Method, RefSelection, FitSelection):
    """Setup PML commands for aligning a pair of selection using  a specified
    alignment method.

    Arguments:
        Method (str): Alignment method. Possible values: align, cealign, super.
        RefSelection (str): Name of reference selection which stays stationary.
        FitSelection (str): Name of selection to align to reference selection.

    Returns:
        str: PML commands for aligning  a pair of selections.

    """

    PMLCmds = []
    if re.match("^align$", Method, re.I):
        PMLCmds.append("""cmd.align("(%s)", "(%s)")""" % (FitSelection, RefSelection))
    elif re.match("^cealign$", Method, re.I):
        PMLCmds.append("""cmd.cealign("(%s)", "(%s)")""" % (RefSelection, FitSelection))
    elif re.match("^super$", Method, re.I):
        PMLCmds.append("""cmd.super("(%s)", "(%s)")""" % (FitSelection, RefSelection))
    else:
        MiscUtil.PrintWarning("PyMOLUtil.SetupPMLForAlignment: Invalid method name: %s" % Method)

    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForBFactorPuttyView(Name, Selection, ColorPalette = "blue_white_red", Enable = True):
    """Setup PML commands for creating a B factor putty view for a specified
    selection. The B factor values must be available for the atoms. The atoms
    are colored using a color spectrum corresponding to a specified color
    palette. Any valid PyMOL color palette name may be used.

    Arguments:
        Name (str): Name of a new PyMOL B factor putty object.
        Selection (str): Name of PyMOL selection.
        ColorPalette (str): Name of color palette to use for color spectrum.
        Enable (bool): Display status of B factor putty object.

    Returns:
        str: PML commands for B factor putty view.

    """
    
    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.spectrum("b", "%s", "(%s)")""" % (ColorPalette, Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append("""cmd.cartoon("putty", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForHydrophobicSurfaceView(Name, Selection, ColorPalette = "RedToWhite", Enable = True):
    """Setup PML commands for creating a hydrophobic surface view for a specified
    selection. The surfaces are colored using a specified color palette. This is only valid
    for amino acids.

    Arguments:
        Name (str): Name of a new PyMOL hydrophobic surface object.
        Selection (str): Name of PyMOL selection.
        ColorPalette (str): Name of color palette to use for coloring surfaces.
            Possible values: RedToWhite or WhiteToGreen for most hydrophobic
            to least hydrophobic amino acids.
        Enable (bool): Display status of surface object.

    Returns:
        str: PML commands for hydrophobic surface view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("surface", "%s")""" % (Name))
    PMLCmds.append("""ColorByHydrophobicity("%s", "%s")""" % (Name, ColorPalette))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForSurfaceView(Name, Selection, Enable = True):
    """Setup PML commands for creating a molecular surface view for a specified
    selection.

    Arguments:
        Name (str): Name of a new PyMOL molecular surface object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of surface object.

    Returns:
        str: PML commands for molecular surface view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("surface", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForBallAndStickView(Name, Selection, Enable = True, SphereScale = 0.3, StickRadius = 0.2):
    """Setup PML commands for creating a ball and stick view for a specified
    selection.

    Arguments:
        Name (str): Name of a new PyMOL ball and stick object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of ball and stick object.
        SphereScale (float): Scaling factor for sphere radii.
        StickScale (float): Scaling factor for stick radii.

    Returns:
        str: PML commands for ball and stick view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "(%s)")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("sphere", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "%s")""" % (Name))
    PMLCmds.append("""cmd.set("sphere_scale", %.1f, "%s")""" % (SphereScale, Name))
    PMLCmds.append("""cmd.set("stick_radius", %.1f, "%s")""" % (StickRadius, Name))
    
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForInorganicView(Name, Selection, Enable = True):
    """Setup PML commands for creating a inorganic view corresponding to
    inorganic residues present in a selection. The inorganic residues are
    identified using inorganic selection operator available in PyMOL. The
    inorganic residues are displayed as 'lines' and 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL inorganic object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of inorganic object.

    Returns:
        str: PML commands for inorganic view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and inorganic)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForSolventView(Name, Selection, Enable = True):
    """Setup PML commands for creating a solvent view corresponding to
    solvent residues present in a selection. The solvent residues are
    identified using solvent selection operator available in PyMOL. The
    solvent residues are displayed as 'nonbonded'.

    Arguments:
        Name (str): Name of a new PyMOL solvent object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of inorganic object.

    Returns:
        str: PML commands for solvent view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and solvent)")""" % (Name, Selection))
    PMLCmds.append("""cmd.show("nonbonded", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolymerChainView(Name, Selection, Enable = True):
    """Setup PML commands for creating a polymer chain view corresponding
    to backbone and sidechain residues in a selection. The polymer chain is
    displayed as 'cartoon'.

    Arguments:
        Name (str): Name of a new PyMOL polymer chain object.
        Selection (str): Name of PyMOL selection.
        Enable (bool): Display status of chain object.

    Returns:
        str: PML commands for polymer chain view.

    """

    PMLCmds = []
    PMLCmds.append("""cmd.create("%s", "((%s) and (backbone or sidechain))")""" % (Name, Selection))
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""util.cbag("%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))
    
    
    PML = "\n".join(PMLCmds)
    
    return PML


def SetupPMLForPolymerComplexView(MoleculeName, PDBFile, Enable = True, ShowSolvent = True, ShowInorganic = True, ShowLines = True):
    """Setup PML commands for creating a polymer complex view for all chains
    in a PDB file. The solvent and inorganic residues are also shown by default.
     The polymer chains are displayed as 'cartoon'. The 'line' display for the
    polymer chains is also shown and may be turned off. The organic residues are
    displayed as 'sticks'. The solvent and inorganic residues are displayed as
    'nonbonded' and 'lines'.
    
    Arguments:
        MoleculeName (str): Name of a new PyMOL molecule object.
        PDBFile (str): Name of PDB file.
        Enable (bool): Display status of chain object.
        ShowSolvent (bool): Display solvent residues.
        ShowInorganic (bool): Display inorganic residues.
        ShowLines (bool): Display lines for polymer chains.

    Returns:
        str: PML commands for polymer complex view.

    """

    PMLCmds = []
    
    PMLCmds.append("""cmd.load("%s", "%s")""" % (PDBFile, MoleculeName))
    PML = _SetupPMLForPolymerComplexView(MoleculeName, Enable, ShowSolvent, ShowInorganic, ShowLines)
    PMLCmds.append(PML)
    
    PML = "\n".join(PMLCmds)
    
    return PML

def SetupPMLForPolymerChainComplexView(ChainComplexName, Selection, ChainName, Enable = True, ShowSolvent = True, ShowInorganic = True, ShowLines = True):
    """Setup PML commands for creating a polymer chain complex view for a specified
    chain in a selection. The solvent and inorganic residues are also shown by
    default. The polymer chain is displayed as 'cartoon'. The 'line' display for the
    polymer chain is also shown and may be turned off. The organic residues are
    displayed as 'sticks'. The solvent and inorganic residues are displayed as
    'nonbonded' and 'lines'.

    Arguments:
        ChainComplexName (str): Name of a new PyMOL polymer chain complex.
        Selection (str): Name of PyMOL selection.
        ChainName (str): Name of a chain.
        Enable (bool): Display status of chain object.
        ShowSolvent (bool): Display solvent residues.
        ShowInorganic (bool): Display inorganic residues.
        ShowLines (bool): Display lines for polymer chain.

    Returns:
        str: PML commands for polymer chain complex view.

    """

    PMLCmds = []
    
    PMLCmds.append("""cmd.create("%s", "(%s and chain %s)")""" % (ChainComplexName, Selection, ChainName))
    PML = _SetupPMLForPolymerComplexView(ChainComplexName, Enable, ShowSolvent, ShowInorganic, ShowLines)
    PMLCmds.append(PML)

    PML = "\n".join(PMLCmds)
    
    return PML

def _SetupPMLForPolymerComplexView(Name, Enable = True, ShowSolvent = True, ShowInorganic = True,  ShowLines = False):
    """Setup PML for creating a polymer complex view."""

    PMLCmds = []
    
    PMLCmds.append("""cmd.hide("everything", "%s")""" % (Name))
    PMLCmds.append("""cmd.show("cartoon", "%s")""" % (Name))
    PMLCmds.append("""util.cba(33, "%s", _self = cmd)""" % (Name))
    PMLCmds.append("""cmd.show("sticks", "(organic and (%s))")""" % (Name))
    if ShowSolvent:
        PMLCmds.append("""cmd.show("nonbonded", "(solvent and (%s))")""" % (Name))
    if ShowInorganic:
        PMLCmds.append("""cmd.show("nonbonded", "(inorganic and (%s))")""" % (Name))
    
    if ShowLines:
        PMLCmds.append("""cmd.show("lines", "%s")""" % (Name))
    else:
        if ShowInorganic:
            PMLCmds.append("""cmd.show("lines", "(inorganic and (%s))")""" % (Name))
    
    PMLCmds.append("""cmd.set_bond("valence", "1", "%s", quiet = 1)""" % (Name))
    PMLCmds.append(SetupPMLForEnableDisable(Name, Enable))

    PML = "\n".join(PMLCmds)
    
    return PML
