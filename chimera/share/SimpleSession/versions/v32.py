# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: v32.py 36625 2012-06-12 20:48:23Z pett $

from v31 import RemapDialog, reportRestoreError, restoreWindowSize, \
	restoreOpenModelsAttrs, noAutoRestore, autoRestorable, \
	registerAfterModelsCB, makeAfterModelsCBs, restoreModelClip, \
	restoreSelections, getColor, findFile, restorePseudoBondGroups, \
	sessionID, idLookup, init, restoreCamera, endRestore, \
	beginRestore, endRestore, restoreColors, restoreSurfaces, restoreVRML, \
	restoreOpenStates, restoreFontInfo, setSessionIDparams, \
	expandSequentialSummary, expandSummary

import globals # so that various version files can easily access same variables
import chimera

def restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo):
	items = []
	sm = globals.sessionMap

	res2mol = []
	atom2mol = []
	for ids, name, cid, display, lineWidth, pointSize, stickScale, \
	pdbHeaders, surfaceOpacity, ballScale, vdwDensity, autochain, \
	ribbonHidesMainchain, hidden in zip(
				expandSummary(molInfo['ids']),
				expandSummary(molInfo['name']),
				expandSummary(molInfo['color']),
				expandSummary(molInfo['display']),
				expandSummary(molInfo['lineWidth']),
				expandSummary(molInfo['pointSize']),
				expandSummary(molInfo['stickScale']),
				molInfo['pdbHeaders'],
				expandSummary(molInfo['surfaceOpacity']),
				expandSummary(molInfo['ballScale']),
				expandSummary(molInfo['vdwDensity']),
				expandSummary(molInfo['autochain']),
				expandSummary(molInfo['ribbonHidesMainchain']),
				expandSummary(molInfo['hidden'])
				):
		m = chimera.Molecule()
		sm[len(items)] = m
		items.append(m)
		m.name = name
		m.setAllPDBHeaders(pdbHeaders)
		from SimpleSession import modelMap, modelOffset
		chimera.openModels.add([m], baseId=ids[0]+modelOffset,
						subid=ids[1], hidden=hidden)
		modelMap.setdefault(ids, []).append(m)
		m.color = getColor(cid)
		m.display = display
		m.lineWidth = lineWidth
		m.pointSize = pointSize
		m.stickScale = stickScale
		m.surfaceOpacity = surfaceOpacity
		m.ballScale = ballScale
		m.vdwDensity = vdwDensity
		m.autochain = autochain
		m.ribbonHidesMainchain = ribbonHidesMainchain
	if molInfo['optional']:
		for attrName, info in molInfo['optional'].items():
			for m, val in zip(items, expandSummary(info)):
				if val is not None:
					setattr(m, attrName, val)

	resStart = len(items)
	for mid, name, chain, pos, insert, rcid, lcid, ss, ribbonDrawMode, \
	ribbonDisplay, label, isHet in zip(
				expandSummary(resInfo['molecule']),
				expandSummary(resInfo['name']),
				expandSummary(resInfo['chain']),
				expandSequentialSummary(resInfo['position']),
				expandSummary(resInfo['insert']),
				expandSummary(resInfo['ribbonColor']),
				expandSummary(resInfo['labelColor']),
				expandSummary(resInfo['ss']),
				expandSummary(resInfo['ribbonDrawMode']),
				expandSummary(resInfo['ribbonDisplay']),
				expandSummary(resInfo['label']),
				expandSummary(resInfo['isHet'])
				):
		m = idLookup(mid)
		r = m.newResidue(name, chain, pos, insert)
		sm[len(items)] = r
		items.append(r)
		r.ribbonColor = getColor(rcid)
		r.labelColor = getColor(lcid)
		r.isHelix, r.isStrand, r.isTurn = ss
		r.ribbonDrawMode = ribbonDrawMode
		r.ribbonDisplay = ribbonDisplay
		r.label = label
		r.isHet = isHet
	if resInfo['optional']:
		residues = items[resStart:]
		for attrName, info in resInfo['optional'].items():
			for r, val in zip(residues, expandSummary(info)):
				if val is not None:
					setattr(r, attrName, val)

	atomStart = len(items)
	for rid, name, element, cid, vcid, lcid, scid, drawMode, display, \
	label, surfaceDisplay, surfaceCategory, surfaceOpacity, radius, vdw, \
	idatmType in zip(
				expandSummary(atomInfo['residue']),
				expandSummary(atomInfo['name']),
				expandSummary(atomInfo['element']),
				expandSummary(atomInfo['color']),
				expandSummary(atomInfo['vdwColor']),
				expandSummary(atomInfo['labelColor']),
				expandSummary(atomInfo['surfaceColor']),
				expandSummary(atomInfo['drawMode']),
				expandSummary(atomInfo['display']),
				expandSummary(atomInfo['label']),
				expandSummary(atomInfo['surfaceDisplay']),
				expandSummary(atomInfo['surfaceCategory']),
				expandSummary(atomInfo['surfaceOpacity']),
				expandSummary(atomInfo['radius']),
				expandSummary(atomInfo['vdw']),
				expandSummary(atomInfo['idatmType'])
				):
		r = idLookup(rid)
		a = r.molecule.newAtom(name, chimera.Element(element))
		sm[len(items)] = a
		items.append(a)
		r.addAtom(a)
		a.color = getColor(cid)
		a.vdwColor = getColor(vcid)
		a.labelColor = getColor(lcid)
		a.surfaceColor = getColor(scid)
		# color inheritance changed starting with v56, so correct for that...
		if a.surfaceColor is None and a.color is not None and a.color != a.molecule.color:
			a.surfaceColor = a.color
		a.drawMode = drawMode
		a.display = display
		a.label = label
		a.surfaceDisplay = surfaceDisplay
		a.surfaceCategory = surfaceCategory
		a.surfaceOpacity = surfaceOpacity
		a.radius = radius
		a.vdw = vdw
		if idatmType:
			a.idatmType = idatmType
	if atomInfo['optional']:
		atoms = items[atomStart:]
		for attrName, info in atomInfo['optional'].items():
			for a, val in zip(atoms, expandSummary(info)):
				if val is not None:
					setattr(a, attrName, val)

	bondStart = len(items)
	for atoms, drawMode, display, radius in zip(
					bondInfo['atoms'],
					expandSummary(bondInfo['drawMode']),
					expandSummary(bondInfo['display']),
					expandSummary(bondInfo['radius'])
					):
		a1, a2 = [idLookup(a) for a in atoms]
		b = a1.molecule.newBond(a1, a2)
		sm[len(items)] = b
		items.append(b)
		b.drawMode = drawMode
		b.display = display
		b.radius = radius
	if bondInfo.get('optional', {}):
		bonds = items[bondStart:]
		for attrName, info in bondInfo['optional'].items():
			for b, val in zip(bonds, expandSummary(info)):
				if val is not None:
					setattr(b, attrName, val)

	from chimera import Point
	for mid, crdSets in crdInfo.items():
		m = idLookup(mid)
		active = crdSets.pop('active')
		for key, crds in crdSets.items():
			coordSet = m.newCoordSet(key, len(crds))
			for aid, crdString in crds:
				idLookup(aid).setCoord(Point(*tuple([float(c)
					for c in crdString.split()])), coordSet)
			if key == active:
				m.activeCoordSet = coordSet
