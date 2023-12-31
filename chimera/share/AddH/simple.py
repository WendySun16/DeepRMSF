from chimera import Coord, Point
from chimera.bondGeom import tetrahedral, planar, linear, single, bondPositions
from chimera.idatm import *
from AddH import newHydrogen, findNearest, roomiest, bondWithHLength, \
							findRotamerNearest
from math import sin, cos, pi, sqrt
from chimera import Element
Element_H = Element(1)

sin5475 = sin(pi * 54.75 / 180.0)
cos5475 = cos(pi * 54.75 / 180.0)

def addHydrogens(atom, bondingInfo, namingSchema, totalHydrogens, idatmType,
							invert, coordinations, altLoc=None):
	away = away2 = planar = None
	geom = bondingInfo.geometry
	substs = bondingInfo.substituents
	if altLoc:
		curBonds = len([nb for nb in atom.neighbors if not nb.altLoc or nb.altLoc == altLoc])
	else:
		curBonds = len(atom.primaryBonds())
	needed = substs - curBonds
	if needed <= 0:
		return
	atPos = atom.xformCoord()
	exclude = coordinations + atom.neighbors
	if atom.altLoc:
		neighbors = [nb for nb in atom.neighbors if not nb.altLoc or nb.altLoc == atom.altLoc]
	elif altLoc:
		neighbors = [nb for nb in atom.neighbors if not nb.altLoc or nb.altLoc == altLoc]
	else:
		neighbors = atom.primaryNeighbors()
	occupancy = None
	if geom == 3:
		if curBonds == 1:
			if altLoc:
				bonded = [nb for nb in atom.neighbors
						if not nb.altLoc or nb.altLoc == altLoc][0]
				grandBonded = [nb for nb in bonded.neighbors
						if not nb.altLoc or nb.altLoc == altLoc]
				grandBonded.remove(atom)
				if len(grandBonded) < 3:
					planar = [a.xformCoord() for a in grandBonded]
				numOcc = 0
				totOcc = 0.0
				for a in grandBonded + [bonded]:
					if a.altLoc == altLoc:
						numOcc += 1
						totOcc += getattr(a, 'occupancy', 0.5)
				occupancy = totOcc / float(numOcc)
			elif len(atom.neighbors) > 1:
				for altLoc in set([nb.altLoc for nb in atom.neighbors if nb.altLoc]):
					addHydrogens(atom, bondingInfo, namingSchema, totalHydrogens, idatmType,
							invert, coordinations, altLoc=altLoc)
				return
			else:
				bonded = atom.neighbors[0]
				if atom.altLoc:
					grandBonded = [nb for nb in bonded.neighbors
							if not nb.altLoc or nb.altLoc == atom.altLoc]
				else:
					grandBonded = bonded.primaryNeighbors()
				grandBonded.remove(atom)
				if len(grandBonded) < 3:
					planar = [a.xformCoord() for a in grandBonded]
		elif curBonds == 2:
			if altLoc:
				bonded = [nb for nb in atom.neighbors
						if not nb.altLoc or nb.altLoc == altLoc]
				grandBonded = []
				for b in bonded:
					grandBonded.extend([nb for nb in b.neighbors
							if not nb.altLoc or nb.altLoc == altLoc])
				grandBonded.remove(atom)
				numOcc = 0
				totOcc = 0.0
				for a in grandBonded + bonded:
					if a.altLoc == altLoc:
						numOcc += 1
						totOcc += getattr(a, 'occupancy', 0.5)
				occupancy = totOcc / float(numOcc)
			elif len(atom.neighbors) > 2:
				for altLoc in set([nb.altLoc for nb in atom.neighbors if nb.altLoc]):
					addHydrogens(atom, bondingInfo, namingSchema, totalHydrogens, idatmType,
							invert, coordinations, altLoc=altLoc)
				return
	if geom == 4 and not neighbors:
		away, d, natom = findNearest(atPos, atom, exclude, 3.5)
		if away:
			away2, d2, natom2 = findRotamerNearest(atPos,
				idatmType[atom], atom, natom, 3.5)
	elif geom == 4 and len(neighbors + coordinations) == 1:
		away, d, natom = findRotamerNearest(atPos,
				idatmType[atom], atom, (neighbors+coordinations)[0], 3.5)
	else:
		away, d, natom = findNearest(atPos, atom, exclude, 3.5)

	bondedPos = []
	for bonded in neighbors:
		bondedPos.append(bonded.xformCoord())

	if coordinations:
		toward = coordinations[0].xformCoord()
		away2 = away
		away = None
	else:
		toward = None
	positions = bondPositions(atPos, geom, bondWithHLength(atom, geom),
		bondedPos, toward=toward, coPlanar=planar, away=away, away2=away2)
	if coordinations:
		coordPos = None
		for pos in positions:
			d = pos.sqdistance(toward)
			if coordPos is None or d < lowest:
				coordPos = pos
				lowest = d
		positions.remove(coordPos)
	if len(positions) > needed:
		positions = roomiest(positions, atom, 3.5, bondingInfo)[:needed]
	for i, pos in enumerate(positions):
		h = newHydrogen(atom, i+1, totalHydrogens, namingSchema,
							invert.apply(pos), bondingInfo, altLoc=altLoc)
		if occupancy is not None:
			h.occupancy = occupancy
