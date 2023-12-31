# --- UCSF Chimera Copyright ---
# Copyright (c) 2000-2006 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: Axes.py 28099 2009-07-14 20:53:13Z pett $

import chimera
from chimera import selection, triggers, SCENE_TOOL_SAVE, SCENE_TOOL_RESTORE
from prefs import prefs, AXIS_SEL_OBJ, AXIS_SEL_ATOMS, OBJ_SEL_AXIS, ATOMS_SEL_AXIS

class GeometryInterface:
	crossingPrompt = "Show angles/distances by choosing two objects from table"
	def __init__(self, parent, status, showCB):
		self.status = status
		self.parent = parent
		self.showCB = showCB
		import Pmw, Tkinter
		row = 0
		parent.columnconfigure(0, weight=1)
		parent.columnconfigure(1, weight=1)
		parent.columnconfigure(2, weight=1)
		Tkinter.Button(parent, text="Define axes...", pady=0,
			command=self._createAxesCB).grid(row=row, column=0)
		Tkinter.Button(parent, text="Define plane...", pady=0,
			command=self._createPlaneCB).grid(row=row, column=1)
		Tkinter.Button(parent, text="Define centroid...", pady=0,
			command=self._createCentroidCB).grid(row=row, column=2)
		row += 1

		from CGLtk.Table import SortableTable
		self.table = SortableTable(parent)
		self.table.grid(row=row, column=0, columnspan=3, sticky="nsew")
		parent.rowconfigure(row, weight=1)
		self.table.addColumn("Name", "name", format="%s ")
		self.table.addColumn(" ", "surfacePiece.color", format=(False, True))
		idCol = self.table.addColumn("ID", "id")
		self.table.addColumn("Shown", "model.display", format=bool)
		self.table.addColumn("Length", "length", format="%4.1f", font="TkFixedFont")
		self.table.addColumn("Radius", "radius", format="%4.1f", font="TkFixedFont")
		self.table.setData(geomManager.items)
		self.table.sortBy(idCol)
		self.table.launch(browseCmd=self._tableCB)
		row += 1

		selGroup = Pmw.Group(parent, tag_text="Choosing in table...")
		selGroup.grid(row=row, column=0, columnspan=3, sticky='ew')
		row += 1
		cbFrame = selGroup.interior()
		cbFrame.columnconfigure(0, weight=1)
		cbFrame.columnconfigure(1, weight=1)
		self.geomSelObjVar = Tkinter.IntVar(parent)
		self.geomSelObjVar.set(prefs[AXIS_SEL_OBJ])
		Tkinter.Checkbutton(cbFrame, variable=self.geomSelObjVar,
			text="selects object").grid(row=0, column=0)
		self.geomSelObjVar.trace_variable("w", self._selModeChangeCB)
		self.geomSelAtomsVar = Tkinter.IntVar(parent)
		self.geomSelAtomsVar.set(prefs[AXIS_SEL_ATOMS])
		Tkinter.Checkbutton(cbFrame, variable=self.geomSelAtomsVar,
			text="selects atoms").grid(row=0, column=1)
		self.geomSelAtomsVar.trace_variable("w", self._selModeChangeCB)

		self.buttons = []
		f = Tkinter.Frame(parent)
		f.grid(row=row, column=0, columnspan=3)
		row += 1
		self.buttons.append(Tkinter.Button(f, text="Delete", pady=0,
			state='disabled', command=self._deleteItemsCB))
		self.buttons[-1].grid(row=0, column=0, sticky='e')

		self.buttons.append(Tkinter.Button(f, text="Rename", pady=0,
			state='disabled', command=self._renameItems))
		self.buttons[-1].grid(row=0, column=1, sticky='w')
		self.renameDialogs = {}

		f = Tkinter.Frame(parent)
		f.grid(row=row, column=0, columnspan=3)
		row += 1
		self.buttons.append(Tkinter.Button(f, text="Report distance", pady=0,
			state='disabled', command=self._atomsDistanceCB))
		self.buttons[-1].grid(row=0, column=0)
		Tkinter.Label(f, text="to selected atoms"
			).grid(row=0, column=1)

		f = Tkinter.Frame(parent)
		f.grid(row=row, column=0, columnspan=3)
		row += 1
		self.alignButton = Tkinter.Button(f, text="Center and align", pady=0,
			state='disabled', command=self._alignCB)
		self.alignButton.grid(row=0, column=0)
		if chimera.tkgui.windowSystem == 'aqua':
			# way too wide otherwise...
			kw = { 'menubutton_width': 4 }
		else:
			kw = {}
		self.alignMenu = Pmw.OptionMenu(f, items=['X', 'Y', 'Z'],
			labelpos='w', label_text="chosen axis or plane normal along", **kw)
		self.alignMenu.grid(row=0, column=1)

		geomManager.registerInterface(self)
		if geomManager.items:
			self._selHandlerID = triggers.addHandler("selection changed",
													self._selChangeCB, None)
			self._modelHandlerID = triggers.addHandler("Model", self._modelChangeCB, None)

	def addItems(self, items):
		if len(items) == len(geomManager.items):
			self._selHandlerID = triggers.addHandler("selection changed",
													self._selChangeCB, None)
			self._modelHandlerID = triggers.addHandler("Model", self._modelChangeCB, None)
		self.table.setData(geomManager.items)
		if hasattr(self, 'createAxesDialog'):
			self.createAxesDialog.addItems(items)

	def feedback(self, sels):
		if not self.parent.winfo_viewable():
			return
		if len(sels) == 2:
			words1 = sels[0].name.split()
			words2 = sels[1].name.split()
			if len(words2) < len(words1):
				words1, words2 = words2, words1
			if len(words1) > 0:
				for i, word in enumerate(words1):
					if word != words2[i]:
						break
			else:
				i = 0
			diff = " ".join(words2[i:])
			infos = geomManager.infoTexts(sels)
			message = "%s/%s %s\n" % (" ".join(words1), diff, ", ".join(infos))
			self.status(message, log=True)
		else:
			self.status(self.crossingPrompt)

	def removeItems(self, items):
		self.table.setData(geomManager.items)
		for item in items:
			if item in self.renameDialogs:
				self.renameDialogs[item].destroy()
				del self.renameDialogs[item]
		if items and not geomManager.items:
			triggers.deleteHandler("selection changed", self._selHandlerID)
			triggers.deleteHandler("Model", self._modelHandlerID)
		if hasattr(self, 'createAxesDialog'):
			self.createAxesDialog.removeItems(items)

	def show(self):
		self.showCB()

	def _alignCB(self):
		# button is only active if exactly one item is selected
		# and that item is alignable
		item = self.table.selected()[0]
		frontPt, backPt = item.getAlignPoints()
		centerPt = frontPt
		frontPt = item.model.openState.xform.apply(frontPt)
		backPt = item.model.openState.xform.apply(backPt)
		xf = chimera.Xform.zAlign(backPt, frontPt)
		alignAxis = self.alignMenu.getvalue()
		if alignAxis == "X":
			xf.premultiply(chimera.Xform.yRotation(90))
		elif alignAxis == "Y":
			xf.premultiply(chimera.Xform.xRotation(90))
		someActive = False
		# make it work like the align command, which transforms
		# active models
		for molId in chimera.openModels.listIds():
			openState = chimera.openModels.openState(*molId)
			if openState.active:
				someActive = True
				openState.globalXform(xf)
		if not someActive:
			self.status("No models are active!", color="red")
			return
		if not item.model.openState.active:
			singular, plural = item.singularPlural()
			self.status("%s model is not active!" % singular.capitalize(),
				color="red")
		cp = item.model.openState.xform.apply(centerPt)
		chimera.viewer.camera.center = (cp.x, cp.y, cp.z)

	def _atomsDistanceCB(self):
		atoms = selection.currentAtoms()
		if not atoms:
			self.status("No atoms selected", color="red")
			return
		items = self.table.selected()
		if not items:
			self.status("No objects chosen", color="red")
			return
		for item in items:
			points = []
			dists = item.pointDistances([a.xformCoord()
							for a in atoms])
			if len(dists) == 1:
				self.status("Distance from %s to %s: %.1f" %
					(atoms[0], item, dists[0]), log=True)
				continue
			import numpy
			dists = numpy.array(dists)
			imin = dists.argmin()
			imax = dists.argmax()
			self.status("Distance from %d atoms to %s: "
				"min: %.1f (%s), mean: %.1f, max: %.1f (%s)\n"
				% (len(atoms), item, dists[imin], atoms[imin],
				dists.mean(), dists[imax], atoms[imax]),
				log=True)

	def _createAxesCB(self):
		if not hasattr(self, 'createAxesDialog'):
			from Axes import CreateAxesDialog
			self.createAxesDialog = CreateAxesDialog()
		self.createAxesDialog.enter()

	def _createCentroidCB(self):
		if not hasattr(self, 'createCentroidDialog'):
			from Centroids import CreateCentroidDialog
			self.createCentroidDialog = CreateCentroidDialog()
		self.createCentroidDialog.enter()

	def _createPlaneCB(self):
		if not hasattr(self, 'createPlaneDialog'):
			from Planes import CreatePlaneDialog
			self.createPlaneDialog = CreatePlaneDialog()
		self.createPlaneDialog.enter()

	def _deleteItemsCB(self):
		items = self.table.selected()
		if not items:
			self.status("No objects chosen", color="red")
			return
		geomManager.removeItems(items)

	def _itemsSelectContents(self, items, add=False):
		contents = set()
		for item in items:
			contents.update(item.selectables)
		if self.geomSelAtomsVar.get():
			if add:
				selection.addCurrent(contents)
			else:
				selection.setCurrent(contents)

	def _lowerCmd(self):
		for rd in self.renameDialogs.values():
			rd.destroy()
		self.renameDialogs.clear()
		self.status("")

	def _modelChangeCB(self, trigName, myData, trigData):
		changed = trigData.modified
		itemModels = set([item.model for item in geomManager.items])
		changedItems = changed & itemModels
		if changedItems:
			self.table.refresh()

	def _raiseCmd(self):
		self.status(self.crossingPrompt)

	def _renameItems(self):
		for item in self.table.selected():
			if item in self.renameDialogs:
				self.renameDialogs[item].enter()
			else:
				self.renameDialogs[item] = RenameDialog(self, item)

	def _renameItemCB(self, item, name):
		item.name = name
		self.table.refresh()
		self.renameDialogs[item].destroy()
		del self.renameDialogs[item]

	def _selChangeCB(self, *args):
		selGraphs = set(selection.currentGraphs())
		itemModels = set([item.model for item in geomManager.items])
		selItems = selGraphs & itemModels
		chosenItems = set([item.model for item in self.table.selected()])
		if chosenItems != selItems:
			items = [item for item in geomManager.items
								if item.model in selItems]
			self.table.select(items)
			if self.geomSelAtomsVar.get():
				self._itemsSelectContents(items, add=True)
			self.feedback(items)

	def _selModeChangeCB(self, *args):
		self._tableCB(self.table.selected(), modeChange=True)

	def _tableCB(self, sels, modeChange=False):
		if self.geomSelAtomsVar.get() or modeChange:
			self._itemsSelectContents(sels)
		# adding a Geometry to the selection while in the midst
		# of trigger processing is problematic, so...
		from chimera import update
		if self.geomSelObjVar.get() and not update.inTriggerProcessing:
			# avoid affecting the selection if models are
			# being destroyed...
			items = set(geomManager.items)
			removable = [a.model
						for a in items.difference(sels)
						if not a.model.__destroyed__]
			if removable:
				selection.removeCurrent(removable)
			addable = [a.model for a in sels if not a.model.__destroyed__]
			if addable:
				selection.addCurrent(addable)
		if modeChange:
			return
		butState = 'normal'
		if not sels:
			butState = 'disabled'
		for but in self.buttons:
			but.configure(state=butState)
		if len(sels) == 1 and sels[0].getAlignPoints()[0] is not None:
			self.alignButton.configure(state='normal')
		else:
			self.alignButton.configure(state='disabled')

		self.feedback(sels)

class Geometry(object):
	@classmethod
	def singularPlural(cls):
		return (getattr(cls, 'singular', cls.__name__).lower(),
			getattr(cls, 'plural', cls.__name__ + 's').lower())

	def __init__(self, number):
		self.number = number

	def getAlignPoints(self):
		# points used to align object with a coordinate axis
		return None, None

	@property
	def id(self):
		let = self.__class__.__name__[0].lower()
		from CGLutil.SortString import SortString
		return SortString(let + str(self.number), cmpVal=(let, self.number))

class GeometrySubmanager(object):
	restoreOrder = None # don't care (last)

	def _getRefInfo(self, items, massWeighting=False, modelOnly=False):
		lowest = items.lowestModel
		sourceModel = lowest
		if modelOnly:
			return sourceModel
		invXForm = lowest.openState.xform
		invXForm.invert()
		coords, weights = items.getCoordsWeights(lowest, invXForm, massWeighting)
		return sourceModel, coords, weights

	def _getSourceModel(self, atoms):
		sel = selection.ItemizedSelection()
		sel.add(atoms)
		mols = sel.molecules()
		if len(mols) == 1:
			return mols[0]
		return None

	def _restoreScene(self, sceneData):
		self.removeItems(self.items)
		if sceneData is not None:
			self._restoreSession(sceneData, scene=True)

	def _sceneData(self, *args, **kw):
		kw['scene'] = True
		return self._sessionData(*args, **kw)

class GeometryManager(object):
	def __init__(self):
		self.managers = []
		self.geomClasses = []
		self.interface = None
		self._suspendQueues = None
		from SimpleSession import SAVE_SESSION
		triggers.addHandler(SAVE_SESSION, self._sessionSave, None)
		triggers.addHandler(SCENE_TOOL_SAVE, self._sceneSave, None)
		triggers.addHandler(SCENE_TOOL_RESTORE, self._restoreScene, None)

	def addInterfaceItems(self, items):
		if self._suspendQueues:
			self._suspendQueues[0].extend(items)
		elif self.interface:
			self.interface.addItems(items)
	
	def angle(self, items):
		from Axes import Axis
		from Planes import Plane
		from Centroids import Centroid
		numAxes = len([i for i in items if isinstance(i, Axis)])
		numPlanes = len([i for i in items if isinstance(i, Plane)])
		numCentroids = len([i for i in items if isinstance(i, Centroid)])
		if numCentroids > 0:
			return None
		if numAxes == 2:
			axis1, axis2 = items
			angle = chimera.angle(axis1.xformDirection(),
									axis2.xformDirection())
			if angle > 90.0:
				angle = 180.0 - angle
		elif numAxes == numPlanes == 1:
			if isinstance(items[0], Axis):
				axis, plane = items
			else:
				plane, axis = items
			angle = chimera.angle(axis.xformDirection(), plane.xformNormal())
			angle = abs(90.0 - angle)
		elif numPlanes == 2:
			plane1, plane2 = items
			angle = chimera.angle(plane1.xformNormal(), plane2.xformNormal())
			if angle > 90.0:
				angle = 180.0 - angle
		else:
			raise ValueError("angle calculation not implemented")
		return angle

	def distance(self, items, **kw):
		from Axes import Axis
		from Planes import Plane
		from Centroids import Centroid
		numAxes = len([i for i in items if isinstance(i, Axis)])
		numPlanes = len([i for i in items if isinstance(i, Plane)])
		numCentroids = len([i for i in items if isinstance(i, Centroid)])
		if numAxes == 2:
			axis1, axis2 = items
			dist = axis1._axisDistance(axis2, **kw)
		elif numAxes == numPlanes == 1:
			if isinstance(items[0], Axis):
				axis, plane = items
			else:
				plane, axis = items
			xformPlane = chimera.Plane(plane.xformOrigin(), plane.xformNormal())
			xc = axis.xformCenter()
			d1, d2 = [xformPlane.distance(pt)
					for pt in [xc + axis.direction * e for e in axis.extents]]
			if d1 * d2 < 0.0:
				dist = 0.0
			else:
				dist = min([abs(d) for d in (d1, d2)])
		elif numPlanes == 2:
			plane1, plane2 = items
			angle = self.angle(items)
			if angle == 0.0:
				o1 = plane1.xformOrigin()
				p2 = chimera.Plane(plane2.xformOrigin(), plane2.xformNormal())
				dist = (o1 - p2.nearest(o1)).length
			else:
				dist = 0.0
		elif numCentroids > 0:
			if isinstance(items[0], Centroid):
				centroid, other = items
			else:
				other, centroid = items
			dist = other.pointDistances(centroid.xformCenter(), **kw)[0]
		else:
			raise ValueError("distance calculation not implemented")
		return dist

	def enableUpdates(self):
		if self._suspendQueues and self.interface:
			addQueue, removeQueue = self._suspendQueues
			if removeQueue:
				self.interface.removeItems(removeQueue)
			if addQueue:
				self.interface.addItems(addQueue)
		self._suspendQueues = None

	def infoTexts(self, items):
			from Axes import Axis
			texts = []
			angle = self.angle(items)
			if angle != None:
				if isinstance(items[0], Axis) and isinstance(items[1], Axis):
					description = "crossing angle"
				else:
					description = "angle"
				texts.append("%s: %.1f degrees" % (description, angle))
			dist = self.distance(items)
			if dist != None:
				texts.append("distance: %.1f angstroms" % dist)
			return texts

	@property
	def items(self):
		items = []
		for manager in self.managers:
			items.extend(manager.items)
		return items

	def registerManager(self, manager, managedClass):
		if not self.managers or manager.restoreOrder == None:
			self.managers.append(manager)
			self.geomClasses.append(managedClass)
		else:
			for i, mgr in enumerate(self.managers):
				if mgr.restoreOrder == None \
				or mgr.restoreOrder > manager.restoreOrder:
					self.managers = self.managers[:i] + [manager] + self.managers[i:]
					self.geomClasses = self.geomClasses[:i] + [managedClass] + self.geomClasses[i:]
					break
			else:
				self.managers.append(manager)
				self.geomClasses.append(managedClass)

	def removeItems(self, items):
		for manager, gclass in zip(self.managers, self.geomClasses):
			subitems = [i for i in items if isinstance(i, gclass)]
			if subitems:
				manager.removeItems(subitems)

	def removeInterfaceItems(self, geoms):
		if self._suspendQueues:
			self._suspendQueues[1].extend(geoms)
		elif self.interface:
			self.interface.removeItems(geoms)

	def registerInterface(self, interface):
		self.interface = interface

	@property
	def selectedItems(self):
		if not self.interface:
			from chimera.selection import currentGraphs
			selGraphs = currentGraphs(asDict=True)
			return [item for item in self.items if item in selGraphs]
		return self.interface.table.selected()

	def showInterface(self):
		if self.interface:
			self.interface.show()
		
	def suspendUpdates(self):
		if self._suspendQueues is None:
			self._suspendQueues = [[], []]

	def _restoreScene(self, trigName, myData, scene):
		self.suspendUpdates()
		info = scene.tool_settings.get("axes/planes/centroids", {})
		for manager in self.managers:
			manager._restoreScene(info.get(manager.__class__.__name__, None))
		self.enableUpdates()

	def _restoreSession(self, geomData):
		from SimpleSession import getColor, idLookup
		self.suspendUpdates()
		for manager in self.managers:
			try:
				data = geomData[manager.__class__.__name__]
			except KeyError:
				continue
			manager._restoreSession(data, fromGeom=True)
		self.enableUpdates()

	def _sceneSave(self, trigName, myData, scene):
		info = scene.tool_settings["axes/planes/centroids"] = {}
		info["version"] = 1
		for manager in self.managers:
			info[manager.__class__.__name__] = manager._sceneData()

	def _sessionSave(self, trigName, myData, sessionFile):
		from SimpleSession import sesRepr
		geomData = {}
		for manager in self.managers:
			geomData[manager.__class__.__name__] = manager._sessionData()
		print>>sessionFile, "geomData = %s" % sesRepr(geomData)
		print>>sessionFile, """
try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")
"""

geomManager = GeometryManager()

from chimera.baseDialog import ModelessDialog
class RenameDialog(ModelessDialog):
	buttons = ('OK', 'Cancel')
	default = 'OK'

	def __init__(self, geomUI, item):
		self.title = "Rename '%s'" % item.name
		self.geomUI = geomUI
		self.item = item
		ModelessDialog.__init__(self)

	def map(self, e=None):
		self.renameOpt._option.focus_set()

	def fillInUI(self, parent):
		from chimera.tkoptions import StringOption
		self.renameOpt = StringOption(parent, 0,
			"Rename %s to" % self.item.name, "", None)

	def Apply(self):
		newName = self.renameOpt.get().strip()
		if not newName:
			self.enter()
			from chimera import UserError
			raise UserError("Must supply a new item name or click Cancel")
		self.geomUI._renameItemCB(self.item, newName)

	def destroy(self):
		self.item = None
		self.geomUI = None
		ModelessDialog.destroy(self)
		
