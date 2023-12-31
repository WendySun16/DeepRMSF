# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: gui.py 38962 2013-07-23 23:04:58Z pett $

import chimera
from chimera.baseDialog import ModelessDialog
from chimera.selection import ItemizedSelection
from chimera.tkoptions import ColorOption, LineWidthOption
from chimera import replyobj, help
import Tkinter
from base import findHBonds, recDistSlop, recAngleSlop
import Pmw
from copy import deepcopy
from prefs import prefs, BOND_COLOR, RELAX_COLOR, LINE_WIDTH

ui = None

def showUI(callback=None):
	global ui

	if not ui:
		ui = HBDialog()
	ui.enter()
	if callback:
		ui.addCallback(callback)

class HBDialog(ModelessDialog):
	title = "H-Bond Parameters"
	buttons = ("OK", "Apply", "Close")
	default = 'OK'
	help = "ContributedSoftware/findhbond/findhbond.html"

	MDResNetMode = False

	SelShorthands = ["any", "cross", "both", "osl"]

	def fillInUI(self, parent):
		# inter/intra-model; slop; colors; models; replace current;
		# label

		# for removing labels in subsequent passes...
		# (use a selection so I don't have to track model closures)
		self.labelSelection = ItemizedSelection()
		self.labelMap = {}

		self.callbacks = []

		row = 0
		cols = 3

		if not self.MDResNetMode:
			self.bondColorOption = ColorOption(parent, row, 'H-bond color',
						prefs[BOND_COLOR], None, noneOkay=0)
			class HBLineWidthOption(LineWidthOption):
				name = "Line width"
				balloon = "Width of pseudobonds (in pixels)"
				default = prefs[LINE_WIDTH]
			self.lineWidth = HBLineWidthOption(parent, row+1,
						None, None, None, width=4)
			df = Tkinter.Frame(parent)
			df.grid(row=row+2, column=0, columnspan=2)
			self.showDistVar = Tkinter.IntVar(parent)
			self.showDistVar.set(False)
			Tkinter.Checkbutton(df, variable=self.showDistVar, text="Label H-bond"
				" with distance", command=self._distShowCB).grid(
				row=0, column=0, columnspan=2)
			from chimera.tkgui import windowSystem
			if windowSystem == 'aqua':
				kw = {}
			else:
				kw = {'padx': 0, 'pady': 0}
			self.dfWidgets = []
			self.dfWidgets.append(Tkinter.Button(df, text="Show",
				command=self._showDistFormatting, **kw))
			self.dfWidgets[-1].grid(row=1, column=0, sticky='e')
			self.dfWidgets.append(Tkinter.Label(df,
				text="distance formatting options"))
			self.dfWidgets[-1].grid(row=1, column=1, sticky='w')
			for dw in self.dfWidgets:
				dw.grid_remove()

			self.findMode = Pmw.RadioSelect(parent,
				buttontype='radiobutton', labelpos='n', pady=0,
				orient='vertical', label_text='Find these bonds:',
				hull_borderwidth=2, hull_relief='ridge')
			self.findMode.add('inter-model')
			self.findMode.add('intra-model')
			self.findMode.add('both')
			self.findMode.invoke('both')
			self.findMode.grid(row=row, column=2, rowspan=3)
			row += 3

		self.relaxParams = RelaxParams(parent, prefs[RELAX_COLOR],
			colorOptions=not self.MDResNetMode)
		self.relaxParams.grid(row=row, columnspan=cols,
							sticky='nsew', padx=2)
		row += 1

		if not self.MDResNetMode:
			self.restrictModels = Tkinter.IntVar(parent)
			self.restrictModels.set(False)
			self.modelFrame = Tkinter.Frame(parent)
			self.molTexts = ["Restrict to models...", "Restrict to models:"]
			self.restrictCheckBox = Tkinter.Checkbutton(self.modelFrame,
				variable=self.restrictModels,
				text=self.molTexts[0], command=self._restrictModelsCB)
			self.restrictCheckBox.grid(row=0, column=0)
			from chimera.widgets import MoleculeScrolledListBox
			self.molList = MoleculeScrolledListBox(self.modelFrame,
							listbox_selectmode="extended")
			self.modelFrame.rowconfigure(0, weight=1)
			self.modelFrame.grid(row=row, column=0, columnspan=cols,
									sticky="nsw")
			parent.rowconfigure(row, weight=1)
			row += 1

		self.inSelection = Tkinter.IntVar(parent)
		self.inSelection.set(False)
		f = Tkinter.Frame(parent)
		f.grid(row=row, column=0, columnspan=cols, sticky="ew")
		Tkinter.Checkbutton(f, variable=self.inSelection,
					text="Only find H-bonds",
					command=self._useSelChange).grid(
					row=0, column=0, sticky="e")
		self.selType = Pmw.OptionMenu(f, items=[
			"with at least one end selected", "with exactly one end selected",
			"with both ends selected", "between selection and atom spec..."],
			initialitem=0, menubutton_state="disabled", command=self._mapEntry)
		self.selType.grid(row=0, column=1, sticky="w")
		self.STAtomSpec = Tkinter.Entry(f, width=10)
		f.columnconfigure(2, weight=1)
		for c in range(cols):
			parent.columnconfigure(c, weight=1)
		row += 1

		if self.MDResNetMode:
			return
		self.includeIntraMol = Tkinter.IntVar(parent)
		self.includeIntraMol.set(True)
		Tkinter.Checkbutton(parent, variable=self.includeIntraMol,
			text="Include intra-molecule H-bonds"
			).grid(row=row, column=0, columnspan=cols, sticky="w")
		row += 1

		self.includeIntraRes = Tkinter.IntVar(parent)
		self.includeIntraRes.set(True)
		Tkinter.Checkbutton(parent, variable=self.includeIntraRes,
			text="Include intra-residue H-bonds"
			).grid(row=row, column=0, columnspan=cols, sticky="w")
		row += 1

		self.revealEnds = Tkinter.IntVar(parent)
		self.revealEnds.set(False)
		Tkinter.Checkbutton(parent, variable=self.revealEnds,
			text="If endpoint atom hidden, show endpoint residue"
			).grid(row=row, column=0, columnspan=cols, sticky="w")
		row += 1

		self.retainCurrent = Tkinter.IntVar(parent)
		self.retainCurrent.set(False)
		Tkinter.Checkbutton(parent, variable=self.retainCurrent,
			text="Retain currently displayed H-bonds").grid(
			row=row, column=0, columnspan=cols, sticky="w")
		row += 1

		self.interSubmodel = Tkinter.IntVar(parent)
		self.interSubmodel.set(False)
		# comment out the below since there seems to be no
		# situations where cross-submodel-same-model hbonds
		# are interesting.  If reading pdbrun files occurs, this
		# may change.
		#Tkinter.Checkbutton(parent, variable=self.interSubmodel,
		#	text="Find H-bonds between submodels\n"
		#	     "having same principal model number"
		#	).grid(row=row, column=0, columnspan=cols, sticky="w")
		#row += 1
		
		self.writeFile = Tkinter.IntVar(parent)
		self.writeFile.set(False)
		Tkinter.Checkbutton(parent, variable=self.writeFile,
			text="Write information to file").grid(row=row,
			column=0, columnspan=cols, sticky="w")
		self.fileName = None
		row += 1

		self.writeLog = Tkinter.IntVar(parent)
		self.writeLog.set(False)
		Tkinter.Checkbutton(parent, variable=self.writeLog,
			text="Write information to reply log").grid(row=row,
			column=0, columnspan=cols, sticky="w")
		row += 1

	def addCallback(self, callback):
		self.callbacks.append(callback)

	def Apply(self):
		states = {}
		for but in self.buttonWidgets.values():
			states[but] = but.cget('state')
			but.config(state='disabled')
		try:
			self._Apply()
		finally:
			for but, state in states.items():
				but.config(state=state)

		for cb in self.callbacks:
			try:
				cb()
			except:
				replyobj.reportException(
						"H-bond callback function")
		self.callbacks = []

	def _Apply(self):
		intramodel = intermodel = 1
		cursel = self.findMode.getcurselection()
		if cursel == 'inter-model':
			intramodel = 0
		elif cursel == 'intra-model':
			intermodel = 0
		
		distSlop = 0.0
		angleSlop = 0.0
		twoColors = False
		relax = self.relaxParams.relaxConstraints
		rc = self.relaxParams.relaxColor
		if relax:
			distSlop = self.relaxParams.relaxDist
			angleSlop = self.relaxParams.relaxAngle
			if self.relaxParams.useRelaxColor:
				twoColors = True
				prefs[RELAX_COLOR] = rc.rgba()
		
		if self.inSelection.get():
			selRestrict = self.SelShorthands[self.selType.index(Pmw.SELECT)]
			if selRestrict == "osl":
				selRestrict = self.STAtomSpec.get()
		else:
			selRestrict = None

		lineWidth = self.lineWidth.get()
		prefs[LINE_WIDTH] = lineWidth
		bc = self.bondColorOption.get()
		prefs[BOND_COLOR] = bc.rgba()

		if self.writeFile.get():
			saveFile = '-'
		else:
			saveFile = None

		if self.restrictModels.get():
			models = self.molList.getvalue()
			if not models:
				raise ValueError, "No restriction models chosen"
		else:
			models = None

		from base import createHBonds
		from Midas import MidasError
		try:
			createHBonds(models=models, intramodel=intramodel,
				intermodel=intermodel, relax=relax, distSlop=distSlop,
				angleSlop=angleSlop, twoColors=twoColors,
				selRestrict=selRestrict, lineWidth=lineWidth,
				saveFile=saveFile, log=self.writeLog.get(),
				retainCurrent=self.retainCurrent.get(),
				intraMol=self.includeIntraMol.get(),
				intraRes=self.includeIntraRes.get(),
				reveal=self.revealEnds.get(), color=bc, slopColor=rc,
				showDist=self.showDistVar.get())
		except MidasError, v:
			from chimera import UserError
			raise UserError(v)
	
	def Close(self):
		self.callbacks = []
		ModelessDialog.Close(self)

	def _distShowCB(self):
		if self.showDistVar.get():
			for dw in self.dfWidgets:
				dw.grid()
		else:
			for dw in self.dfWidgets:
				dw.grid_remove()

	def _mapEntry(self, val):
		import Pmw
		if self.selType.index(val) == self.selType.index(Pmw.END):
			self.STAtomSpec.grid(row=0, column=2, sticky="ew")
		else:
			self.STAtomSpec.grid_forget()

	def _restrictModelsCB(self):
		restrict = self.restrictModels.get()
		self.restrictCheckBox.config(text=self.molTexts[restrict])
		if restrict:
			self.molList.grid(row=0, column=1, sticky="ns")
		else:
			self.molList.grid_forget()

	def _showDistFormatting(self):
		from StructMeasure.gui import DISTANCES, StructMeasure
		from chimera.dialogs import display
		dlg = display(StructMeasure.name)
		dlg.notebook.selectpage(DISTANCES)

	def _useSelChange(self):
		if self.inSelection.get():
			self.selType.component("menubutton").config(
							state="normal")
		else:
			self.selType.component("menubutton").config(
							state="disabled")

class RelaxParams(Pmw.Group):
	def __init__(self, parent, relaxColor, colorOptions=True):
		self._relaxConstraints = Tkinter.IntVar(parent)
		self._relaxConstraints.set(True)
		Pmw.Group.__init__(self, parent, ring_relief='ridge',
					tag_pyclass=Tkinter.Checkbutton,
					tag_text='Relax H-bond constraints',
					tag_variable=self._relaxConstraints)
		self._relaxLabel = Tkinter.Label(self.interior(),
			text="Relax constraints by: ")
		self._relaxLabel.grid(row=0, column=0, rowspan=2, sticky='e')
		self._relaxDist = Pmw.EntryField(self.interior(),
			labelpos='e', label_text='angstroms', validate={
			'validator': 'real', 'min': 0.0 },
			value=str(recDistSlop), entry_width=5)
		self._relaxDist.grid(row=0, column=1, sticky='w')
		self._relaxAngle = Pmw.EntryField(self.interior(),
			labelpos='e', label_text='degrees', validate={
			'validator': 'real', 'min': 0.0 },
			value=str(recAngleSlop), entry_width=5)
		self._relaxAngle.grid(row=1, column=1, sticky='w')
		if not colorOptions:
			return
		self._useRelaxColor = Tkinter.IntVar(parent)
		self._useRelaxColor.set(False)
		f = Tkinter.Frame(self.interior())
		f.grid(row=2, column=0, columnspan=2)
		Tkinter.Checkbutton(f, text="Color H-bonds not meeting "
			"precise criteria differently:",
			variable=self._useRelaxColor).grid(row=0,
							column=0, sticky='e')
		self._relaxColor = ColorOption(f, 0, "", relaxColor, None,
							startCol=1, noneOkay=0)

	def enable(self):
		self.configure(tag_state=Tkinter.NORMAL)
		self._relaxLabel.config(state=Tkinter.NORMAL)
		#self._relaxDist.config(entry_state=Tkinter.NORMAL,
		#		label_state=Tkinter.NORMAL)
		#self._relaxAngle.config(entry_state=Tkinter.NORMAL,
		#		label_state=Tkinter.NORMAL)
		# not sure why the above doesn't work, so:
		for w in (self._relaxDist, self._relaxAngle):
			for c in ('entry', 'label'):
				wc = w.component(c)
				wc.config(state=Tkinter.NORMAL)
		# TODO: color

	def disable(self):
		self.configure(tag_state=Tkinter.DISABLED)
		self._relaxLabel.config(state=Tkinter.DISABLED)
		#self._relaxDist.config(entry_state=Tkinter.DISABLED,
		#		label_state=Tkinter.DISABLED)
		#self._relaxAngle.config(entry_state=Tkinter.DISABLED,
		#		label_state=Tkinter.DISABLED)
		# not sure why the above doesn't work, so:
		for w in (self._relaxDist, self._relaxAngle):
			for c in ('entry', 'label'):
				wc = w.component(c)
				wc.config(state=Tkinter.DISABLED)
		# TODO: color

	def __getattr__(self, attrName):
		if attrName == "relaxConstraints":
			return self._relaxConstraints.get()
		if attrName == "relaxDist":
			self._relaxDist.invoke()
			if not self._relaxDist.valid():
				from chimera import UserError
				raise UserError('Invalid "relaxed" distance'
					' value; must be non-negative number')
			return float(self._relaxDist.getvalue())
		if attrName == "relaxAngle":
			self._relaxAngle.invoke()
			if not self._relaxAngle.valid():
				from chimera import UserError
				raise UserError('Invalid "relaxed" angle'
					' value; must be non-negative number')
			return float(self._relaxAngle.getvalue())
		if attrName == "useRelaxColor":
			return self._useRelaxColor.get()
		if attrName == "relaxColor":
			return self._relaxColor.get()
		return Pmw.Group.__getattr__(self, attrName)

	def __setattr__(self, attrName, value):
		if attrName == "relaxConstraints":
			return self.__dict__['_relaxConstraints'].set(value)
		if attrName == "relaxDist":
			self.___dict__['relaxDist'].setvalue(str(value))
		if attrName == "relaxAngle":
			self.___dict__['relaxAngle'].setvalue(str(value))
		if attrName == "useRelaxColor":
			return self.___dict__['useRelaxColor'].set(value)
		if attrName == "relaxColor":
			return self.___dict__['relaxColor'].set(value)
		self.__dict__[attrName] = value
