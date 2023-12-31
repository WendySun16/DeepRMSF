# Copyright (c) 2000 by the Regents of the University of California.
# All rights reserved.  See http://www.cgl.ucsf.edu/chimera/ for
# license details.
#
# $Id: chainPicker.py 41943 2018-10-18 18:07:23Z pett $

import chimera
from chimera.baseDialog import ModelessDialog
import Tkinter
import Pmw
from operator import add
from base import getPhysicalChains

class ChainPicker(ModelessDialog):
	title = 'Select Chains'
	buttons = ('OK', 'Apply', 'Close',)
	default = 'OK'
	oneshot = True
	help="UsersGuide/modelpanel.html#selchains"

	def __init__(self, models):
		self.models = models
		ModelessDialog.__init__(self)
	
	def fillInUI(self, parent):
		Tkinter.Label(parent, text="Select the following chains:"
					).grid(row=0, column=0, columnspan=2)
		parent.columnconfigure(0, weight=1)
		parent.columnconfigure(1, weight=1)
		parent.rowconfigure(1, weight=1)
		
		# chains by ID
		idGroup = Pmw.Group(parent, tag_text="by ID")
		idGroup.grid(row=1, column=0, sticky='nsew')
		idGroup.interior().columnconfigure(0, weight=1)
		idGroup.interior().rowconfigure(0, weight=1)
		idGroupFrame = Pmw.ScrolledFrame(idGroup.interior())
		idGroupFrame.grid(sticky="nsew")
		idRow = 0

		# physical chains
		physGroup = Pmw.Group(parent, tag_text="by connectivity")
		physGroup.grid(row=1, column=1, sticky='nsew')
		physGroup.interior().columnconfigure(0, weight=1)
		physGroup.interior().rowconfigure(0, weight=1)
		physGroupFrame = Pmw.ScrolledFrame(physGroup.interior())
		physGroupFrame.grid(sticky="nsew")
		physRow = 0

		# find the damn things
		for model in self.models:
			if not hasattr(model, 'bonds'):
				continue

			physical = getPhysicalChains(model)

			idCodes = {}
			for res in model.residues:
				idCodes[res.id.chainId] = 1

			Tkinter.Label(idGroupFrame.interior(), text="Model %s" %
					model.oslIdent()[1:]).grid(row=idRow,
					column=0)
			idRow = idRow + 1
			Tkinter.Label(physGroupFrame.interior(), text="Model %s" %
					model.oslIdent()[1:]).grid(row=physRow,
					column=0)
			physRow = physRow + 1

			def _idSort(id1, id2):
				# put single-letter codes first
				if len(id1) == 1 and len(id2) == 1:
					return cmp(id1, id2)
				if len(id1) == 1:
					return -1
				if len(id2) == 1:
					return 1
				return cmp(id1, id2)

			codeList = idCodes.keys()
			codeList.sort(_idSort)
			self._codeVars = {}
			for code in codeList:
				var = Tkinter.IntVar(parent)
				var.set(0)
				self._codeVars[(model, code)] = var
				if code == ' ':
					codeText = '(blank)'
				else:
					try:
						seq = model.sequence(code)
					except (AssertionError, KeyError):
						codeText = code
					else:
						description = getattr(seq, 'descriptiveName', None)
						if description:
							codeText = "%s: %s" % (code, description)
						else:
							codeText = code
				Tkinter.Checkbutton(idGroupFrame.interior(), pady=0,
					variable=var, text=codeText).grid(
					row=idRow, column=0, sticky='w')
				idRow = idRow + 1

			self._physVars = []
			for resList in physical[1:]:
				# find smallest/largest res number
				lowID = resList[0].id
				highID = resList[0].id
				for res in resList:
					if res.id.position > highID.position:
						highID = res.id
					elif res.id.position < lowID.position:
						lowID = res.id
				
				var = Tkinter.IntVar(parent)
				var.set(0)
				self._physVars.append([model, var, resList])
				if lowID == highID:
					text = lowID
				else:
					text = "%s-%s" % (lowID, highID)
				Tkinter.Checkbutton(physGroupFrame.interior(),
					pady=0, variable=var, text=text).grid(
					row=physRow, column=0, sticky='w')
				physRow = physRow + 1
			if physical[0]:
				var = Tkinter.IntVar(parent)
				var.set(0)
				self._physVars.append([model, var, physical[0]])
				Tkinter.Checkbutton(physGroupFrame.interior(),
					text="misc.", pady=0,
					variable=var).grid(row=physRow,
					column=0, sticky='w')
				physRow = physRow + 1

			if model != self.models[-1]:
				# add a "newline"
				Tkinter.Label(idGroupFrame.interior(), text=" "
						).grid(row=idRow, column=0)
				idRow = idRow + 1
				Tkinter.Label(physGroupFrame.interior(), text=" "
						).grid(row=physRow, column=0)
				physRow = physRow + 1

	def Apply(self):
		from operator import add
		from chimera.selection import ItemizedSelection, \
							OSLSelection, REPLACE
		sel = ItemizedSelection()
		chainOSL = ""
		for codeInfo, var in self._codeVars.items():
			model, code = codeInfo
			if var.get():
				chainOSL = chainOSL + model.oslIdent() \
								+ ":." + code
		if chainOSL:
			sel.merge(REPLACE, OSLSelection(chainOSL))
		
		for model, var, resList in self._physVars:
			if var.get():
				sel.add(resList)
				
		sel.addImplied(vertices=0)
		chimera.tkgui.selectionOperation(sel)
