# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: __init__.py 39985 2014-07-15 18:35:35Z pett $

import os.path
import Tkinter
import chimera
from chimera.tkoptions import InputFileOption
from types import StringTypes

formatName = "MMTK"

class ParamGUI:
	labels = ["NetCDF"]

	def __init__(self, parent):
		from Trajectory.prefs import prefs, INPUT_FILES
		inputPrefs = prefs[INPUT_FILES].setdefault('MMTK', {})
		self.options = []
		for i in range(len(self.labels)):
			label = self.labels[i]
			self.options.append(InputFileOption(parent, i,
				label, inputPrefs.get(label, True),
				None, title="Choose %s File" % label,
				historyID="MMTK %s" % label))
		parent.columnconfigure(1, weight=1)

	def loadEnsemble(self, startFrame, endFrame, callback):
		args = []
		from Trajectory.prefs import prefs, INPUT_FILES
		# need to change a _copy_ of the dictionary, otherwise
		# when we try to save the "original" dictionary will also
		# have our changes and no save will occur
		from copy import deepcopy
		inputPrefs = deepcopy(prefs[INPUT_FILES])
		for i in range(len(self.labels)):
			path = self.options[i].get()
			label = self.labels[i]
			if not os.path.exists(path):
				raise ValueError, \
					"%s file '%s' does not exist!" % (
								label, path)
			inputPrefs['MMTK'][label] = path
			args.append(path)
		prefs[INPUT_FILES] = inputPrefs

		loadEnsemble(args, startFrame, endFrame, callback)

def loadEnsemble(inputs, startFrame, endFrame, callback, relativeTo=None):
	from MMTKNetcdf import MMTKNetcdf
	if relativeTo:
		import os.path
		for i, f in enumerate(inputs):
			if type(f) not in StringTypes or os.path.isabs(f):
				continue
			inputs[i] = os.path.join(relativeTo, f)
	ensemble = MMTKNetcdf(*tuple(inputs + [startFrame, endFrame]))
	from chimera import replyobj
	replyobj.status("Creating interface", blankAfter=0)
	try:
		callback(ensemble, keepLongBonds=True)
	finally:
		replyobj.status("Interface created")
