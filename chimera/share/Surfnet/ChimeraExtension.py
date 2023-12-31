# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---

import chimera.extension

class InterfaceSurfnetEMO(chimera.extension.EMO):
	def name(self):
		return 'Surfnet'
	def description(self):
		return 'display surface cavities between two sets of atoms'
	def categories(self):
		return ['Surface/Binding Analysis']
	def activate(self):
		self.module('choose').InterfaceSurfnetCB()
		return None

chimera.extension.manager.registerExtension(InterfaceSurfnetEMO(__file__))
