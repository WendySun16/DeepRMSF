# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: prefs.py 30896 2010-07-07 00:32:20Z pett $

from chimera import preferences

HELIX_COLOR = "helix color"
SHEET_COLOR = "sheet color"
COIL_COLOR = "other color"

defaultColors = {
	HELIX_COLOR: 'orange red',
	SHEET_COLOR: 'purple',
	COIL_COLOR: 'gray'
}
options = {}
options.update(defaultColors)
prefs = preferences.addCategory("ColorSS", preferences.HiddenCategory,
							optDict=options)
