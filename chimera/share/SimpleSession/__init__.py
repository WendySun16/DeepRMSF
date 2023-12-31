# --- UCSF Chimera Copyright ---
# Copyright (c) 2000 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: __init__.py 41991 2019-01-16 21:23:46Z pett $

"""Save/restore most aspects of a Chimera modeling session"""

from save import saveSession, sessionID, noAutoRestore, autoRestorable, \
	sesRepr, SAVE_SESSION, BEGIN_RESTORE_SESSION, END_RESTORE_SESSION, \
	registerAttribute, colorID, isRegisteredAttribute
RESTORE_SESSION = END_RESTORE_SESSION
