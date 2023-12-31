# --- UCSF Chimera Copyright ---
# Copyright (c) 2000-2012 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: __init__.py 39983 2014-07-14 20:58:21Z gregc $

from HTMLParser import HTMLParser

from send_to_chimera import getWebFileLock, getWebFile

import chimera
from chimera.baseDialog import ModelessDialog
from chimera.baseDialog import ModalDialog
from chimera import replyobj
import os, sys

from tcl_sockets import Server

class PUPPET_ERROR(Exception):
    pass   # "Unable to initialize Puppet to listen for requests!"

class AC_KEY_ERROR(Exception):
    pass   # "Error while accessing web info file"

class NO_AUTH(Exception):
    pass   # "Bad Credentials"

class NO_URL(Exception):
    pass   # "Error while retrieving file from web"

class CANCEL_FETCH(Exception):
    pass ## when the user decides to not supply authentication info
         ## i.e. presses 'Cancel'

puppet = None
gui_opener = None

_consider_yourself_warned = False

class PasswdAsker(ModalDialog):
    title = "Authentication Required"
    buttons = ('OK', 'Cancel')
    highlight = 'OK'

    def __init__(self, host, realm):
        self.realm = realm
        self.host = host
        ModalDialog.__init__(self, oneshot=1)

    def fillInUI(self, parent):
        import Tkinter

        self.authFrame = Tkinter.Frame(parent)
        self.authFrame.pack(side='top', fill='both')

        self.authInfo = Tkinter.Text(self.authFrame, relief='flat',
                                     wrap='word', pady=5, height=3,
                                     width=40)
        self.authInfo.insert(0.0,
                "Authentication required to access '%s' on host '%s'. " \
                "Please enter a valid username and password. " % (self.realm, self.host))
        self.authInfo.configure(state="disabled", background = self.authFrame.cget('background'))
        self.authInfo.pack(side='top', fill='x')


        import Pmw
        self.usernameEntry = Pmw.EntryField(self.authFrame, labelpos='w',
                                            label_text = "username: ")
        self.usernameEntry.pack(side='top', fill='x', pady=2)
        self.usernameEntry.component('entry').focus()

        self.passwdEntry = Pmw.EntryField(self.authFrame, labelpos='w',
                                          label_text = "password: ",
                                          entry_show = "*",
                                          command = self.OK)
        self.passwdEntry.pack(side='top', fill='x', pady=2)

        self.rememberVar = Tkinter.IntVar()
        self.rememberButton = Tkinter.Checkbutton(self.authFrame,
                    text="Remember my username and password",
                    variable = self.rememberVar, command=self.rememberPassword)
        self.rememberButton.pack(side='top', anchor='w', pady=3)

    def rememberPassword(self):
        pass

    def OK(self):
        self.usernameVal = self.usernameEntry.get()
        self.passwdVal = self.passwdEntry.get()

        ModalDialog.Cancel(self, value=1)

    def getAuthInfo(self):
        return (self.usernameVal, self.passwdVal)

    def getRememberVar(self):
        return self.rememberVar.get()


class WarnUserDialog(ModalDialog):

    title = "Warning!!"
    buttons = ('Yes', 'No', "Details")
    provideStatus = False
    help = "ContributedSoftware/webdata/webdata.html"


    def __init__(self, details_text):
        self.details_text = details_text

        ModalDialog.__init__(self, oneshot=1)
        self.details_shown = False

    def fillInUI(self, parent):
        import Tkinter, Pmw

        self.parent = parent

        warn_t = Tkinter.Text(parent, padx=10, pady=5,
                width=40, height=6, relief='flat', wrap='word')

        warning_text = "Some files can contain code which could harm " \
                       "your computer.  If you do not trust the source " \
                       "of this file, you should not open it in Chimera. " \
                       "\n\n" \
                       "Would you like to open this file in Chimera?"

        warn_t.insert('0.0', warning_text)

        warn_t.configure(state='disabled', background=self.uiMaster().cget('background'))
        warn_t.grid(row=0, column=0, sticky="nsew")
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)

        from waprefs import WA_WARN_ONCE, WA_WARN_ALWAYS, WA_WARN_NEVER

        self.warn_levels = [WA_WARN_ONCE, WA_WARN_ALWAYS, WA_WARN_NEVER]

        from waprefs import WEBACCESS_PREF, WA_WARN_LEVEL
        cur_ask_lvl = self.warn_levels[ (chimera.preferences.get(WEBACCESS_PREF, WA_WARN_LEVEL))-1 ]


        askLevel = Pmw.OptionMenu(parent,
                labelpos='w', label_text="Warn me about these files:",
                command=self.processAskSel, items=self.warn_levels,
                initialitem=cur_ask_lvl)

        askLevel.grid(row=1, column=0, sticky='nw', padx=10, pady=5)
        parent.rowconfigure(1, weight=1)

        ##################################
        #self.ask_again_var = Tkinter.IntVar(parent)
        #self.ask_again_var.set(0)
        #
        #ask_again_cb = Tkinter.Checkbutton(parent, variable=self.ask_again_var,
        #                                   padx=10, pady=10,
        #                                   text="Do not ask again this session"
        #                                   )
        #ask_again_cb.grid(row=1, column=0, sticky="nw")
        #parent.rowconfigure(1, weight=1)
        ##################################

        from chimera.HtmlText import HtmlText
        self.code_t = Pmw.ScrolledText(parent,
                text_pyclass = HtmlText, text_relief = 'sunken',
                text_wrap = 'word', text_width=40, text_height=8)

        #self.code_t = Pmw.ScrolledText(parent,
        #                               text_relief='sunken',
        #                               text_wrap='word',
        #                               text_width=40,
        #                               text_height=8,
        #                               )

        self.code_t.configure(text_state='disabled', hscrollmode='dynamic',
                vscrollmode='dynamic',
                text_background=parent.cget('background'))

        self.code_t.settext(self.details_text)

        #self.code_t.grid(row=1, column=0, sticky="s")
        #parent.rowconfigure(1, weight=1)

        #self.code_t.forget()

    def processAskSel(self, val):
        from waprefs import WEBACCESS_PREF, WA_WARN_LEVEL
        chimera.preferences.set(WEBACCESS_PREF, WA_WARN_LEVEL, self.warn_levels.index(val) + 1, asSaved=1)
        chimera.preferences.save()


    def Yes(self):
        ModalDialog.Cancel(self, value='yes')

    def No(self):
        ModalDialog.Cancel(self, value='no')

    def Details(self):
        if not self.details_shown:
            self.code_t.grid(row=2, column=0, sticky='s')
            self.parent.rowconfigure(2, weight=1)
            self.buttonWidgets['Details'].configure(text="Hide Details")
            self.details_shown = True
        else:
            self.code_t.grid_forget()
            self.buttonWidgets['Details'].configure(text="Details")
            self.details_shown = False


class ListenStatusDialog(ModelessDialog):

    title = "DBPuppet Status"
    buttons = ('OK')
    highlight = 'OK'

    def __init__(self, text):
        self.text = text
        ModelessDialog.__init__(self, oneshot=1)

    def fillInUI(self, parent):
        import Tkinter
        l = Tkinter.Label(parent, text = "%s" % self.text)
        l.grid(row=0, column=0, sticky="nsew")
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)


import urllib
class GUI_URLopener(urllib.FancyURLopener):

    def __init__(self, *args, **kwargs):
        self.clear_cache = False

        urllib.FancyURLopener.__init__(self, *args, **kwargs)

    def prompt_user_passwd(self, host, realm):
        pa = PasswdAsker(host, realm)
        res = pa.run(chimera.tkgui.app)

        if res:
            username, passwd = pa.getAuthInfo()
            if not pa.getRememberVar():
                self.clear_cache = True
        else:
            raise CANCEL_FETCH("canceled fetch")
            #username, passwd = None, None

        return username, passwd

    def get_user_passwd(self, host, realm, clear_cache = 0):
        u, p = urllib.FancyURLopener.get_user_passwd(self, host, realm, clear_cache)
        if self.clear_cache:
            self.reset_cache()
            self.clear_cache = False

        return u, p

    def reset_cache(self):
        self.auth_cache = {}

    def http_error_401(self, url, fp, errcode, errmsg, headers, data=None):
        """Error 401 -- authentication required.
        See this URL for a description of the basic authentication scheme:
        http://www.ics.uci.edu/pub/ietf/http/draft-ietf-http-v10-spec-00.txt"""

        #print "401 CALLED"
        #print "HEADERS IS ", headers

        if not 'www-authenticate' in headers:
            urllib.URLopener.http_error_default(self, url, fp,
                                                errcode, errmsg, headers)
        stuff = headers['www-authenticate']
        #print "STUFF IS ", `stuff`
        import re

        for s in stuff.split(", "):
            match = re.match('[ \t]*([^ \t]+)[ \t]+realm="([^"]*)"', s)
            if not match:
                continue

            scheme, realm = match.groups()

            if scheme.lower() != 'basic':
                continue

            name = 'retry_' + self.type + '_basic_auth'
            if data is None:
                return getattr(self, name)(url, realm)
            else:
                return getattr(self, name)(url, realm, data)

        urllib.URLopener.http_error_default(self, url, fp,
                                            errcode, errmsg, headers)

    def http_error_404(self, url, fp, errcode, errmsg, headers, data=None):
        raise IOError, "Received HTTP error 404 ('%s')" % errmsg


class DBPuppet:

    def __init__(self, name):
        self.server = None

        ## for security purposes. path to a file that will hold a 'key' which
        ## will allow Helper application to communicate with Chimera
        self.listen_port = None
        self.name = name

        # make sure that self.keys is initialized
        import random
        self.keys = []
        for r in range(3):
            self.keys.append(random.randint(0, 1000000))

    def printMsg(self, msg):
        sys.__stdout__.write("[%s]:  %s" % (self.name, msg))
        sys.__stdout__.flush()

    def portOnTop(self, triggerName=None, closure=None, data=None):
        if not self.listen_port:
            return

        try:
            lock = getWebFileLock()
            if lock is None:
                return
            with lock:
                self.remove_webinfo(self.listen_port, lock=lock)
                self.write_webinfo(self.listen_port, lock=lock)
        except (AC_KEY_ERROR, OSError) as what:
            replyobj.error("Web Data:  %s\n" % what)
            self.stop_listening()

    def stop_listening(self, triggerName=None, closure=None, data=None):
        """make the DBPuppet extension stop listening to requests from browsers"""

        try:
            self.remove_webinfo(self.listen_port, removeFile=True)
        finally:
            self.server.close()
            self.server = None
            self.listen_port = None

            chimera.triggers.deleteHandler("Chimera exit", self.exitHandler)
            chimera.triggers.deleteHandler("Focus In", self.focusHandler)

    def write_webinfo(self, port, lock=None):
        if lock is None:
            lock = getWebFileLock()
            if lock is None:
                return
        with lock:
            keyfile = getWebFile('a')
            if not keyfile:
                keyfile = self._new_file()
            with keyfile:
                keyfile.write("%s,%s\n" % (port, os.getpid()))

                try:
                    keyfile.close()
                except IOError, e:
                    replyobj.warning('Unable to update web info file: %s' % e)

    def remove_webinfo(self, port, removeFile=False, lock=None):
        ## clean up after yourself!!!
        if lock is None:
            lock = getWebFileLock()
            if lock is None:
                return
        with lock:
            keyfile = getWebFile('r+')
            if not keyfile:
                if not chimera.debug:
                    return
                raise AC_KEY_ERROR("can't open web info file for reading and writing")

            with keyfile:
                new_contents = []
                found_me = False
                my_pid = os.getpid()
                port = int(port)

                for l in keyfile.readlines():
                    try:
                        the_port, the_pid = l.split(",")
                        the_port = int(the_port)
                        the_pid = int(the_pid)
                    except ValueError:
                        ## must be the first line of keys
                        if len(l.split()) == 3:
                            new_contents.append(l)
                        continue

                    ## is a line of format port,pid
                    if the_pid == my_pid:
                        ## this entry is the one from *this* Chimera
                        ## or an older chimera with the same pid
                        found_me = True
                    else:
                        ## i.e. this entry is not the one from *this* Chimera
                        if self.check_pid_exists(the_pid):
                            new_contents.append(l)

                if not found_me:
                    return

                if len(new_contents) <= 1 and removeFile:
                    keyfile.close()
                    try:
                        os.remove(keyfile.name)
                    except OSError:
                        if not chimera.debug:
                            return
                        replyobj.error("error removing webinfo file %s" % keyfile.name)
                    return

                keyfile.seek(0)
                keyfile.writelines(new_contents)
                keyfile.truncate()
                try:
                    keyfile.close()
                except IOError, e:
                    replyobj.warning('Unable to update web info file: %s' % e)

    import sys
    if sys.platform.startswith('win'):
        def check_pid_exists(self, pid):
            ## On Windows, assume it exists
            return True
    else:
        def check_pid_exists(self, pid):
            try:
                os.kill(pid, 0)
            except OSError, what:
                import errno
                if what.errno == errno.ESRCH:
                    ## the process no longer exists...
                    return False
                elif what.errno == errno.EPERM:
                    ## someone else's process - shouldn't get here!!
                    return False
            # play it safe and assume it exists
            return True

    def _new_file(self):
        ## return new keyfile
        ## limit file to owner readable/writable. This will ensure that
        ## chimera and client are being run by the same user.
        old_mask = os.umask(0o177)
        keyfile = getWebFile('w', create=True)
        os.umask(old_mask)
        if not keyfile:
            raise AC_KEY_ERROR("can't create web info file")

        ## Write some random numbers into the key file.  The client will
        ## have to supply these numbers if it wants to do business.
        keystring = ' '.join(str(k) for k in self.keys)
        keyfile.write(keystring + "\n")
        return keyfile

    def start_listening(self):
        """make the DBPuppet extension listen to requests from browsers"""

        if self.server:
            ## already listening
            return

        if chimera.nogui:
            ## Since DBPuppet uses TCL ports to listen, it's not
            ## available in nogui mode where there is no tcl event loop
            raise PUPPET_ERROR("DBPuppet extension not available in nogui mode")

        try:
            lock = getWebFileLock()
            if lock is None:
                from send_to_chimera import _web_path
                replyobj.warning("Unable to receive web data: %s"
                        % _web_path[len('error:'):])
                return
            with lock:
                pass
        except (OSError, IOError):
            replyobj.warning("Unable to lock web data")
            return
        with lock:
            keyfile = getWebFile('r')
            if keyfile:
                with keyfile:
                    try:
                        keys = [int(k) for k in keyfile.readline().split()]
                        if len(keys) != 3:
                            raise ValueError
                        self.keys = keys
                    except ValueError:
                        ## seems as if the key file is corrupted, remove it...
                        keyfile.close()
                        try:
                            os.remove(keyfile.name)
                        except OSError:
                            replyobj.error("Web data key file is corrupted, "
                                           "and can't be removed. Please remove"
                                           "file %s" % keyfile.name)
                            return

        ## new server socket
        from Tkinter import TclError
        try:
            self.server = Server(chimera.tkgui.app, socket_connection_cb, self.keys)
            self.listen_port = self.server.getPortNo()
        except TclError:
            ## port is probably in use, probably by another Chimera on the machine.
            self.server = None
            raise PUPPET_ERROR("port %s appears to be in use" % self.listen_port)

        self.write_webinfo(self.listen_port)
        self.exitHandler = chimera.triggers.addHandler("Chimera exit", self.stop_listening, None)
        self.focusHandler= chimera.triggers.addHandler("Focus In", self.portOnTop, None)


class StripHTMLParser(HTMLParser):
    """This is a very general class that is written to get all the information
    between <pre> </pre> tags in an HTML document. It is needed because the URLs
    that were being fetched from ModBase were actual HTML documents and not just
    text files..."""

    def __init__(self):
        self.in_pre = 0
        HTMLParser.__init__(self)
        self.content=''

    def handle_starttag(self, tag, attrs):
        if tag == 'pre':
            self.in_pre = 1
    def handle_endtag(self, tag):
        if tag == 'pre':
            self.in_pre = 0
    def handle_data(self, data):
        if self.in_pre:
            self.content += data
    def get_content(self):
        return self.content


def stripHTML(file):
    """pass in a file in HTML, and it will extract all the information
    between the <pre> and </pre> tags. Possibly only useful for ModBase
    """

    shp = StripHTMLParser()
    f = open(file, 'r')
    shp.feed(f.read())
    f.close()

    f = open(file, 'w')
    f.write((shp.get_content()).strip())
    f.close()

def dangerous(loc):
    """return None if not dangerous, otherwise return
    the prefix or suffix that mached"""

    def findSuffix(loc, suffixes):
        """return matched suffix if any"""
        base, ext = os.path.splitext(loc)
        ext = ext.lower()
        return ext in suffixes and ext or None

    def findPrefix(loc, prefixes):
        """return matched prefix if any"""
        tmp = loc.split(':')
        if len(tmp) <= 1:
            return None
        return tmp[0] in prefixes and tmp[0] or None

    # loop through chimera.fileInfo for dangerous ones
    # and check those prefixes and suffices
    # TODO: cache this somehow
    prefixes = []
    suffixes = []
    infoTypes = chimera.fileInfo.types()
    for t in infoTypes:
        if not chimera.fileInfo.dangerous(t):
            continue
        prefixes.extend(chimera.fileInfo.prefixes(t))
        suffixes.extend(chimera.fileInfo.extensions(t))
    return findPrefix(loc, prefixes) or findSuffix(loc, suffixes)

def getURL(url, filename, useTempDir=None):
    """pass in a the url and the name of the file you'd like to save it to, and
    this function will return the location of that file"""

    tmp_path = getTmpFilename(filename, useTempDir)

    global gui_opener
    if not gui_opener:
        gui_opener = GUI_URLopener()

    try:
        f, h = gui_opener.retrieve(url=url, filename=tmp_path)
    except IOError, what:
        if 401 in what.args:
            gui_opener.reset_cache()
            raise NO_AUTH("Not authorized")
        else:
            raise NO_URL(what)

    return tmp_path

def getTmpFilename(filename, useTempDir):
    import tempfile

    orig_base = filename.split(".", 1)[0]
    orig_exts = filename.split(".", 1)[1]

    tmp_name = filename

    if useTempDir:
        tmp_dir = useTempDir
    else:
        tmp_dir = tempfile.gettempdir()

    count = 1

    while os.path.exists(os.path.join(tmp_dir, tmp_name)):
        tmp_name = orig_base + "-%d" % count + ".%s" % orig_exts
        count = count + 1

    return os.path.join(tmp_dir, tmp_name)


def doMidasCommand(cmd):
    """execute a midas command in Chimera"""

    from Midas import midas_text, MidasError
    try:
        midas_text.makeCommand(cmd)
    except MidasError, e:
        chimera.replyobj.error("%s" % e)

def activate_puppet(state):
    if chimera.nogui:
        # not supported
        return
    run_puppet(state.get())

def run_puppet(on):
    global puppet

    if on:
        ##if you have not already created an instance
        if not puppet:
            ## create an instance
            puppet = DBPuppet("Database Puppet")

            ## start listening
            try:
                puppet.start_listening()
            except AC_KEY_ERROR, what:
                replyobj.error("Web Data:  %s\n" % what)
                puppet = None

            except PUPPET_ERROR, what:
                replyobj.error("%s:  %s" % (PUPPET_ERROR, what))
                puppet = None

            else:
                replyobj.status("Chimera is now accepting web data...")

        else:
            ## you already have an instance
            replyobj.status("Chimera is already accepting web data!!")

    else:
        if not puppet:
            replyobj.status("Chimera isn't accepting web data yet!!")
            return
        else:
            try:
                puppet.stop_listening()
            except AC_KEY_ERROR, what:
                replyobj.error("Web Data:  %s\n" % what)

            puppet = None
            replyobj.status("Chimera is no longer accepting web data...")


def ConfigBrowser():
    import BrowserCfgGUI
    chimera.dialogs.register("Browser Configuration", BrowserCfgGUI.BrowserCfgGUI, replace=1)
    chimera.dialogs.display("Browser Configuration")


def verify_connection(socket, keys):
    if socket.from_host != '127.0.0.1':
        replyobj.error("Attempted connection from host %s. Closing connection." % socket.from_host)
        return 0

    socket.write("CHIMERA\n")

    key = socket.readline()

    k = key.split()
    if len(k) != 3:
        replyobj.error("ERROR: misformatted key, closing connection....")
        socket.write("NO\n")
        return 0

    ## make sure the keys supplied by the client (read over the socket stream)
    ## match the keys you wrote

    if [int(i) for i in k] != keys:
        replyobj.error("ERROR: keys don't match, closing connection....")
        socket.write("NO\n")
        return 0
    else:
        socket.write("OK\n")
        return 1

def determineFileType(f):
    module = None

    for line in f.readlines():
        if line.find("ChimeraPuppet") >= 0:
            module = line[line.find("\"")+1:line.rfind("\"")]
            break

    if not module:
        raise chimera.UserError, "Couldn't find xml \"type\" attr. in %s" % f.name

    return module

def socket_connection_cb(socket, keys):
    if not verify_connection(socket, keys):
        socket.close()
        return

    ## just reading path, not contents
    data = socket.readline()
    data = data.strip()

    ## open what was sent over socket (i.e path to dloaded file),
    ## not what you just wrote

    if needToWarn():
        if dangerous(data):
            warning_dlg = WarnUserDialog(
                    "No details available for this file.\n\n"
                    "If you do not trust the source of "
                    "this file, save it to your computer "
                    "as a text file and review the contents "
                    "before opening in Chimera.")
            res = warning_dlg.run(chimera.tkgui.app)
            if res == 'no':
                socket.write("OPEN OK\n")
                socket.close()
                return

    try:
        chimera.openModels.open(data, prefixableType=True)
    except:
        socket.write("OPEN ERROR\n")
        socket.close()
        raise
    else:
        socket.write("OPEN OK\n")
        socket.close()

def _openWebData(file_loc):
    f = open(file_loc, 'r')
    module = determineFileType(f)
    f.close()

    try:
      #print "trying to import %s" % module
      mod = __import__(module, globals(), locals(), [], 1)
    except ImportError, what:
      chimera.replyobj.error("Module Missing: %s" % what)
      return []

    ## need to do this to account for the paths returned by
    ## mkstemp on windows - "\\" is displayed as '\' when printed
    ## This is passed into module.handle_file, and python tries to escape
    ## whatever is after the '\', as is C:\documents and settings\administrator
    new_loc = file_loc.replace("\\", "/")
    cls = getattr(mod, module)
    cls().handle_file(new_loc)
    return []


def needToWarn():
    from waprefs import WEBACCESS_PREF, WA_WARN_LEVEL
    current_warning_pref = chimera.preferences.get(WEBACCESS_PREF, WA_WARN_LEVEL)

    if current_warning_pref == 1:
        if not getWarned():
            setWarned(True)
            return True
        else:
            return False
    elif current_warning_pref == 2:
        return True
    elif current_warning_pref == 3:
        return False


def getWarned():
    return _consider_yourself_warned


def setWarned(val):
    global _consider_yourself_warned
    _consider_yourself_warned = val


from chimera import fileInfo
fileInfo.register('Chimera web data', _openWebData,
        ['.chimerax'], ['chimerax'], mime=['application/x-chimerax'],
        canDecompress=False, dangerous=True, category=fileInfo.SCRIPT)


def _onTopEvent(app, event):
    chimera.triggers.activateTrigger('Focus In', None)


if not chimera.nogui:
    chimera.triggers.addTrigger('Focus In')
    chimera.tkgui.app.bind('<Enter>', lambda event, w=chimera.tkgui.app: _onTopEvent(w, event))
