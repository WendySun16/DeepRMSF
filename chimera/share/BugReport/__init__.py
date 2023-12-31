# vi:sw=4:
import chimera, sys
import httplib, mimetypes, urlparse

ADDL_INFO =  "(Describe the actions that caused this problem to occur here)"

def canonPlatName(plat):
    plat_names = {
        'linux2'    : "Linux",
        'osf1V5'    : "Alpha Tru64 Unix",
        'win32'     : "Windows",
        'irix6-n32' : "SGI IRIX",
        'darwin'    : "Mac OS X"
    }

    try:
        return plat_names[plat]
    except KeyError:
        return plat

def encode_multipart_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded
    as files
    Return (content_type, body) ready for httplib.HTTP instance
    """
    BOUNDARY = '----------ThIs_Is_tHe_bouNdaRY_$'
    CRLF = '\r\n'
    L = []
    for (key, value) in fields:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        if isinstance(value, unicode):
            value = value.encode('utf8')
        L.append(value)
    for (key, filename, value) in files:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
        L.append('Content-Type: %s' % get_content_type(filename))
        L.append('Content-Length: %d' % len(value))
        L.append('')
        L.append(value)
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body

def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'


class BugReport:

    def __init__(self, name='', email='', description='', info='', filename='',
            version='', platform='', includeSysInfo=True):
        self.name        = name
        self.email       = email
        self.description = description
        self.info        = info
        self.filename    = filename

        self.version     = version or chimera.version.release
        import platform as pyplat
        self.platform    = platform or ("%s%s (%s) %s" % (
                        pyplat.system(), "64" if sys.maxsize > 2 ** 32 else "",
                        chimera.opengl_platform(), chimera.operating_system()))

        if self.description:
            self.description = self.description.replace("&", "&amp;")
            self.description = self.description.replace("<", "&lt;")
            self.description = self.description.replace(">", "&gt;")
            self.description = '<br>\n'.join(self.description.split('\n'))

        if self.info:
            self.info = self.info.replace("&", "&amp;")
            self.info = self.info.replace("<", "&lt;")
            self.info = self.info.replace(">", "&gt;")
            self.info = '<br>\n'.join(self.info.split('\n'))

        if not self.description:
            self.description = \
                        "<font color=\"blue\"><i>%s</i></font>\n" % ADDL_INFO
        
        if includeSysInfo:
            self.addSysInfo()

        self.addGraphicsSettings()

    def addGraphicsSettings(self):
        if self.info:
            self.info += "<hr>\n"
        self.info += graphicsSettings()

    def addSysInfo(self):
        if self.info:
            self.info += "<hr>\n"
        self.info += systemInfo()

    def getValue(self, val):
        try:
            return getattr(self, val)
        except AttributeError:
            return None

    def printStuff(self):
        ## for debugging

        print "NAME ", self.name
        print "EMAIL", self.email
        print "DESC ", self.description
        print "INFO ", self.info
        print "VERS ", self.version
        print "PLAT ", self.platform

def graphicsSettings():
        import chimera, Lighting
        v = chimera.viewer
        c = v.camera
        om = chimera.openModels
        cofrm = {om.Fixed:'fixed',
                 om.CenterOfModels:'center of models',
                 om.Independent:'independent',
                 om.CenterOfView:'center of view',
                 om.FrontCenter:'front center'}.get(om.cofrMethod, 'unknown')
        lod = chimera.LODControl.get()
        try:
                graphicsInfo = (
                    "<i>Multisampling</i>:   %s<br>\n" % (chimera.multisample,) +
                    "<i>Shadows</i>: %s<br>\n" % (v.showShadows,) +
                    "<i>Shadow texture size</i>: %s<br>\n" % (v.shadowTextureSize,) +
                    "<i>Silhouettes</i>:  %s<br>\n" % (v.showSilhouette,) +
                    "<i>Depth cue</i>:  %s<br>\n" % (v.depthCue,) +
                    "<i>Subdivision quality</i>:  %.2f<br>\n" % (lod.quality,) +
                    "<i>Single-layer transparency</i>:  %s<br>\n" % (v.singleLayerTransparency,) +
                    "<i>Transparent background</i>:  %s<br>\n" % (chimera.bgopacity,) +
                    "<i>Shaders supported</i>:  %s<br>\n" % (v.haveShaderSupport(),) +
                    "<i>Using shader</i>:  %s<br>\n" % (v.haveShader(),) +
                    "<i>Window size</i>:  %d %d<br>\n" % tuple(v.windowSize) +
                    "<i>Camera mode</i>:  %s<br>\n" % (c.mode(),) +
                    "<i>Orthographic projection</i>:  %s<br>\n" % (c.ortho,) +
                    "<i>Center of rotation</i>:  %s<br>\n" % (cofrm,) +
                    "<i>Near/far clipping</i>:  %s<br>\n" % (v.clipping,) +
                    "<i>Key light</i>:  %s<br>\n" % (v.keyLight != None,) +
                    "<i>Fill light</i>:  %s<br>\n" % (v.fillLight != None,) +
                    "<i>Back light</i>:  %s<br>\n" % (v.backLight != None,) +
                    "<i>Ambient light</i>:  %.2f<br>\n" % (v.ambient,) +
                    "<i>Specular sharpness</i>:  %.2f<br>\n" % (Lighting.sharpness(),) +
                    "<i>Specular reflectivity</i>:  %.2f<br>\n" % (Lighting.reflectivity(),))

                graphicsInfo += '<hr>\n' + openGLCapabilities()
        except:
                graphicsInfo = "<i>Error gathering graphics info</i>"
        return graphicsInfo

def openGLCapabilities():
        import OpenGLDebug, libgfxinfo as gi
        fni = OpenGLDebug.Features.items()
        fni.sort()
        flist = []
        for n,i in fni:
            s = gi.hasInfo(i)
            if s & 0x2:
                state = 'native'
            elif s & 0x1:
                state = 'extension'
            else:
                state = 'not supported'
            if s & 0x4:
                state += ', disabled'
            flist.append('<i>%s</i>: %s<br>' % (n,state))
        cap = '\n'.join(flist)
        return cap

def systemInfo():
        from chimera import SubprocessMonitor as SM
        from xml.sax.saxutils import escape
        openglInfo = (
            "<i>OpenGL Vendor</i>:   %s<br>\n"
            "<i>OpenGL Renderer</i>: %s<br>\n"
            "<i>OpenGL Version</i>:  %s<br>\n" % (
                escape(chimera.opengl_vendor()),
                escape(chimera.opengl_renderer()),
                escape(chimera.opengl_version()))
        )
        import os
        rootdir = os.getenv('CHIMERA')
        if not rootdir:
            bindir = ''
        else:
            bindir = os.path.join(rootdir, 'bin')
        if sys.platform == 'win32':
            prog = [os.path.join(bindir, 'machinfo.exe')]
        elif sys.platform.startswith('linux'):
            prog = [os.path.join(bindir, 'linuxdist')]
        elif sys.platform.startswith('darwin'):
            prog = [
                '/usr/sbin/system_profiler',
                '-detailLevel',
                'mini',
                'SPHardwareDataType',
                'SPSoftwareDataType',
                'SPDisplaysDataType'
            ]
        else:
            prog = []

        if not prog:
            osInfo = []
        else:
            try:
                subproc = SM.Popen(prog, stdin=None, stdout=SM.PIPE,
                                                                stderr=None)
                osInfo = subproc.stdout.readlines()
            except OSError:
                osInfo = []
            for i, line in enumerate(osInfo):
                line = line.strip().decode('utf-8', 'replace')
                colon = line.find(':')
                if colon == -1 or colon + 1 == len(line):
                        osInfo[i] = escape(line)
                        continue
                osInfo[i] = '<i>%s</i>%s' % (escape(line[:colon]), escape(line[colon:]))

        info = "%s%s" % (openglInfo, '<br>'.join(osInfo))
        return info


from chimera import dialogs, replyobj
from BugReportGUI import BugReportGUI, BugNotification


def displayDialog(wait=False):
    if dialogs.find(BugReportGUI.name):
        replyobj.status("Bug report already in progress!",
                        color="red", blankAfter=15)
        return None
    else:
        br_gui = dialogs.display(BugReportGUI.name, wait)
        return br_gui


def bugNotify(explanation, description):
    BugNotification(explanation, description)


dialogs.register(BugReportGUI.name, BugReportGUI)
