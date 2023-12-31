import chimera
import chimera.extension

class MovieRecorder_EMO(chimera.extension.EMO):
    def name(self):
        return 'Movie Recorder'
    def description(self):
        return 'Make a movie from the contents of the graphics window'
    def categories(self):
        return ['Utilities']
    def icon(self):
        return self.path('film.png')
    def activate(self):
        self.module().showMRDialog()        

chimera.extension.manager.registerExtension(MovieRecorder_EMO(__file__))


def doMovieCmd(cmdName, args):
    from MovieRecorder import moviecmd
    moviecmd.movie_command(cmdName, args)

from Midas.midas_text import addCommand
addCommand("movie", doMovieCmd, help=True, changesDisplay=False)
