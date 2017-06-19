from ..imports import *

class Spectrograph(Talker):
    def __init__(self):

        Talker.__init__(self)
        self.setupDetector()
        self.setupDisperser()
        self.setupDirectories()
        self.setupExtraction()
