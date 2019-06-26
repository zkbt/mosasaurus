from ...Night import Night

class Visit(Night):
    '''
    For Hubble observations, it makes more sense to organize data into visits.
    However, a visit is basically the same structure as a night: a group
    of exposures that are connected in time and with the same basic setup.
    '''

    def createVisitLog(self, remake=False):
        self.createNightlyLog(remake=remake)
