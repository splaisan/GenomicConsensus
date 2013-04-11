

class WorkChunk(object):
    """
    A chunk of the reference
    """
    def __init__(self, window, hasCoverage):
        self.window      = window
        self.hasCoverage = hasCoverage
