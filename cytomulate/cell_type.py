
class CellType:
    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.parent = None
        self.children = []
        self.marker_pattern = []
        self.expression_level = []
        self.gating_markers = []
