from anytree import NodeMixin

class Atom(NodeMixin):
    """Creates nodes for the tree."""

    def __init__(self, name, direction=None, parent=None, children=None):
        """Every node has the atom type, and direction relative to the previous atom."""

        super(Atom, self).__init__()
        self.name = name
        self.direction = direction
        self.parent = parent
        if children:
            self.children = children