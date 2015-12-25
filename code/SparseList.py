class SparseListOfList():
    def __init__(self):
        self.d = dict()

    def __setitem__(self, index, value):
        self.d[index] = value

    def __getitem__(self, index):
        try:
            return self.d[index]
        except KeyError:
            x = []
            self.d[index] = x
            return x