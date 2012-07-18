import sys
from viewer import DendroViewer

filename = False
if len(sys.argv) == 2:
    filename = sys.argv[1]
v = DendroViewer(filename)
v.run()

