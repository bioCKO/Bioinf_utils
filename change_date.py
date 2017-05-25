__author__ = 'mjohnpayne'

import os.path, time, sys

file = sys.argv[1]

print "last modified: %s" % time.ctime(os.path.getmtime(file))
print "created: %s" % time.ctime(os.path.getctime(file))

