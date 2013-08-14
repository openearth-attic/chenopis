# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import argparse
import sys
import re

# <codecell>

fortran_files = [os.path.join('../src', filename) for filename in os.listdir('../src') if filename.endswith('.F') or filename.endswith('.F90')]

# <codecell>

preprocess_lines = []
preprocess_re = re.compile(r'^!(?P<tag>[A-Z]{2,5})(?P<code>.*)')

def upgrade_preprocess(filename):
    lines = open(filename).readlines()
    newlines = []
    tag = None
    for i, line in enumerate(lines):
        match = preprocess_re.search(line)
        if match:
            newtag = tag is None
            differenttag = not newtag and tag != match.group('tag')
            if differenttag:
                newlines.append("#endif\n")
            if newtag or differenttag:
                newlines.append("#ifdef HAVE_{}\n".format(match.group('tag')))
            tag = match.group('tag')
            # replace tag by spaces of the same length
            newlines.append(' '*(len(tag)+1) + match.group('code') + '\n')
        else:
            endtag = tag is not None
            if endtag:
                newlines.append("#endif\n")
            tag = None
            newlines.append(line)
    newfile = open(filename, 'w')
    newfile.writelines(newlines)
    newfile.close()
    
                

# <codecell>

for filename in fortran_files:
    upgrade_preprocess(filename)

# <codecell>


