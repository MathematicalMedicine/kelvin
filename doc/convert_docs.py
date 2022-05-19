#!/usr/bin/env python
# Copyright (C) 2015, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

# Author: Jo Valentine-Cooper <jvc@mathmed.org>

"""
Assembles our Markdown documentation into HTML files.

"""

import markdown
from pathlib import Path
import re

# in the Orig Usage docs, there's some column number lines that get inserted
# before certain code blocks. these take care of them.
class AddColumnHeaders(markdown.extensions.Extension):
    def __init__(self, *args, **kwargs):
        try:
            self.starters = kwargs['starters']
            del kwargs['starters']
        except KeyError:
            self.starters=None
        
        markdown.extensions.Extension.__init__(self, *args, **kwargs)
    
    def extendMarkdown(self, md, md_globals):
        md.preprocessors.add('codecolumnheader',
                ColumnHeaderAdder(starters=self.starters), "_begin")

class ColumnHeaderAdder(markdown.preprocessors.Preprocessor):
    """Adds a column number header to codeblocks beginning with the lines
    specified by the 'starters' kwarg.
    
    """
    
    # FIXME: A better, VASTLY cleaner way to do this would be to use a
    # POSTprocessor, and toy with the output CodeHilite creates. But that's a
    # bit too much work for today. :D
    
    def __init__(self, *args, **kwargs):
        try:
            self.starters = kwargs['starters']
            del kwargs['starters']
        except KeyError:
            pass
        
        markdown.preprocessors.Preprocessor.__init__(self, *args, **kwargs)
    
    def run(self, lines):
        new_lines = []
        for line in lines:
            # first check: is one of our 'code' blocks starting at this line?
            for starter in self.starters:
                if line.startswith("    {}".format(starter)):
                    # yes it is, so we create a header line.
                    
                    # start with two spaces to accommodate line numbers
                    colline = ["<pre class='colnum'><span>  ",]
                    
                    # we need to figure out how much to pad the values for our
                    # header line (since it ultimately comes out as
                    # preformatted text). there's a regex trick we can do to
                    # split and keep whitespace, but it gets us a bit of an
                    # awkward-looking list - something like this:
                    # ['', '   ', 'val', '   ', 'val2', '  ', 'val3']
                    # what we really want is:
                    # ['val   ', 'val2  ', 'val3']
                    # (that first ['', '    '] is the four spaces that marked
                    # this as a code block; we don't want to keep them)
                    # so a little postprocessing is done by doing that split
                    # and then iterating over it.
                    cols = iter(re.split(r'(\s+)', line))
                    pair_cols = [col + next(cols, '') for col in cols]
                    pair_cols.pop(0)
                    
                    # now that that's done, we go ahead and assemble the header
                    # line and add it to our output
                    colnum = 1
                    for col in pair_cols:
                        colline.append("{colnum: <{width}}".format(
                                colnum=colnum,
                                width=max(len(col), len(str(colnum)))))
                        colnum = colnum + 1
                    colline.append("</span></pre>")
                    new_lines.append("".join(colline))
                    new_lines.append("")
                            # space needed between our html and the codeblock;
                            # rendering gets screwed up otherwise
            
            # we don't want to get rid of any lines, just add them, so we
            # always append any line we look at no matter what. :)
            new_lines.append(line)
        return new_lines



def convert_docs():
    """
    Converts all our Markdown-formatted documentation files into HTML.
    
    """
    
    # markdown extensions and config
    colheaders = AddColumnHeaders(starters=
            ['212 101 0   0   1   1   3   4   1   3   2   2   1   2   1   1',
            'CHR  MARKER       KOSAMBI'])
    mdconf = {
            'extensions': ['toc', 'extra', 'codehilite', colheaders],
            'extension_configs': { 'codehilite': { 'linenums': True } },
            'encoding': 'utf-8'
            }
    
    # where all our files are
    docs = Path(__file__).parent
    
    # pull in header and footer
    header = open(docs / "raw" / "docsheader.html", mode="r", encoding="utf-8")
    header_html = header.read()
    header.close()
    
    footer = open(docs / "raw" / "docsfooter.html", mode="r", encoding="utf-8")
    footer_html = footer.read()
    footer.close()
    
    # convert all docs
    for mdfile in docs.glob('*.md'):
        with open(mdfile.with_suffix('.html'), mode="w", encoding="utf-8",
                errors="xmlcharrefreplace") as content, \
                open(mdfile, mode="r", encoding="utf-8",
                errors="xmlcharrefreplace") as mdsource:
            content.write(header_html)
            content.write(markdown.markdown(mdsource.read(), **mdconf))
            content.write(footer_html)

if __name__ == "__main__":
    convert_docs()

