from __future__ import print_function

import os
import re
from glob import glob


def parser_reference(file_handle, hyper=False):
    """
    parser_reference(file_handle: file, hyper: bool) -> dict

    Parse reference information for each ref in ref.html page

    Arguments
    ---------
    file_handle: file
        file handle object for parsing reference information
    hyper: bool
        whether keep the hyper link at the end of each reference
        False for not reserve, otherwise True, default is False

    Reason:
        <a></a> are not allowed to nested in html

    Returns
    -------
    {int : str}
        int is the index of reference on the ref page
        str is the ref information
    """
    line = file_handle.readline()
    contents = {}
    while line:
        if line.startswith("<tr>"):
            index_part, cont_part = line.split(']</td><td>')
            index = int(index_part.split('[')[1])
            cont = cont_part.split('</td>')[0]
            if hyper is False:
                cont = cont.split('<a class')[0].strip()
            contents[index] = cont
        line = file_handle.readline()
    return contents


def glob_sci_html(path):
    """
    glob_sci_html(path: str) -> None

    Seach over given path and tweak all the qualified files

    Argument
    --------
    path: the path to grob the desired html

    This is the function change the html contents

    """
    # file_path = path + 'sci_doc_*.html'
    file_name = glob(path + 'sci_doc_*.html')
    bib_path = path + 'sci_doc_zref.html'
    with open(bib_path) as f:
        contents = parser_reference(f, hyper=False)
        # print(contents)
    for i in file_name:
        new_html = ''
        flag = True # True for new html, False for tweaked one
                    # The flag here is prevent duplicate changes
        with open(i, 'r') as f:
            line = f.readline()
            while line:
                position = line.find("id=\"id") # find the hyper in html
                if position > -1:
                    # line = re.sub(pattern, 'tooltip', line)
                    ref_position = line.find('reference internal')
                    if ref_position == -1:
                        flag = False
                        break
                    line = line.replace('reference internal', 'tooltip')
                    line = add_hover_box(line, contents)
                    new_html += line
                else:
                    new_html += line
                line = f.readline()
        # print(new_html)
        if flag:
            with open(i, 'w') as f:
                f.write(new_html)
        # break


def add_hover_box(line, bib_contents):
    """
    add_hover_box(line: str, bib_contents: str) -> str

    Add the bib_contests to the hover window of a specific hyperlink

    Arguments
    ---------
    line: the line contains reference hyperlink
    bib_contents: the ref info will show up in the hover window

    Returns
    -------
    The new string with a '<span class="tooltiptext">bib info</span>'
    """
    extra_text = "<span class=\"tooltiptext\">{}</span>"
    r_index = 0
    newline = ''
    while True:
        l_index = line.find('[')
        if l_index == -1:
            newline += line[:]
            break
        else:
            r_index = line.find(']')
            number = int(line[l_index+1: r_index])
            contents = bib_contents[number]
            newline += line[:r_index + 1] + extra_text.format(contents)
            line = line[r_index + 1:]
    return newline

def _test_content():
    """test parser function
    """
    path = relative_path + '/_build/html/'
    file_path = path + 'sci_doc_zref.html'
    with open(file_path) as f:
        print('open successly')
        contents = parser_reference(f, hyper=False)
        print(contents)



if __name__ == "__main__":
    relative_path = os.path.dirname(os.path.realpath(__file__))
    path = relative_path + '/_build/html/'
    file_path = path + 'sci_doc_zref.html'
    glob_sci_html(path)
    # _test_content()
