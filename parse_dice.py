#!/usr/bin/env python
'''
This code is meant to extract the comments from the Dice code and write them to
.rst files that can be prossessed by Sphinx. (For python 2.X)

Must be executed from the main Dice directory and will create a temporary
directory for the .rst files at the same level as the main Dice directory so
they will not be changed when the user switches branches to the gh-pages branch.

Make sure that your comments follow the proper indentation of the tree that
appears in the documentation of a function or class with look messy and jumbled.

Skeleton Format Example for Function
------------------------------------
bool SomeFunc(int Input1, std::vector<int>& Input2) {
    /*!
    Function description


    .. note:: Some note

    :Arguments:

        int Input1:
            description
        std::vector<int>& Input2:
            description

    //(if you return an object)
    :Returns:

        bool Output1:
            description

    */

Documentation of classes should follow the same pattern of indentation and
should include descriptions of private/public objects and just names of
functions. Document any important functions in the class as shown above and
they should be picked up by the parser automatically.
'''
import os, sys

########################
### Global Variables ###
########################
comm_start = '/*!'
comm_end = '*/'
# Get all possible data type combinations (not all are used)
data_types = ["bool", "char", "CItype", "double", "int", "string", "void",
"oneInt","twoInt"]
pr = ["*","&"]
data_types += ["vector<%s>" % dt for dt in data_types]
data_types += [dt+"%s" % p for dt in data_types for p in pr]
data_types += [dt+"*" for dt in data_types]

obj_types = ["class"]
func_direc = ".. cpp:function:: "
class_direc = ".. cpp:class:: "

#################
### Functions ###
#################


def SearchForCommStart(content, line, start_seq):
    for l in range(line, line + 5):
        if (len(content[l].split()) > 0 and content[l].split()[0] == start_seq):
            return True
    return False


def ParseFile(file_name, comments):
    '''
    Reads files and appends comments to comments objects.
    '''

    in_comment = False
    with open(file_name, "r") as f:
        content = f.readlines()

    for l in range(len(content)):
        # Get Function/Class Call
        if (len(content[l].split()) != 0 and
            (content[l].split()[0] in data_types or
             content[l].split()[0] in obj_types) and
                content[l].split()[-1] == "{"):

            if (content[l].count("(") - content[l].count(")") == 0):
                if (SearchForCommStart(content, l, comm_start)):
                    if (content[l].split()[0] in data_types):
                        comments.append(func_direc +
                                        content[l].translate(None, "{") + "\n\n")
                    elif (content[l].split()[0] in obj_types):
                        comments.append(class_direc +
                                        content[l].translate(None, "{")[5:] +
                                        "\n\n")

            else:
                open_count = content[l].count("(")
                close_count = content[l].count(")")
                function_call = content[l]

                # Search backwards for start of call
                for i in range(1, 20):
                    open_count += content[l - i].count("(")
                    close_count += content[l - i].count(")")

                    function_call = content[l - i] + function_call

                    if (open_count == close_count):
                        function_call = function_call.translate(None, "{")
                        function_call = function_call.translate(None, "\n")
                        if (SearchForCommStart(content, l, comm_start)):
                            comments.append(func_direc + function_call + "\n\n")
                        break

        # Get Comments
        if (len(content[l].split()) > 0 and content[l].split()[0] == comm_start):
            in_comment = True

        if (in_comment == False):
            continue
        if (in_comment == True and len(content[l].split()) == 0):
            comments.append("\n")
            continue

        if (content[l].split()[0] == comm_end):
            in_comment = False
            comments.append("\n\n")

        if (in_comment):
            comments.append(content[l])


def WriteRST(comments, input_file, prefix="../manual/", ext=".rst"):
    '''
    Writes comments to file with same name, but an rst extension.
    '''

    page_title = input_file[:-4].title() + "\n"
    for c in range(len(page_title) - 1):
        page_title += "*"
    page_title += "\n\n"

    with open(prefix + input_file[:-4] + ext, "w") as f:
        f.write(page_title)

        for line in comments:
            if (len(line.split()) == 0):
                f.write("\n")
                continue
            if (line.split()[0].strip() == comm_start):
                continue
            f.write(line)


############
### Main ###
############
if __name__ == '__main__':
    # Get files to parse
    files_to_parse = []
    for file in os.listdir(os.getcwd()):
        if file.endswith(".cpp"):
            files_to_parse.append(file)

    # Make directory for manual rst files        
    if (not os.path.exists("./../manual")):
        os.makedirs("./../manual")

    # Read/Write File(s)
    for i in range(len(files_to_parse)):
        comments = []
        print files_to_parse[i]
        ParseFile(files_to_parse[i], comments)
        WriteRST(comments, files_to_parse[i])
