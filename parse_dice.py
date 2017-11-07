#!/usr/bin/env python
'''
This code is meant to extract the comments from the Dice code and write them to
.rst files that can be prossessed by Sphinx.

Must be executed from the main Dice directory and will create a temporary
directory for the .rst files at the same level as the main Dice directory so
they will not be changed when the user switches branches to the gh-pages branch.

Skeleton Format Example:
'''
import os, sys

########################
### Global Variables ###
########################
comm_start = '/*!'
comm_end = '*/'
data_types = ["bool","char","CItype","double","int","string","void"]
obj_types = ["class"]
func_direc = ".. cpp:function::"

#################
### Functions ###
#################

def ParseFile(file_name, comments):
    '''
    Reads files and appends comments to comments objects.
    '''

    in_comment = False
    with open(file_name,"r") as f:
        content = f.readlines()

    for l in range(len(content)):
        # if (len(content[l].split()) == 0 and in_comment): continue

        if (len(content[l].split()) != 0 and
            content[l].split()[0] in data_types and
            content[l].split()[-1] == "{"):
            if (content[l].count("(") - content[l].count(")") == 0):
                comments.append(func_direc + content[l].translate(None,"{") +\
                 "\n\n")

            else:
                open_count = content[l].count("(")
                close_count = content[l].count(")")
                function_call = content[l]

                for i in range(1,20):
                    open_count += content[l-i].count("(")
                    close_count += content[l-i].count(")")

                    function_call = content[l-i] + function_call

                    if (open_count == close_count):
                        function_call = function_call.translate(None,"{")
                        comments.append(func_direc + function_call + "\n\n")
                        break

        if (in_comment == False): continue
        if (in_comment == True): comments.append("\n"); continue

        if (content[l].split()[0] == comm_start): in_comment = True
        if (content[l].split()[0] == comm_end): in_comment = False
        if (in_comment):
            comments.append(content[l])



def WriteRST(comments,input_file, ext=".rst"):
    '''
    Writes comments to file with same name, but and rst extension.
    '''

    page_title = input_file[:-4].title() + "\n"
    for c in range(len(page_title)-1):
        page_title += "*"
    page_title += "\n\n"

    with open(input_file[:-4] + ext, "w") as f:
        f.write(page_title)

        for line in comments:
            if (len(line.split()) == 0 ): f.write("\n"); continue
            if (line.split()[0].strip() == comm_start): continue
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

    # Read/Write File(s)
    for i in range(len(files_to_parse)):
        comments = []
        ParseFile(files_to_parse[0],comments)
        WriteRST(comments,files_to_parse[0])
