# Python script to prepend a string onto the words in an input string.
import string
import sys

def prepend(line, directory):
    newline = ''
    for word in string.split(line, ' '):
        if word:
            newline = newline + directory + word + ' '
    print(newline)

if __name__ == '__main__':
    exec(sys.argv[-1])
