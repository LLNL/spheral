# Python script
import string, sys

def concat(filename):
    newline = ''
    f = open(filename)
    lines = f.readlines()
    for line in lines:
	newline = newline + line[:-1]
    print(newline)

if __name__ == '__main__':
    exec(sys.argv[-1])
