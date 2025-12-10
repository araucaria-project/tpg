#!/usr/bin/env python3
#script to generate uid / uobi

import sys
import uuid

def main():
    if len(sys.argv) < 2:
        txt =  str(uuid.uuid4())[:8]
        print(txt)
    else:
        i = 0
        n = int(sys.argv[1])
        while i < n:
            txt = str(uuid.uuid4())[:8]
            print(txt)
            i = i + 1

if __name__ == '__main__':
    main()
