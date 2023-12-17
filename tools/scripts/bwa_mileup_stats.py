#!/usr/bin/env python3
import re
import sys
import math

def nearest(x):
    return round(x * 100) / 100

def main():
    filename = sys.argv[1]

    with open(filename, 'r') as f:
        for line in f:
            pos = line.split()
            mismatch = 0
            ratio = 0

            rd = pos[4]
            rd = re.sub(r'\^\S', '', rd)
            rd = re.sub(r'\$', '', rd)

            if int(pos[3]) > 0:
                while True:
                    m = re.search(r'\+(\d+)[ATCGNatcgn]+', rd)
                    if not m:
                        break
                    l = int(m.group(1))
                    rd = re.sub(r'\+' + str(l) + r'\S{' + str(l) + r'}', '', rd)
                    mismatch += 1

                while True:
                    m = re.search(r'\-(\d+)[ATCGNatcgn]+', rd)
                    if not m:
                        break
                    l = int(m.group(1))
                    rd = re.sub(r'\-' + str(l) + r'\S{' + str(l) + r'}', '', rd)

                while len(rd) > 0:
                    i = rd[0]
                    rd = rd[1:]
                    if i in ['.', ',']:
                        continue
                    mismatch += 1

                if mismatch > int(pos[3]):
                    sys.stderr.write(f"{pos[0]}\t{pos[1]}\t{pos[2]}\t{pos[3]}\t{mismatch}\tmismatch more than coverage!\n")

                ratio = nearest(100 * (1 - mismatch / int(pos[3])))

            print(f"{pos[0]}\t{pos[1]}\t{pos[2]}\t{pos[3]}\t{mismatch}\t{ratio}")

if __name__ == "__main__":
    main()
