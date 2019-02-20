#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pylab as mpl

def main():

    print('hello plotboxes')

    ids = list()
    xs = list()
    ys = list()
    flags = list()
    # x;y;z;t;flags;corep;pog;boxids
    boxes = open('output.txt', 'r')
    for line in boxes:
        counter = 0
        for value in line.split(";"):
            #print(value)
            if(counter == 0):
                xs.append(float(value))
            if(counter == 1):
                ys.append(float(value))
            if(counter == 4):
                flags.append(int(value))
            if(counter == 7):
                ids.append(int(value))
            counter +=1
    print('number of stations: ' + str(len(ids)))

    distinct_ids = list(set(ids))
    print('boxes: ' + str(distinct_ids) + ' len: ' + str(len(distinct_ids)))

    x = np.array(xs)
    y = np.array(ys)
    ids = np.array(ids)
    boxesN = [0 for i in range(len(distinct_ids))]
    for i in range(len(distinct_ids)):
        currentN = 0
        for j in ids:
            if(distinct_ids[i] == j):
                currentN +=1
                boxesN[i] = currentN
        print('size of box: ' + str(i) + ' is: ' + str(boxesN[i]))
        print(distinct_ids[i])
        I = np.where(ids==distinct_ids[i])
        print(I)
        mpl.plot(x[I],y[I],'.')

    print('max x: ' + str(np.max(xs)) + ' min x: ' + str(np.min(xs)))
    mpl.xlim(np.max(xs),np.min(xs))
    print('max y: ' + str(np.max(ys)) + ' min y: ' + str(np.min(ys)))
    mpl.ylim(np.max(ys),np.min(ys))
    mpl.gca().set_aspect(1)
    mpl.show()

if __name__ == "__main__":
   main()
