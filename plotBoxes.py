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
    # id;x;y;
    boxes = open('/home/louiseo/Documents/SCT/myrepo/splitBoxes.txt', 'r')
    for line in boxes:
        counter = 0
        for value in line.split(";"):
            #print(value)
            if(counter == 0):
                ids.append(int(value))
            if(counter == 1):
                xs.append(float(value))
            if(counter == 2):
                ys.append(float(value))
            counter +=1
    print('number of stations: ' + str(len(ids)))

    distinct_ids = list(set(ids))
    print('boxes: ' + str(distinct_ids) + ' len: ' + str(len(distinct_ids)))

    boxesN = [0 for i in range(len(distinct_ids))]
    for i in distinct_ids:
        currentN = 0
        for j in ids:
            if(i == j):
                currentN +=1
                boxesN[i] = currentN
        print('size of box: ' + str(i) + ' is: ' + str(boxesN[i]))

    x_array = np.asarray(xs)
    y_array = np.asarray(ys)

    totalSizeSoFar = 0
    for size in boxesN:
        #print(str(len(x_array[totalSizeSoFar:(totalSizeSoFar+size)])))
        #print(str(len(y_array[totalSizeSoFar:(totalSizeSoFar+size)])))
        x = np.array(x_array[totalSizeSoFar:(totalSizeSoFar+size)])
        y = np.array(y_array[totalSizeSoFar:(totalSizeSoFar+size)])
        mpl.plot(x,y,'.')
        totalSizeSoFar += int(size)

    xsct = list()
    ysct = list()
    # id;x;y;
    sctbox = open('/home/louiseo/Documents/SCT/myrepo/sctBox.txt', 'r')
    for line in sctbox:
        counter = 0
        for value in line.split(";"):
            #print(value)
            if(counter == 1):
                xsct.append(float(value))
            if(counter == 2):
                ysct.append(float(value))
            counter +=1
    print('number of stations (sct): ' + str(len(xsct)))

    xy_sct = zip(xsct,ysct)
    xy_s = zip(xs,ys)
    diff = list(set(xy_s)-set(xy_sct))
    print('diff len: ' + str(len(diff)))
    x,y = zip(*diff)
    mpl.plot(x,y,'kx')
    #xo_array = np.asarray(xsct)
    #yo_array = np.asarray(ysct)
    #mpl.plot(xo_array,yo_array,'ko')

    print('max x: ' + str(np.max(xs)) + ' min x: ' + str(np.min(xs)))
    mpl.xlim(np.max(xs),np.min(xs))
    print('max y: ' + str(np.max(ys)) + ' min y: ' + str(np.min(ys)))
    mpl.ylim(np.max(ys),np.min(ys))
    mpl.gca().set_aspect(1)
    mpl.show()

if __name__ == "__main__":
   main()
