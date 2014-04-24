#!/usr/bin/python

import os, sys
import getopt
import random, math
import numpy as np
import os, sys
import fileinput
import time
import copy
import urllib
import logging
from collections import namedtuple, defaultdict
from math import sqrt
from heapq import heapify, heappush, heappop
from itertools import count


if sys.platform[:3] == 'win':
    timefunc = time.clock
else:
    timefunc = time.time
    
trace = lambda *args: None


def read_coords(input_data):
    Point = namedtuple("Point", ['x', 'y'])
    INF = 100000000
    best_cost = 0
    lines = input_data.split('\n')
    try:
        nodeCount = int(lines[0])
    except:
        nodeCount = 99
    points = []
    x = []
    y = []
    for i in range(1, nodeCount+1):
        line = lines[i]
        parts = line.split()
        points.append(Point(float(parts[0]), float(parts[1])))
    for i in points:
        x.append(i.x)
        y.append(i.y)
    return points, x, y


def cartesian_matrix(points):
    adj_matrix = [[0]*len(points) for i in range(len(points))]
    for i, p1 in enumerate(points):
        for j, p2 in enumerate(points):
            px = p2.x - p1.x
            py = p2.y - p1.y
            d = math.sqrt(px**2 + py**2)
            adj_matrix[i][j] = d
    print(len(adj_matrix))
    adj_matrix = np.array(adj_matrix)
    return adj_matrix


def all_pairs(size,shuffle=random.sample):
    '''generates all i,j pairs for i,j from 0-size uses shuffle to randomise (if provided)'''
    r1=range(size)
    r2=range(size)
    if shuffle:
        shuffle(r1, size-1)
        shuffle(r2, size-1)
    for i in r1:
        for j in r2:
            yield (i,j)
            

def reversed_sections(tour):
    '''generator to return all possible variations where the section between two cities are swapped
       E:
       --- The strategy used by the author in reversed_sections looks like a Branch-and-Bound?
    
    '''
    for i,j in all_pairs(len(tour)):                    
        if i != j:
            copy=tour[:]                                
            if i < j:
                copy[i:j+1]=reversed(tour[i:j+1])
            else:
                copy[i+1:]=reversed(tour[:j])
                copy[:j]=reversed(tour[i+1:])
            if copy != tour:
                yield copy


def swapped_cities(tour):
    '''generator to create all possible variations where two cities have been swapped'''
    for i,j in all_pairs(len(tour)):
        if i < j:
            copy=tour[:]
            copy[i],copy[j]=tour[j],tour[i]
            yield copy


def tour_length(matrix,tour):
    '''total up the total length of the tour based on the distance matrix'''
    tot_length=0
    num_cities=len(tour)
    for i in range(num_cities):
        j=(i+1)%num_cities
        city_i=tour[i]
        city_j=tour[j]
        tot_length+=matrix[city_i,city_j]
    return tot_length


def init_random_tour(tour_length):
   tour=list(range(tour_length))
   tour = random.sample(tour, tour_length)
   return tour


def init_notrandom_tour(tour_length):
   tour=list(range(tour_length))
   return tour


def prim(matrix, start):
        Path, Queue = {}, [(0, None, start)]
        while Queue:
            _, pulled_neighbour, u_current = heappop(Queue)
            if u_current in Path: continue
            Path[u_current] = pulled_neighbour
            for pushed_neighbour, weight in enumerate(matrix[u_current]):
                heappush(Queue, (weight, u_current, pushed_neighbour))
        return Path        


def mtsp(matrix, root):
    Tree, Cycle = defaultdict(list), []  
    for child, parent in prim(matrix, root).items():
        Tree[parent].append(child)
    def walk(root):
        Cycle.append(root)
        for visit in Tree[root]: walk(visit)
    walk(root)
    return Cycle


def run_hillclimb(init_function,move_operator,objective_function,max_iterations):
    from hillclimb import hillclimb_and_restart
    iterations,score,best=hillclimb_and_restart(init_function,move_operator,objective_function,max_iterations)
    return iterations,score,best


def timer_tsp(func, *pargs, _reps=10, **kargs):
    trace(func, pargs, kargs, _reps)
    start = timefunc()
    follow_up = []
    score = 0
    for i in range(_reps):
        iterations,score1,best = func(*pargs, **kargs)
        print(score1)
        follow_up.append((iterations,score1,best))
        score = score1+score
    score = score/_reps
    elapsed = (timefunc() - start)/_reps
    return elapsed, score, follow_up, _reps

def figures(G, points, x, y, solution):
    xs = []
    ys = []
    
    for p in solution:
        xs.append(points[p].x)
        ys.append(points[p].y)
    pylab.figure()
    pylab.plot(x,y, 'ro', markersize = 6)
    pylab.plot(xs, ys, 'b')
    
    pylab.show()


def usage():
    print("usage: python %s [-t] [-o <output image file>] [-v] [-m reversed_sections|swapped_cities] -n <max iterations> <city file>" % sys.argv[0])

def test():
    '''
    '''
    city_file = None
    input_data = None
    cmd = ("curl -L" if sys.platform == "darwin" else "wget")
    url = 'wget https://raw.githubusercontent.com/lilspikey/TSP/master/city100.txt'
    cmd += " {0}".format(url)
    if not os.path.isfile(os.getcwd()+"/"+url.split("/")[-1]):
        city_fileoriginal_download = os.popen(cmd).read()
    city_fileoriginal = os.getcwd()+'/'+url.split("/")[-1]
    os.system("sed -i 's/,/ /g' "+city_fileoriginal)
    with open(city_fileoriginal) as city_file:
        city_file = ' '.join(city_file.readlines())
        coords,x,y=read_coords(city_file)
        matrix=cartesian_matrix(coords)
        for index, test in enumerate((init_notrandom_tour, init_random_tour, mtsp)):
            if index < 2:
                init_function=lambda: test(len(coords))
            else:
                init_function=lambda: test(matrix,random.choice(list(range(len(coords)))))
            objective_function=lambda tour: -tour_length(matrix,tour)
            for move_operator in (reversed_sections, swapped_cities):
                max_iterations = 50
                logging.info('using move_operator: %s and max iterations %i'%(move_operator,max_iterations))
                print('move_operator: ',move_operator.__name__, ' '*2+test.__name__+' '*4+'-'*33)
                elapsed, score, follow_up, _reps = timer_tsp(run_hillclimb,init_function,
                                    move_operator,objective_function,max_iterations)
                print('%-9s: %.5f sec => score was %.5f after %i iterations and %i reps'%(test.__name__, elapsed, score, max_iterations, _reps))
        print()
        print()        
    if os.path.isfile(os.getcwd()+"/"+url.split("/")[-1]):
        os.remove(os.getcwd()+"/"+url.split("/")[-1])


def main():
    if sys.version_info[0] < 3:
        raise "Must be using Python 3"
        sys.exit(1)
    try:
        options, args = getopt.getopt(sys.argv[1:], "hto:vm:n:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    out_file_name=None
    max_iterations=None
    #verbose=None
    move_operator=reversed_sections
    format='%(asctime)s %(levelname)s %(message)s'
    
    for option,arg in options:
        if option == '-v':
            logging.basicConfig(format=format, level=logging.INFO)
            logging.info("verbose mode")
        elif option == '-t':
            test()
            sys.exit()
        elif option == '-h':
            usage()
            sys.exit()
        elif option == '-o':
            out_file_name=arg
        elif option == '-n':
            max_iterations=int(arg)
        elif option == '-m':
            if arg == 'swapped_cities':
                move_operator=swapped_cities
            elif arg == 'reversed_sections':
                move_operator=reversed_sections
    
    logging.basicConfig(format=format)
    logging.warning("no info!")
    if max_iterations is None:
        usage();
        sys.exit(2)
    
    if out_file_name and not out_file_name.endswith(".png"):
        usage()
        print("output image file name must end in .png")
        sys.exit(2)
    
    city_fileoriginal=args[0]
    with open(city_fileoriginal, 'r') as city_file:
        city_file = ' '.join(city_file.readlines())
        if int(city_file.split('\n')[0]) <= 100:
            try:
                coords, x, y =read_coords(city_file)
                init_function=lambda: init_random_tour(len(coords))
                matrix=cartesian_matrix(coords)
                objective_function=lambda tour: -tour_length(matrix,tour)
                logging.info('using move_operator: %s'%move_operator)
                iterations,score,best=run_hillclimb(init_function,move_operator,objective_function,max_iterations)
                print(iterations,score,best)
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
                city_file.close()
                sys.exit(2)
        else:
            print('WARNING: data is too long (len > 100) - danger of out of memory')


if __name__ == "__main__":
    main()            



