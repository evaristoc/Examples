#!/usr/bin/python
# -*- coding: utf-8 -*-


import random
import copy
import sys, os
from heapq import heapify, heappush, heappop



#### Global DATA STRUCTURES

# Generic formulation for the Capacitated Facility Location Problem
# This example was converted from the ZIMPL examples directory. See:
#
# https://code.google.com/p/python-zibopt/source/browse/trunk/examples/facility-location.py?r=170
#
# for an example to solve the same problem using zibopt library

FACILITIES = 1, 2, 3, 4
CUSTOMERS  = range(1, 10)

# Costs for opening a facility
FIXED_COST = {1:500, 2:600, 3:700, 4:800}

# Capacity of a facility at each site
CAPACITY = {1:40, 2:55, 3:73, 4:90}

# Demand from each customer
DEMAND = {1:10, 2:14, 3:17, 4:8, 5:9, 6:12, 7:11, 8:15, 9:16}

# Transportation cost from each facility to each customer
TRANSPORTATION = {                                            \
  1 : {1:55, 2: 4, 3:17, 4:33, 5:47, 6:98, 7:19, 8:10, 9: 6}, \
  2 : {1:42, 2:12, 3: 4, 4:23, 5:16, 6:78, 7:47, 8: 9, 9:82}, \
  3 : {1:17, 2:34, 3:65, 4:25, 5: 7, 6:67, 7:45, 8:13, 9:54}, \
  4 : {1:60, 2: 8, 3:79, 4:24, 5:28, 6:19, 7:62, 8:18, 9:45}  \
}




#### Global DATA STRUCTURES
Facility = {'facility_index':None, 'distcost':None}
Customer = {'index':None, 'energydemand':None, 'facilities':None}
Facility_Global = {'cindex':None, 'energy_facility':None, 'setup_cost':None}


def solve_it():
    
#### RAW DATA treatment
    inf = float("inf")
    # parse the input

    facilities_global = []    
    for i in FACILITIES:
        fg = copy.copy(Facility_Global)
        fg['index'] = i - 1 
        fg['status'] = 0
        fg['setup_cost']=FIXED_COST[i]
        fg['energy_facility'] = CAPACITY[i]
        facilities_global.append(fg)
    
    customers_list = []
    for i in CUSTOMERS:
        c = copy.copy(Customer)
        c['index']=i
        c['energydemand']=DEMAND[i]
        ftc = []
        for j in FACILITIES:
            a = copy.deepcopy(Facility)
            a['facility_index'] = j
            a['distcost'] = TRANSPORTATION[j][i]
            ftc.append(a)
        c['facilities'] = ftc
        customers_list.append(c)

    
    customer_order = []
    # OBS: customer order is a dense nested list in COLUMN-like form where rows are ordered according to energy
    for c in customers_list:
        dist_temp = []
        for d in c['facilities']:
            dist_temp.append(d['distcost'])
        customer_order.append((c['index'], dist_temp, c['energydemand']))
    customer_order = sorted(customer_order, key = lambda x: x[2])


    # values of index and energy per customer...
    customer_indexed = tuple([(c[0], c[2]) for c in customer_order])


    c_ro_temp = []
    for c in customer_order:
        c_ro_temp.append(c[1])
    # OBS: customer reorder is the distances to facilities of ordered customers, presented in a ROW-like form (see the zip(...))
    customer_reorder = []
    for cro in zip(*c_ro_temp):
        customer_reorder.append(cro)

    # Now we assign each row to a reference of each facility and re-sort the table...
    # I am including a set of those points already seen...
    c_ro_temp = []
    for ii, ff in enumerate(facilities_global):
        c_ro_temp.append((ff['index'], ff['energy_facility'], ff['setup_cost'], customer_reorder[ii], [set(), set(), set()]))
    distance_matrix = sorted(c_ro_temp, key = lambda y: (-y[1], y[2], sum(y[3])/len(y[3])))


#### FUNCTIONS        
    def system_cost_comparison(ca, pa, cind = len(customers_list)-1, dist_mat = distance_matrix):
        for v in pa:
            if v == 0:
                return
        ca_cost = 0
        ca_ffs = set(ca)
        for cc, ff in enumerate(ca):
            # find assigned customers
            if ff != 0:
                # keep track of unique facilities
                ca_cost += dist_mat[ff-1][3][cc]
            if cc == cind: break
        ca_cost_dist = ca_cost
        
        for ff in ca_ffs:
            ca_cost += dist_mat[ff-1][2]
        return ca_cost



    # a greedy-like code to fill in a facility
    def heaped_localknapsack(rowindex, current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost, dist_mat = distance_matrix):
        # if the facility already saw all...
        if len(dist_mat[rowindex][4][0] | dist_mat[rowindex][4][1] | dist_mat[rowindex][4][2]) == len(customers_list):
            # there is no reason to try again or making any update: just go to the next
            return current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost
        Q = []
        highdef = -inf
        defender_list = set()
        current = copy.deepcopy(current_status[rowindex]) #??
        # would save in the priority queue in the following order:
        # 1) a status: 0 if mine, 1 if not seen, 2 if other facility AND not lost
        # 2) negative of distance
        # 3) energy of customer
        # 4) relative index of customer (being used for all now...)
        # 5) absolute index of customer
        for ii, ff in enumerate(current_assignment):
            if dist_mat[rowindex][4][1] == set():
                # initialise with an ideal priority
                heappush(Q, (0, -dist_mat[rowindex][3][ii], customer_indexed[ii][1], ii, customer_indexed[ii][0]))
            else:
                # it's empty
                if ff == 0:
                    heappush(Q, (1, -dist_mat[rowindex][3][ii], customer_indexed[ii][1], ii, customer_indexed[ii][0]))
                # it's already mine
                elif ff == rowindex-1:
                    heappush(Q, (0, -dist_mat[rowindex][3][ii], customer_indexed[ii][1], ii, customer_indexed[ii][0]))               
                # this is OK: belongs to another facility and NOT PREVIOUSLY LOST
                elif ii not in dist_mat[rowindex][4][0]:
                    heappush(Q, (2, -dist_mat[rowindex][3][ii], customer_indexed[ii][1], ii, customer_indexed[ii][0]))
        
        # how much available capacity the facility has?
        c = dist_mat[rowindex][1]
        r = c
        findex = dist_mat[rowindex][0]+1
        def filling(Q, r, filled = False, c = dist_mat[rowindex][1]):
            my_customers = []
            while Q:
                i, d, enerdemand, relindex, cindex = heappop(Q)
                d = -d
                if enerdemand > c:
                    # save relative index of customer as too heavy to be considered for this facility
                    dist_mat[rowindex][4][2].add(relindex)
                elif enerdemand <= r:
                    r -= enerdemand
                    proposed_assignment[relindex] = rowindex+1
                    my_customers.append((d, enerdemand, relindex, cindex))
                    # saved relative index of customer as seeing
                    dist_mat[rowindex][4][1].add(relindex)
                    # if empty, occupy
                    # update for the previously empty spaces
                    if current_assignment[relindex] == 0:
                        #current[i] = (d, enerdemand, cindex)
                        current_assignment[relindex] = proposed_assignment[relindex]
                        current_status[rowindex][relindex] = (d, enerdemand, cindex)
                else:            
                    break
            return my_customers, r
        my_customers, r = filling(Q, c)
        # run until exhausting the possibilities with customers (empty Q)
        while my_customers:
            if proposed_assignment == current_assignment:
                print("equal")
                break
            # if still different...
            #... find defenders ONE by ONE
            # order the fight according to those customers the challenger facility would be willing to lose first..
            else:
                if len(my_customers) > 1:
                    my_customers = sorted(my_customers, key = lambda z: (z[0], -z[1]))
                print('num customers for facility ', rowindex, ': ',len(my_customers))
                for sc in my_customers:
                    # find the defender of that customer
                    d, enerdemand, relindex, cindex = sc 
                    defender_facility_row = current_assignment[relindex]-1
                    if defender_facility_row == rowindex: continue
                    print('a defender was found ', defender_facility_row, 'against ', rowindex,'for ', cindex)
                    #fight
                    #what we are looking for is an heuristic that take in consideration the "power" that each facility to keep the customer...
                    def_n = 0
                    cha_n = 0
                    # by correcting the following I got in trouble!!!
                    for nn in current_assignment:
                            if nn == defender_facility_row+1:
                                def_n += 1
                    for hh, nn in enumerate(proposed_assignment):
                            if nn == rowindex+1:
                                cha_n += 1
                    # finding an heuristic to make the comparisons is KEY!!
                    # Best so far...
                    ##heuristic01
                    #def_cost = (dist_mat[defender_facility_row][2]/def_n + dist_mat[defender_facility_row][3][relindex])
                    #cha_cost = (dist_mat[rowindex][2]/cha_n + dist_mat[rowindex][3][relindex]) 
                   ##heuristic02
                    def_cost = (dist_mat[defender_facility_row][2]/def_n + dist_mat[defender_facility_row][3][relindex])*dist_mat[defender_facility_row][3][relindex]/(dist_mat[defender_facility_row][2]+1)
                    cha_cost = (dist_mat[rowindex][2]/cha_n + dist_mat[rowindex][3][relindex])*dist_mat[rowindex][3][relindex]/(dist_mat[rowindex][2]+1)
                    ##heuristic03 
                    #def_cost = (dist_mat[defender_facility_row][2]/def_n + dist_mat[defender_facility_row][3][relindex])*dist_mat[defender_facility_row][3][relindex]/(dist_mat[defender_facility_row][2]+1)/def_n
                    #cha_cost = (dist_mat[rowindex][2]/cha_n + dist_mat[rowindex][3][relindex])*dist_mat[rowindex][3][relindex]/(dist_mat[rowindex][2]+1)/cha_n
                    # if winning...
                    if cha_cost < def_cost:
                        #update
                        #challenger
                        current_assignment[relindex] = proposed_assignment[relindex]
                        # evaluate how the changes goes for the totality of the system...
                        ok = False
                        for i, v in enumerate(current_assignment):
                            if v == 0: break
                            if i == len(current_assignment) - 1: ok = True
                        if ok == True:
                            if system_cost_comparison(current_assignment, proposed_assignment) < ca_base_cost:
                                ca_base = copy.deepcopy(current_assignment)
                                print('win ',ca_base_cost, system_cost_comparison(current_assignment, [1]*len(customers_list)))
                                ca_base_cost = system_cost_comparison(current_assignment, [1]*len(customers_list))
                        current_status[rowindex][relindex] = (d, enerdemand, cindex)
                        #defender
                        current_status[defender_facility_row][relindex] = 0
                        dist_mat[defender_facility_row][4][0].add(relindex)
                        defender_list.add(defender_facility_row)
                        highdef = max(highdef, defender_facility_row, rowindex)
                    # if losing...
                    else:
                        if current_status[rowindex][relindex] != 0:
                            current_status[rowindex][relindex] == 0
                        dist_mat[rowindex][4][0].add(relindex)
                        proposed_assignment[relindex] = current_assignment[relindex]
                        ok = False
                        for i, v in enumerate(current_assignment):
                            if v == 0: break
                            if i == len(current_assignment) - 1: ok = True
                        if ok == True:
                            if system_cost_comparison(current_assignment, proposed_assignment) < ca_base_cost:
                                ca_base = copy.deepcopy(current_assignment)
                                print('lose ', ca_base_cost, system_cost_comparison(current_assignment, [1]*len(customers_list)))
                                ca_base_cost = system_cost_comparison(current_assignment, [1]*len(customers_list))
                        r += enerdemand
                # continue adding customers if needed...
                my_customers, r = filling(Q, r, False)
                # if at this point, my_customers is empty, means that there is nothing else to fill in and that others need update...
                if my_customers == []:
                    break

        # after finishing all customers, update ALL defender facilities...
        print(defender_list)
        if defender_list != set():
            #for lfr in defender_list:
            for lfr in range(0, highdef+1):
                current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost = heaped_localknapsack(lfr, current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost)
        
        return current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost


####### initialising containers and variables

    current_assignment = [0]*len(customers_list)
    proposed_assignment = copy.deepcopy(current_assignment)
    current_status = [[0]*len(customers_list) for f in facilities_global]
    all_assigned = False
    not_assigning = 0
    not_assignment = [0]*len(customers_list)
    free_customers = []
    ca_base = []
    ca_base_cost = inf

###### running the function for each customer...
    for fi in range(len(facilities_global)):
        current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost = heaped_localknapsack(fi, current_assignment, proposed_assignment, current_status, ca_base, ca_base_cost)


##### printing results...
    print('best cost ==> ',system_cost_comparison(ca_base, proposed_assignment))    
    print('best case ((indexed) customers using (value) facilities: ')
    print('==> ',ca_base)


if __name__ == '__main__':
    print(solve_it())

       