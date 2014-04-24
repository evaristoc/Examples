import os, sys, errno
import collections
import itertools
import operator
import copy

def silentremove(filename):
# delete existing files every time it starts
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured




def creating_crovalfiles(data, eths, count):
# create new files and count the number of cases
### data is file
### eths is str
### count is int
	data.seek(0)
	count_train = 0
	count_val = 0
	count_test = 0
	training_data = open(os.getcwd()+'training.tab', 'a')
	val_data = open(os.getcwd()+'validation.tab', 'a')
	test_data = open(os.getcwd()+'test.tab', 'a')
	current_count = 0
	for line in data:
		lins = line.split('\t')
		if lins[0] == eths:
			current_count += 1
			prop = current_count/count
			if prop <= .6:
				training_data.write(line)
				count_train += 1
			if prop > .6 and prop <= .8:
				val_data.write(line)
				count_val += 1
			if prop > .8:
				test_data.write(line)
				count_test += 1
	training_data.close()
	val_data.close()
	test_data.close()
	return count_train, count_val, count_test




def creating_dict(data, eth, variables):
	'''
	- creating_dict will create several dictionaries, each with all indices of the variables as key
	- changed into a 1 sub-dictionary for all complete values (if 0 then -1; if 1 then 1) and another for blanks, instead of separated dicts (duplicated info...)
	- lose information about how many true and how many false as I group then in a sum, but save memory
	- IMPORTANT: this is a sub-level: each dictionary will be saved in a position of a list which is the value part of another super-dict which keys are the ethnicities!!!
	- see first_statistics(...) below, which transform the data and assigns each sub-dictionary to the corresponding place in the list associated to each ethnicity
	'''
### data is file
### eth is str
### variable is list
# initialise collections
	data.seek(0)
	varval_count = [-99 for yy in range(len(variables))]
	blank_count = [-99 for yy in range(len(variables))]
	cases_count = 0
#	dummy = raw_data.readline(); del(dummy)
# create a dictionary from data of the counts of 1, 0 and NA for each ethnic group
	for line in data:
		lins = line.split('\t')
		if lins[0] == eth:
			cases_count += 1
			for iter,v in enumerate(lins[1:], start = 1):
				if v == '1':
					if varval_count[iter] == -99: varval_count[iter] = 0	
					varval_count[iter] = varval_count[iter] + 1			# count cases by var
				if v == '0':
					if varval_count[iter] == -99: varval_count[iter] = 0
					varval_count[iter] -= 1
				if v == '-99':
					if blank_count[0] == -99: blank_count[0] = 0
					if blank_count[iter] == -99: blank_count[iter] = 0
					blank_count[0] = blank_count[0] + 1
					blank_count[iter] = blank_count[iter] + 1
				# if iter == 1:
					# print(eth, cases_count, varval_count[1], blank_count[1])
	data.seek(0)
# apply an imputation for of blanks when it is just 10% of data by proportionally adding 1's and 0's according to current observation
# change counts by percentages (keep a total count variable as cases_count)
	for iter in range(len(variables)):
		if (cases_count - blank_count[iter]) == 0 or cases_count == 0: continue
		if blank_count[iter] != -99:
			varval_count[iter] = varval_count[iter]/(cases_count - blank_count[iter]) # freq calculated against the number of valid cases
			blank_count[iter] = blank_count[iter]/cases_count						# freq calculated against the total of cases
		else:
			varval_count[iter] = varval_count[iter]/(cases_count)
		# if iter == 1:
			# print(eth, cases_count, varval_count[1], blank_count[1])
#	return blank_count
	return cases_count, varval_count, blank_count



def first_statistics(data, eths, variables, count_train = 0, count_val = 0, count_test = 0, create_files = False):
	'''
	- each dictionary created by creating_dict(...) will be saved in a position of a list which is the value part of another super-dict which keys are the ethnicities!!!
	- isinstance(training_eths['CEU'][2],dict)
	- (xxx_eths{eth:[count(int), f_count{variable_index:proportion}, t_count{variable_index:proportion}, blank_count{variable_index:proportion}]})
	'''
### data is file
### eths is dict
### variables is list
### count_xxx is int as keyarg
### create_files is boolean as keyarg
	eths_genfreq = {x:0 for x in eths.keys()}
	for k in eths_genfreq.keys():
		eths_genfreq[k] = creating_dict(data, k, variables)
		if create_files:
			count_train_cases, count_val_cases, count_test_cases =  creating_crovalfiles(data, k, eths_genfreq[k][0])
			count_train = count_train + count_train_cases
			count_val = count_val + count_val_cases
			count_test = count_test_cases
	return eths_genfreq, count_train, count_val, count_test

	
	
	
def creating_markers(specificity, eths_genfreq, variables):
### specificity is float
### eths_genfreq is dict
### variables is list
# variable initialisations
	y =[-99 for yy in range(len(variables))]
	eths_marker = {x:copy.copy(y) for x in eths_genfreq.keys()}
# for each sub-group of genes, create sorted list (deeper collection) based on its frequency
# use specificity to determine how specefic to that etchnic group the included genes must be
# identify those genes that are exclusively different for each ethnic group using set()
	for k in eths_marker.keys():
		count_marker_elements = 0
		for iter, var in enumerate(eths_genfreq[k][1]):
			if abs(var) >= specificity and var != -99:
				if iter < 10:
					print(k, iter, var)
				eths_marker[k][iter] = var
				if iter < 10:
					print(k, iter, eths_marker[k][iter])
				count_marker_elements += 1
		eths_marker[k][0] = count_marker_elements
		print(k, eths_marker[k][0])
	for i in eths_marker.keys():
		print(i, eths_marker[i][0])	
	return eths_marker

# for i in eths_marker.keys():
	# print(i, eths_marker[i][0])
	
# for i in eths_marker.keys():
	# for index, val in enumerate(eths_marker[i]):
		# if index < 100:
			# print(i, index, val)

			
def creating_markers2(specificity, eths_genfreq, variables):

# variable initialisations
	y =[-99 for yy in range(len(variables))]
	eths_marker = {x:copy.copy(y) for x in eths_genfreq.keys()}
# for each sub-group of genes, create sorted list (deeper collection) based on its frequency
# use specificity to determine how specefic to that etchnic group the included genes must be
# identify those genes that are different for each ethnic group based on its sign
	for k in eths_marker.keys():
		count_marker_elements = 0
		for iter, var in enumerate(eths_genfreq[k][1]):
			if abs(var) >= specificity and var != -99:
#				if iter < 10:
#					print(k, iter, var)
				eths_marker[k][iter] = var
#				if iter < 10:
#					print(k, iter, eths_marker[k][iter])
				count_marker_elements += 1
		eths_marker[k][0] = count_marker_elements
#		print(k, eths_marker[k][0])
	eths_marker_temp = copy.deepcopy(eths_marker)
	for k in eths_marker.keys():
		for l in eths_marker_temp.keys():
			if l != k:
				for iter, v in enumerate(eths_marker[k][1:], start = 1):
					if iter == 0:
						print('iter = 0')
					if v == -99:
						continue
					elif eths_marker_temp[l][iter] == -99:
						continue
					elif (v >= 0 and eths_marker_temp[l][iter] >= 0) or (v < 0 and eths_marker_temp[l][iter] < 0):
							eths_marker[k][iter] = -99
							eths_marker[k][0] -= 1
#			print(k, eths_marker[k][0])
	return eths_marker			

	
	
	
def cross_val():
	pass



def cmp_trainingvsblind(lins, marker):
# lins is list
# marker is list
# a function to compare genes of blinds that are equal to markers excluding comparisons for NA in blind
# similarity/count is my similarity function resembling a Pearson correlation: -1 for negatively correlated, 0 for no correlation, 1 for strong correlation
# also includes a error evaluation based on a range between 0 and 2: the closer to 0, the lower the error
	assert len(lins)==len(marker)
	completecases = 0
	similarity = 0
	error = 0
	x_val = 0
	for index, freq in enumerate(marker):
		if lins[index] != '-99' and marker[index] != -99:
			if lins[index] == '0': x_val = -1
			if lins[index] == '1': x_val = 1
			completecases += 1
			similarity = similarity+x_val*freq
			error = abs(freq - x_val)
	return similarity/completecases, error/completecases, similarity, error, completecases


def marker_in_blind(comparison_data, eths_marker, count_file):
	eths_comparison = {x:[] for x in range(1,count_file+1)}
	#dummy = blind_data.readline(); del(dummy)
	for k in eths_marker.keys():
		comparison_data.seek(0)
		case_number = 1
		for line in comparison_data:
			if case_number > count_file:
				break
			lins = line.split('\t')
			true_gen = 'true_gen', lins[0], cmp_trainingvsblind(lins, eths_marker[k]), eths_marker[k][0], 'CLASS', k
			eths_comparison[case_number].append(true_gen)
			case_number += 1
	return eths_comparison	


def xfrange(start, stop, step):
# a generator/iterator that travels along a real number, not along integers only (as range does)
    while start < stop:
        yield start
        start += step

#MAIN

filenames = ['training.tab', 'validation.tab', 'test.tab']
for filename in filenames:
	silentremove(os.getcwd()+'\\'+filename)


raw_data = open(os.getcwd()+'\\genestrain.tab', 'r')
eths = {}
line = raw_data.readline()
while 1:
	if not raw_data.readline(): break
	if raw_data.readline().split('\t')[0] != '':
		eths[raw_data.readline().split('\t')[0]] = [0,0,0,0]

raw_data.seek(0)

variables = 0
variables = [x for x in raw_data.readline().split('\t')]
raw_data.seek(0)
#raw_eths, count_cases_train, count_cases_val, count_cases_test = first_statistics(raw_data, eths, variables)
raw_eths, count_cases_train, count_cases_val, count_cases_test = first_statistics(raw_data, eths, variables, create_files = True)
print(count_cases_train, count_cases_val, count_cases_test)

training_data = open(os.getcwd()+'\\training.tab', 'r')
training_eths, dummy1, dummy2, dummy3 = first_statistics(training_data, eths, variables)


val_data = open(os.getcwd()+'\\genesblind.tab', 'r')
# the following dict collect the values of the chosen ethnicity...

# ... ordered according to markers that were built using a range of specificities (the for-loop)
for specificity in xfrange(.5, .60, .20):
# we first create the markers, using the training dataset...
	eths_marker = creating_markers2(specificity, training_eths, variables)
# then we compare our markers and see if the val, test or blind files have the values for a particular distribution...
	eths_val_final = marker_in_blind(val_data, eths_marker, count_cases_val)
# 
	for i in eths_val_final.keys():
		print(i)
		print(eths_val_final[i])

#raw_data.close(); training_data.close(); val_data.close()