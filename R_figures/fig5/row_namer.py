#!/bin/python

import csv

def sample_add(file):
	with open(file) as f:
		test = csv.reader(f)
		for row in test:
			a = row
                        a.append(file)
			with open("add_samp" + file, 'a+') as new:
                        	writer = csv.writer(new)
				writer.writerows([a])

f=open("names.txt", "r")
for i in f:
	b = i
	b = b.strip()
	sample_add(b)

#f=open("guru99.txt", "a+")


#f= open("guru99.txt","w+")
