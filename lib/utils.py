#!/usr/bin/env python
import os
import sys


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def welcome_msg(title, author):
	'Print Welcome message'
	print '============================================================'
	print '\t\t\t\t\t\t\t\t'
	print '\t:-) >>------->MODE-TASK<-------<< (-:\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t   :-) >>--->',title,'<---<< (-:\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t\t\t\t\t\t\t\t'
	print '\tAuthor(s):', author,'\t\t\t\t'
	print '\tResearch Unit in Bioinformatics (RUBi)\t\t'
	print '\tRhodes University, 2017\t\t\t\t'
	print '\tDistributed under GNU GPL 3.0\t\t\t'
	print '\t\t\t\t\t\t\t\t'
	print '\thttps://github.com/michaelglenister/NMA-TASK\t'
	print '\t\t\t\t\t\t\t\t'
	print '============================================================'
	return;
