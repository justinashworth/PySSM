#!/usr/bin/python

# variable and function library for PySSM

# Justin Ashworth, University of Washington, 2006
from AshworthUtil import *

dna_filter = re.compile( "[^acgtACGT]" )
ENSE_Intron = re.compile( "(ENSE|Intron)" )

bases_uppercase = ['A','C','G','T']

def sort_2D_1up(a,b):
	if a[0] < b[0]: return -1
	if a[0] > b[0]: return 1
	return 0

# sort number strings as ints
def as_int( a, b ):
	if int(a) < int(b): return -1
	if int(a) > int(b): return 1
	return 0

def pdb_prefix( filename ):
	filename = filename.split('.')
	if filename[-1] == 'pdb':
		return filename[0].split('_')
	return []

def dict_rvs( dict ):
	rvs = {}
	for key,val in dict.items():
		rvs[ val ] = key
	return rvs

def make_colors( difs, hitseq, annotations, positions, positive, negative, colormode ):

	spread = positive - negative
	colors = []
	for i in range( len(difs) ):
		colors.append( "#fff" )
		base = hitseq[i]
#		print base,
		dif = difs[i]
#		print dif,
		# translate internal unsigned index (i) to user-specific site position (integer)
		position = None
		for pos,index in positions.seqmap.items():
			if i == index:
				position = pos
				break
		if position == None:
			print 'error, site position integer not specified by internal site index'
		if annotations.has_key(position):
			if base.upper() in annotations[position] or base.lower() in annotations[position]:
				colors[i] = "#fff"

		elif colormode == "rb":

			if dif > positive: colors[i] = "#f00"
			elif dif < negative: colors[i] = "#00f"
			else:
				blue = int( 255 * (positive-dif)/spread )
				if blue < 0: blue = 0
				red = int( 255 * (dif-negative)/spread )
				if red < 0: red = 0
				red = dec2hex(   0 + red )
				blue = dec2hex(  0 + blue )
				colors[i] = "#%s00%s" % ( red, blue )
#				print colors[i]

		elif colormode == "rgb":

			if dif > positive: colors[i] = "#f00"
			elif dif < negative: colors[i] = "#0f0"
			else:
				green = int( 205 * (7*(positive-dif)/spread/4-0.75) )
				if green < 0: green = 0
				red = int(   205 * (7*(dif-negative)/spread/4-0.75) )
				if red < 0: red = 0
				blue = int(  205 - 1.5*red - 1.5*green )
				if blue < 0: blue = 0
				red = dec2hex(   50 + red )
				green = dec2hex( 50 + green )
				blue = dec2hex(  50 + blue )
				colors[i] = "#%s%s%s" % ( red, green, blue )

	return colors

def make_color_key( length, colormode ):

	color_key = []
	for i in range( length ):

		if colormode == "rb":
			blue = int( 255*( 1 - float(i)/float(length) ) )
			if blue < 0: blue = 0
			red = int( 255*( float(i)/float(length) ) )
			if red < 0: red = 0
			red = dec2hex(   0 + red )
			blue = dec2hex(  0 + blue )
			color_key.append( "#%s00%s" % ( red, blue ) )

		elif colormode == "rgb":
			green = int( 205*( 1 - 2*float(i)/float(length) ) )
			if green < 0: green = 0
			red = int(   205*( 2*float(i)/float(length) - 1 ) )
			if red < 0: red = 0
			blue = int(  205 - 1.5*red - 1.5*green )
			if blue < 0: blue = 0
			red = dec2hex(   50 + red )
			green = dec2hex( 50 + green )
			blue = dec2hex(  50 + blue )
			color_key.append( "#%s%s%s" % ( red, green, blue ) )
	return color_key

def lowercase_substitutions( seq, wt ):
	outseq = []
	for i in range(len(seq)):
		outseq.append(seq[i])
		if i <= len(wt) and seq[i].lower() != wt[i] and seq[i].upper() != wt[i]:
			outseq[i] = outseq[i].lower()
	return string.join( outseq, '' )

def make_lowercase_annotations( seq, annotations, positions ):
	newseq = ''
	for i in range( len(seq) ):
		base = seq[i]
		# translate internal unsigned index (i) to user-specific site position (integer)
#		position = positions.list[i][ base.lower() ].position
		position = None
		for pos,index in positions.seqmap.items():
			if i == index:
				position = pos
				break
		if position == None:
			print 'error, site position integer not specified by internal site index'
#		print seq,i,position
		if annotations.has_key(position):
			if base.upper() in annotations[position] or base.lower() in annotations[position]:
				newseq += base.lower()
				continue
		newseq += base
	return newseq
