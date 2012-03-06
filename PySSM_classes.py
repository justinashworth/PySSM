#!/usr/bin/python

# Justin Ashworth, University of Washington, 2006

# non-Tkinter classes used by PySSM

from PySSM_lib import *
from Tkinter import END
import re,string

# a contiguous length of basepairs with at least one score value
class Pocket:

	def __init__( self, position, seq, bkt, ddg, hbb, filename ):
		self.position = position
		self.seq = seq

		self.energies = {}
		self.energies['score'] = bkt
		self.energies['ddG'] = ddg
		self.energies['HBsp'] = hbb

		# note that 'energies' should NEVER be changed, and that ONLY 'scores' should be used for scoring.  If weighting is to be done, do it ONLY on 'scores'
		self.scores = {}
		self.scores['score'] = bkt
		self.scores['ddG'] = ddg
		self.scores['HBsp'] = hbb

		self.filename = filename

	def __str__( self ):
		short_file = self.filename.split('/')[-1]
		return '%4s %5s%10.2f%8.2f%7.2f  %s' % \
			( self.position, self.seq,
			  self.energies['score'], self.energies['ddG'], self.energies['HBsp'],
				short_file )

	def __cmp__( self, other ):
		poscmp = int(self.position) - int(other.position)
		if poscmp != 0: return poscmp
		else: return self.seq == other.seq

	# the original values should never be written over!  Use only the '_weighted' terms to change the apparent values for the pocket
	def multiply_by_weight( self, weight ):

		# can 'map' work on dictionaries?
		for term in self.scores:
			self.scores[term] = self.scores[term] * weight

	def add_weight( self, weight ):

		for term in self.scores:
			self.scores[term] = self.scores[term] + weight

	def reset_scores( self ):

		for term in self.energies:
			self.scores[term] = self.energies[term]

# contains a list of position hashes containing Pockets
# this is a normal list with some member data
# any functions that act on this list should be member functions of this class
class Positions:

	def __init__( self ):

		self.list = [] # list of dictionaries of Pocket objects at each position
		self.seqmap = {}
		# all pockets must be the same length and come from the same pdbroot
		self.pocketlength = 0
		self.pdbroot = ''

	# watch out, currently can penalize the same pocket more than once!
	# for adding single-nucleotide scores to any pocket containing the nucleotide
	def add_weight_matrix( self, weight_matrix ):

		seqmap_rvs = dict_rvs( self.seqmap )

		# find the exact pockets to change and then lookup/add by explicit key
		for position, weights in weight_matrix.items():
			for i in range( self.pocketlength ):
				offset = i - self.pocketlength + 1

				# get indices that are within range to include the weighted position
				if self.seqmap.has_key( str( int(position) + offset ) ):
					index = self.seqmap[ str( int(position) + offset ) ]
				else: continue # boundaries

				for weight in weights:
					base = weight[0]
					ipos = -1*offset
					seqs = []
					all_combinations( 0, self.pocketlength, [], seqs )
					seqs = only_combos_with( base, ipos, seqs )

					for seq in seqs:
						self.list[ index ][ seq ].add_weight( float(weight[1]) )

	def add_wt_bonuses( self, wildtype, bonus, mode = 'paranoid' ):

		# no other mode exists right now
		if mode != 'paranoid': return
		# 'paranoid' means only give bonuses to complete matches
		# the alternative would be to give non-correlative ('single-bp') bonuses

		# assumes pocket list is ordered by position
		for i in range( len(self.list) ):
			wtseq = string.join( wildtype[ i : i+self.pocketlength ], '' )
#			print i, wtseq
			self.list[ i ][ wtseq ].add_weight( bonus );

#	# if anything scored better than wildtype, give it the wildtype score
	def make_wt_best( self, wildtype, best_by = 0 ):

		for i in range( len(self.list) ):
			wtseq = string.join( wildtype[ i : i+self.pocketlength ], '' )
			wtscores = self.list[i][wtseq].scores
			for pocket in self.list[i].values():
				for term,score in wtscores.items():
					if pocket.seq == wtseq: continue
					if pocket.scores[term] < score - best_by:
						pocket.scores[term] = score - best_by

	# not currently being used
	def indices_containing( self, position ):

		seqmap_rvs = dict_rvs( self.seqmap )
		lower = int(position) + self.pocketlength - 1
		upper = int(position) - self.pocketlength + 1

		return [ i for i in range(len(self.list)) if \
		         lower >= int(seqmap_rvs[i]) and \
						 upper <= int(seqmap_rvs[i]) ]

	def invert_scores( self ):
		for i in range( len(self.list) ):
			self.list[i] = self.list[i] * -1.0

	#   -3-|6|
	#  -2-|5|
	# -1-|4|-7-
	# AAAAAAAAA
	# 123456789
	#
	# 9 positions represented by 7 overlapping pocket regions of length 3. the positions represented by pockets -1-, -7- must be upweighted by a factor 9/6, -2- and |6| by a factor 9/8, so that they contribute as much to the score as do 'saturated' pockets (-3-,|4|,|5|).
	def upweight_ends( self ):

		saturated = self.pocketlength * self.pocketlength

		for i in range( len(self.list) ):
			d = i+1
			if d >= self.pocketlength:
				if i < len(self.list) - self.pocketlength + 1: continue
				else: d = len(self.list) - i

			unsaturated = 0
			for j in range( self.pocketlength ):
				unsaturated += d
				if d != self.pocketlength: d += 1
			weight = float(saturated) / float(unsaturated)
	#		print i, d, saturated, unsaturated, weight
			for pocket in self.list[i].values():
				pocket.multiply_by_weight( weight )

	def reset_scores( self ):

		for position in self.list:
			for pocket in position.values():
				pocket.reset_scores()

# a gene, composed of separate contiguous nucleotide Segments
class Gene:

	def __init__( self, name, seq = '' ):

		self.length = 0
		self.name = name
		self.segment_seqs = []
		self.num_exons = 0
		if seq != '': self.add_seg( seq )

	def add_seg( self, seq ):
		self.segment_seqs.append( Segment( seq, self.name, self.num_exons ) )
		if self.segment_seqs[-1].length == 0:
			# Segment class found this segment to contain no sequence information
			self.segment_seqs.pop()
		else:
			self.length += self.segment_seqs[-1].length
			self.num_exons += self.segment_seqs[-1].num_exons

	def full_sequence( self ):
		full_seq = []
		separator = ['...']
		for seg in segment_seqs:
			full_seq.extend( seg.seq )
			if seg != segment_seqs[-1]: full_seq.append( separator )
		return full_seq

	def __str__( self ):

		return "%s, %s bp" % ( self.name, self.length )

# a Segment is a continuous searchable length of dna sequence, comprised of some combination of contiguous introns and exons
class Segment:

	def __init__( self, seq, gene_name, last_exon ):

		self.last_exon = last_exon
		seq = ENSE_Intron.split( seq )
		seq = filter( lambda x: regex_remove( ENSE_Intron, x ), seq )
		seq = map( lambda x: regex_filter( dna_filter, x ), seq )
		# seq is now a list of the introns and exons in the order encountered
		self.num_exons = 0
		for block in seq:
			if block == '': continue
			if block[0] in bases_uppercase: self.num_exons += 1

		exon_lab = '%i' % ( last_exon + 1 )
		if self.num_exons > 1: exon_lab += '-%i' % ( last_exon + self.num_exons )
		self.name = '%s_%s' % ( gene_name, exon_lab )

		# all pockets stored in lower-case, so the actual sequenced searched should be made lower-case in advance to allow direct lookup
		self.seq = string.lower( reduce( lambda x, y: x+y, seq ) )
		self.length = len( self.seq )

	def __str__( self ):
		return '%s  %44s  %i bp' % \
			( stringf_left( self.name, 20 ), self.seq[0:20]+'...'+self.seq[-20:],
			  len(self.seq) )


# a nucleotide subsequence with its score and info about its origin
class Hit:

	def __init__( self, gene_header, position, site, score, tag, dif_from_wt ):

		self.gene = gene_header
		self.position = position
		self.site = site
		self.score = score
		self.tag = tag
		self.dif_from_wt = dif_from_wt

	def __str__( self ):

		if self.dif_from_wt < 0: dif_from_wt = "-"
		else: dif_from_wt = str( self.dif_from_wt )

		return "%10.1f, %s, %s, %3s, %s, %8i" % \
		  ( self.score, self.site, self.tag, dif_from_wt, self.gene, self.position+1 )

	def compare_to( self, comparison ):

		if self.tag == "rvs":
			marked = mark_mismatches( rvs_comp_str(self.site), comparison )
		else:
			marked = mark_mismatches( self.site, comparison )

		if self.dif_from_wt < 0: dif_from_wt = "-"
		else: dif_from_wt = str( self.dif_from_wt )

		return "%10.1f, %s, %s, %3s, %s, %8i" % \
		  ( self.score, marked, self.tag, dif_from_wt, self.gene, self.position+1 )

	def __cmp__( self, other ):

		score_cmp = self.score - other.score
		if score_cmp != 0: return score_cmp
		else:
			sdwt = self.dif_from_wt
			odwt = other.dif_from_wt
			# negative values == undefined -> bottom of sort (instead of top)
			if sdwt < 0:
				if odwt < 0: return 0
				else: return -1
			elif odwt < 0: return 1
			else: return sdwt - odwt

# maintains an ordered list (fixed length) of the top Hits, with feedback
class HitManager:

	def __init__( self, textbox, maxhits, wildtype ):

		self.hitlist = []
		self.textbox = textbox
		self.maxhits = maxhits
		self.fail = 9999.

		self.wthit_exists = 0
		if wildtype != "":
			self.wildtype = wildtype
			self.wthit_exists = 1

		self.hitlist.append( Hit( "> dummyhit", 1, "placeholder",
		                          9999, "dum", -1 ) )
		print "\nBest hit so far:\n"
		self.textbox.insert( END, "Best hit so far:" )

	def add_hit( self, hit ):

		i = len( self.hitlist ) - 1
		if len( self.hitlist ) == self.maxhits and \
			           hit.score >= self.hitlist[i].score:
			return

		else:

			while i >= 0:
				if hit.score >= self.hitlist[i].score:
					self.hitlist.insert( i+1, hit )
					break
				elif i == 0:
					self.hitlist.insert( i, hit )
					print hit
					self.textbox.insert( END, str(hit) )
					self.textbox.see( END )
					self.textbox.update_idletasks()
				i -= 1

			if len( self.hitlist ) > self.maxhits: self.hitlist.pop()
			self.fail = self.hitlist[ len( self.hitlist ) - 1 ].score

	def add_wt_hit( self, hit ):

		self.wthit = hit
		self.wthit_exists = 1

	def display( self, textbox, comparison = [] ):

		self.textbox.delete( 0, END )
		self.textbox.insert( END, "Top %s hits:" %( self.maxhits ) )

		if self.wthit_exists == 1:
			self.textbox.insert( END, string.upper(str(self.wthit)) )

		for hit in self.hitlist:
			textbox.insert( END, hit.compare_to(comparison) )

	def write_file( self, logpath, comparison = [] ):

		logfile = open( logpath, "w" )
		seq_hdr = string.join( [' ' for i in range(len(self.hitlist[0].site)) ], '' )
		logfile.write('%10s, %s  %s, %s, %s\n' % \
		  ( 'score', seq_hdr, 'dir', 'mut', 'gene+location' ) )
		for hit in self.hitlist:
			logfile.write( hit.compare_to(comparison) + "\n" )
		logfile.close()

# (tuple) best-case early rejection with priority
class ScoreMap:

	def __init__( self, positions, scoreterm ):

		self.best_cases = []

		i = -1
		tot = 0
		for position in positions.list: # for all positions
			i += 1
			best = 9999.
			for pocket in position.values(): # for seqs (pockets) at this position
				if pocket.scores[ scoreterm ] < best: # find best
					best = pocket.scores[ scoreterm ]
			tot += best # accumulate best case scenario
			self.best_cases.append( ( i, best ) ) # remember best energy at this pos

		# sort positions so that most influential position is checked first
		self.best_cases.sort( sort_2D_1up )

		# replace best E at pos with best-case E given prioritized order
		for i in range( len( self.best_cases ) ):
			tot -= self.best_cases[i][1]
			self.best_cases[i] = ( self.best_cases[i][0], tot )

