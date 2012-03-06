#!/usr/bin/env python
# version 1 2009-07-28
# vi: set wrap:

### guidelines ###
# Make new windows inherit from Toplevel (derived Toplevel classes)
# use 'self' for class and function vars/objs ONLY if necessary

# This is a position-specific search matrix (PSSM) sequence-searching program written in Python by Justin Ashworth (ashwortj@u.washington.edu), University of Washington, 2006, with scoring functions carried over from a perl script version written previously by Justin Ashworth and Umut Ulge.

# The general strategy is to search on a sliding window basis, prioritizing the first position searched in each window by the relative weight of the score for that position.  i.e., the positions that are worth most are checked first, in a sequence-discontinuous manner.  This allows the early rejection of windows that are guaranteed to fail after a given number (n<length) positions checked.  The rejection threshold is a function of the currently worst hit in a capped list of hits, and the algorithm speeds up as the average score of the best hits increases, as this leads to earlier rejection.  Thus, a great number of irrelevant pairwise comparisons are avoided.  The relative efficiency of this approach vs. a dynamic programming one is debatable.  A Needleman-Wunch-Smith-Waterman style approach would calculate/lookup PSSM pairwise comparisons only once for each pair, but many of these comparisons will be irrelevant.  However the sliding window repeatedly looks up (but does not recalculate) some of the same pair scores many times.

# The scores are derived from real or theoretical free energies, where negative energies are favorable.  The direct use of energy as the match score is preferable to the use of probability for two reasons: 1) the raw data are usually energies, and 2) energies are additive, while probabilities are multiplicative.  The use of a favorable negative score for matches presents no apparent limitations on the search algorithm.

# The total energy (score) of a site is composed of summed and/or averaged energy terms for each position or subset of positions in the site.  The energy term(s) associated with any identity (i.e. one nucleotide) or combinations of identities (i.e. a pocket of 3 nucleotides) in the site are stored are looked up from the user's data (a comprehensive position-specific weight matrix).

# On scoring using overlapping pockets: Nucleotide specificity is definitely not divisible into independent single identities - neighbors matter (there is even experimental evidence of this).  Therefore since specificity at adjacent nucleotides is likely to covariate, it useful to model and score covariate substitutions.  Scoring by more than one nucleotide at a time allows correlative scoring.  The way this is currently handled is that pockets of length n are scored at each position in the query sequence, allowing overlap with the assumption that scoring each nucleotide more than once is okay, as long as every nucleotide is scored the same number of times (which will happen to be equal to the pocket length).  If the first pocket scored is from positions 1 to pocketlength (and the last ends on sitelenth-pocketlength+1) then the scores must be upweighted for positions that are less than 'pocketlength' from the beginning or end of the query, because these positions will be scored less times than the rest.

# Pmw 'megawidgets' (and all other non-standard modules) are intentionally avoided in this code to prevent dependency issues.

# Note: Python 2.5 or newer is strongly recommended, as previous versions have problems releasing allocated memory once it is no longer needed.

# Moreover not all of the memory problems are fixed in this program yet.  For example, memory allocation associated with objects/functions containing Toplevel() widgets is not deallocated when the Toplevel window is closed.  Fortunately the memory IS released back to the system when the program is closed.

import os, re, operator, string, glob
from Tkinter import *
import tkFont
import FileDialog

# basic variables and functions are in PySSM_lib.py
from PySSM_lib import *
# non-Tkinter class are in PySSM_classes.py
from PySSM_classes import *

from optparse import OptionParser
p = OptionParser()
p.add_option( '-s', '--script', help='script mode (no GUI)' )
p.add_option( '-t', '--term', help='score term (key)', default='bk_tot' )
p.add_option( '-i', '--invert', help='invert scores (currently only works on loaded pssm files', default=False, action='store_true', dest='invert' )
options,args = p.parse_args()

#pre-compiled regular expressions
re_pdb = re.compile( "pdb" )
re_gz = re.compile( "gz" )
re_energies = re.compile( "energies" )
re_hit = re.compile( "hit" )

class App:

	def __init__( self, master ):

		# global vars
		# should eliminate these whenever logically sound
		# many globals can be left undefined, so that __dict__ can be used to check for their existence (hence initialization)
		self.bestfiles = {}
		self.pdbroot = ""
		self.wildtype = ""
		self.annotations = {}
		self.geneswindows = []

		self.maxhits = IntVar()
		self.maxhits.set( 100 )
		self.maxhits_graphic = 10
		self.master = master

		self.delim = StringVar()
		self.delim.set( ":" )

		self.colormode = StringVar()
		self.colormode.set( "rb" )

		self.invert_scores_flag = IntVar()
		self.invert_scores_flag.set( 0 )
		self.upweight_ends_flag = IntVar()
		self.upweight_ends_flag.set( 1 )

		self.wt_bonus_flag = IntVar()
		self.wt_bonus_flag.set( 0 )
		self.wt_bonus = DoubleVar()
		self.wt_bonus.set( "-2" )
		self.wt_best_flag = IntVar()
		self.wt_best_flag.set( 0 )
		self.wt_best = DoubleVar()
		self.wt_best.set( "-0.5" )

		self.init_cutoffs()

		self.temp = 1 # Boltzman temperature

		master.configure( background = "#ddd" )
		mainmenu = Titlemenu( master, self )
		self.mainframe = Frame( master, background = "#ddd", width = 50 )
		self.buttonframe = Frame( master, background = "#ddd" )
		self.mainframe.grid(   row = 0, column = 0 )
		self.buttonframe.grid( row = 1, column = 0, sticky = E+W, ipady = 5 )
		self.main_widgets()

	def main_widgets( self ):

		self.status = Label( self.mainframe, text = "",
		                     background = "#ddd",
		                     justify = LEFT )

		self.message = Text( self.mainframe, height = 20, width = 60,
		                     background = "#fff" )
		messagesb = Scrollbar( self.mainframe, orient = VERTICAL )
		self.message.configure( yscrollcommand = messagesb.set )
		messagesb.configure( command = self.message.yview )

		self.status.grid(  row = 0, column = 0, sticky = W )
		self.message.grid( row = 1, column = 0 )
		messagesb.grid( row = 1, column = 1, sticky = N+S )

		bclear = Button( self.buttonframe,
		                 text = "Clear",
		                 highlightthickness = 0,
		                 command = self.global_reset )

#		quicklabel = Label( self.buttonframe,
#		                    text = "Quickstart:",
#		                    background = "#ddd" )

#		self.energyfile = File_match_menu( self.buttonframe, 'score', self, 24, True )

		bclear.grid(          row = 0, column = 0, sticky = W )
#		quicklabel.grid(      row = 0, column = 1, sticky = E )
#		self.energyfile.grid( row = 0, column = 3, sticky = E )

	def demo( self ):

		self.get_scores( 'file', 'msoB.energies' )
		self.get_from_file( 'matrix', 'msoB.cen4' )
		self.get_from_file( 'wt', "msoB.wt" )
		geneswindow = GenesWindow( self, "XSCID.seq" )
		self.geneswindows.append( geneswindow )

	def create_data_window( self ):

		self.dataw = Toplevel( self.master, background = "#ddd" )
		self.dataw.title = "Data window" # ineffective
		dataw_head = Frame( self.dataw, background = "#ddd" )
		dataw_head.grid( row = 0, column = 0, sticky = W )
		dataw_box = Frame( self.dataw, background = "#ddd" )

		# sorting radio buttons
		# this variable also controls which energy term is used during search
		self.sortkey = StringVar()
		sort_tot = Radiobutton( dataw_head, text = "tot",
		                        variable = self.sortkey, value = "score",
		                        background = "#ddd",
		                        highlightthickness = 0,
		                        activebackground = "#ddd",
		                        activeforeground = "#257",
		                        command = self.energy_show )

		sort_ddg = Radiobutton( dataw_head, text = "DDG",
		                        variable = self.sortkey, value = "ddG",
		                        background = "#ddd",
		                        highlightthickness = 0,
		                        activebackground = "#ddd",
		                        activeforeground = "#257",
		                        command = self.energy_show )

		sort_hbb = Radiobutton( dataw_head, text = "hb",
		                        variable = self.sortkey, value = "HBsp",
		                        background = "#ddd",
		                        highlightthickness = 0,
		                        activebackground = "#ddd",
		                        activeforeground = "#257",
		                        command = self.energy_show )
		sort_tot.select()

		self.sortmode = StringVar()
		modes = ['global','position','unique']
		sortmode_chks = []
		for mode in modes:
			sortmode_chks.append( Radiobutton( dataw_head, text = mode,
		                        variable = self.sortmode, value = mode,
		                        background = "#ddd",
		                        highlightthickness = 0,
		                        activebackground = "#ddd",
		                        activeforeground = "#257",
		                        command = self.energy_show )
			                     )
		sortmode_chks[1].select()

		sort_tot.grid( row = 0, column = 0, ipadx = 8 )
		sort_ddg.grid( row = 0, column = 1, ipadx = 8 )
		sort_hbb.grid( row = 0, column = 2, ipadx = 8 )
		column = 2
		for chk in sortmode_chks:
			column += 2
			chk.grid( row = 0, column = column, ipadx = 8 )

		self.tbox1 = Listbox( dataw_box, font = "Courier -12",
		                      height = 40, width = 60, background = "#fff" )
		tbox1sb = Scrollbar( dataw_box, orient = VERTICAL )
		self.tbox1.configure( yscrollcommand = tbox1sb.set )
		tbox1sb.configure( command = self.tbox1.yview )
		self.tbox1.bind( "<Double-Button-1>", self.dblclk_pymol )

		dataw_box.grid(  row = 1, column = 0 )
		self.tbox1.grid( row = 0, column = 0, ipadx = 5, ipady = 5 )
		tbox1sb.grid(    row = 0, column = 1, ipadx = 2, sticky = N+S )

	def test_bind_listbox( self, event ):

		self.message.insert( END, self.tbox1.get( ACTIVE ) + "\n" )

	def dblclk_pymol( self, event ):

		row = self.tbox1.get( ACTIVE )
		keys = row.split()
		file = self.bestfiles[ self.pdbroot + "_" + keys[0] + "_" + keys[1] ]
		if os.path.exists( file ):
			print file
			self.message.insert( END, "loading file %s in PyMOL.\n" % file )
			# the following doesn't work yet
			if sys.platform=='win32':
				os.system( "c:\\Program Files\\DeLano Scientific\\PyMOL\\PyMOL.exe %s" % file )
			else: os.system( "pymol %s" % file )
		else:
			self.message.insert( END,
			                     "File %s not found in current directory!" % file )

	def create_settings_dialog( self, master ):

		settingsw = Toplevel( self.master, background = "#ddd" )

		rownum = 0
		maxlab = Label( settingsw, text = "# Hits:",
		                     background = "#ddd" )
		usermaxhits = Entry( settingsw, width = 10,
		                     text = self.maxhits.get(),
		                     textvariable = self.maxhits,
		                     background = "#fff" )
		maxlab.grid( row = rownum, column = 0 )
		usermaxhits.grid( row = rownum, column = 1, columnspan = 2 )
		rownum += 1

		hdrlab = Label( settingsw, text = "Delimiter:",
		                     background = "#ddd" )
		userheaderdelim = Entry( settingsw, width = 10,
		                         text = self.delim.get(),
		                         textvariable = self.delim,
		                         background = "#fff" )
		hdrlab.grid( row = rownum, column = 0 )
		userheaderdelim.grid( row = rownum, column = 1, columnspan = 2 )
		rownum += 1

#		inv_scrs_chk = Checkbutton( settingsw, background = "#ddd",
#		                            highlightthickness = 0,
#		                            activebackground = "#ddd",
#		                            text = "Invert scores",
#		                            variable = self.invert_scores_flag )

#		inv_scrs_chk.grid( row = rownum, columnspan = 3 )
#		rownum += 1

		cutlb = Label( settingsw, text = "Energy Cutoffs (color):",
		               background = "#ddd" )
		cutlb.grid( row = rownum, columnspan = 3 )
		rownum += 1

		cutentries = []
		for scoreterm, cutoff in sorted( self.cutoffs.items() ):
			cutentries.append( Label( settingsw, text = scoreterm,
			                          background = "#ddd" ) )
			cutentries[-1].grid( row = rownum, column = 0 )
			for i in range(2):
				val = cutoff[i]
				cutentries.append( Entry( settingsw, width = 6, text = val.get(),
		                              textvariable = val, background = "#fff" ) )
				cutentries[-1].grid( row = rownum, column = i+1 )
			rownum += 1

		upw_ends_chk = Checkbutton( settingsw, background = "#ddd",
		                            highlightthickness = 0,
		                            activebackground = "#ddd",
		                            text = "Upweight ends",
		                            variable = self.upweight_ends_flag )

		upw_ends_chk.grid( row = rownum, columnspan = 3 )
		rownum += 1

		wt_bonus_chk = Checkbutton( settingsw, background = "#ddd",
		                            highlightthickness = 0,
		                            activebackground = "#ddd",
		                            text = "WT bonus of",
		                            variable = self.wt_bonus_flag )
		wt_bonus_chk.grid( row = rownum, columnspan = 2 )
		self.wt_bonus_entry = Entry( settingsw, width = 4,
		                             text = self.wt_bonus.get(),
		                             textvariable = self.wt_bonus,
		                             background = "#fff" )
		self.wt_bonus_entry.grid( row = rownum, column = 2 )
		rownum += 1

		wt_best_chk = Checkbutton( settingsw, background = "#ddd",
		                           highlightthickness = 0,
		                           activebackground = "#ddd",
		                           text = "WT best by",
		                           variable = self.wt_best_flag )
		wt_best_chk.grid( row = rownum, columnspan = 2 )
		self.wt_bonus_entry = Entry( settingsw, width = 4,
		                             text = self.wt_best.get(),
		                             textvariable = self.wt_best,
		                             background = "#fff" )
		self.wt_bonus_entry.grid( row = rownum, column = 2 )

		rownum += 1

		# make this a radio menu
#		self.clrmode = Label( settingsw, text = "Color mode:" )
#		self.usercolormode = Entry( settingsw, text = self.colormode.get(),
#		                              textvariable = self.colormode )

		settings_exit = Button( settingsw, text = "done",
		                        command = settingsw.destroy )

		settings_exit.grid( row = rownum, columnspan = 3 )

	def global_reset( self ):

		self.window_reset()
		self.energy_reset()

	def window_reset( self ):
		# unfortunately I still can't get memory deallocation from destruction of the Toplevel classes

		for geneswindow in self.geneswindows:
			geneswindow.destroy()
			del geneswindow
		self.geneswindows = []

		if self.__dict__.has_key( 'dataw' ): self.dataw.destroy()

	def energy_reset( self ):

		self.pdbroot = ""
#		self.energyfile['state'] = NORMAL

	def loadfile( self, master ):

		dialog = FileDialog.LoadFileDialog( master, title = "Load score file" )
		filename = dialog.go()
		if filename != None:
			self.get_scores( 'file', filename )

	def loadpath( self, master, mode = '' ):

		dialog = FileDialog.FileDialog( master, title = "Load pdb path" )
		filepath = dialog.go()
		if filepath != None:
			if not os.path.isdir( filepath ):
				self.status['text'] = "%s is not a directory!" % filepath
				print "%s is not a directory!" % filepath
			else: self.get_scores( mode, filepath )

	def loadpssm( self, master ):

		dialog = FileDialog.LoadFileDialog( master, title = "Load PSSM" )
		filename = dialog.go()
		if filename != None:
			self.get_scores( 'pssm', filename )

	def loadmatrix( self, master ):

		dialog = FileDialog.LoadFileDialog( master, title = "Load matrix" )
		filename = dialog.go()
		if filename != None:
			self.load_matrix( filename )

	def loadannotations( self, master ):

		dialog = FileDialog.LoadFileDialog( master, title = "Load lowercase annotations" )
		filename = dialog.go()
		if filename != None:
			self.load_lowercase_annotations( filename )

	def loadseq( self, master ):
		dialog = FileDialog.LoadFileDialog( master, title = "Load sequence file" )
		filename = dialog.go()
		if filename != None:
			self.geneswindows.append( GenesWindow( self, filename ) )

	def loadwt( self, master, app ):
		dialog = FileDialog.LoadFileDialog( master, title = "Load wildtype file" )
		filename = dialog.go()
		if filename != None:
			app.get_from_file( 'wt', filename )

	def load_matrix( self, source ):
		file = open( source, 'r' )
		type = file.readline().split()[1]
		file.close()

		self.aux_matrix = load_simple_matrix( source )
		self.message.insert( END, 'loaded auxiliary weight matrix from %s\n'
		                     % source )
		for position, weights in self.aux_matrix.items():
			self.message.insert( END, position + str(weights) + '\n' )

	def get_scores( self, mode, source, term = 'bk_tot' ):

		self.energy_reset()
		self.positions = Positions() # postponed definition of object global
		if re.search('pssm$',source): mode = 'pssm'

		if mode == 'file':
			# for quickstart, load the first "energies" file found in .
			if source == "quickfile":
				for file in os.listdir( "." ):
					if re_energies.search( file ):
						source = file
						break

			if os.path.exists( source ):
				print "Getting previously-tabulated energies...\n"
				self.message.insert( END, "opening score file: %s...\n" % source )
				print "opening %s\n" % source

				get_energies_from_file( source, self.positions, self.bestfiles )

		elif mode == 'pssm':
			if os.path.exists( source ):
				print "Loading PSSM...\n"
				self.message.insert( END, "opening pssm file: %s...\n" % source )
				get_energies_from_pssm( source, self.positions )

		else: # mode == 'pdb' or mode == '' or mode == 'best'

			if not os.path.exists( source ): source = '.'
			self.message.insert( END, "opening pdb path: '%s' ...\n" % source )

			get_energies_from_pdbs( source, self.positions, self.bestfiles, term, mode )

			if len(self.positions.list) == 0:
				self.message.insert( END, "Failed to load any pdb's from %s!" % source )
				return

		self.pdbroot = self.positions.pdbroot
#		self.energyfile['state'] = DISABLED

		self.create_data_window()
		self.energy_show()

		self.status['text'] = "Energies obtained"

		if mode == 'file': return
		rawfile = "energies"
		if self.pdbroot != '': rawfile = self.pdbroot + '.' + rawfile
		write_raw( self.positions, rawfile )

	def load_lowercase_annotations( self, filename ):
		for line in open(filename):
			if line.startswith('#'): continue
			s = line.split()
			if len(s) < 2: continue
			position = s[0]
			base = s[1]
			if not self.annotations.has_key(position): self.annotations[position] = []
			self.annotations[position].append(base)

		print 'loaded lowercase annotations:'
		for p,bases in self.annotations.items():
			print p,
			for b in bases: print b,
			print

	def energy_show( self ):

		self.tbox1.delete( 0, END )

		allpockets = []
		for position in self.positions.list:
			allpockets.extend( position.values() )
		allpockets.sort( lambda x,y: sort_pockets( x, y, self.sortkey.get(),
		                                           self.sortmode.get() ) )
		for pocket in allpockets:
			self.tbox1.insert( END, '%s' % str( pocket ) )

	def get_from_file( self, type, filename ):
		if type == 'score': self.get_scores( 'file', filename )
		if type == 'wt': self.get_wt( filename )
		elif type == 'matrix': self.load_matrix( filename )
		elif type == 'seq':
			self.geneswindows.append( GenesWindow( self, filename ) )
		elif type == 'annotations': self.load_lowercase_annotations( filename )

	def get_wt( self, filename ):

		self.message.insert( END, "opening wildtype file: %s...\n" % filename )
		file = open( filename, "r" )

		for line in file:
			self.wildtype = self.validate_site( string.strip( line ) )
			break # just one line for now

		self.message.insert( END, str( len(self.wildtype) ) +
		                     " bp wildtype site loaded: " + self.wildtype + "\n" )

	def validate_site( self, inseq ):

		site = ""
		for i in range( len(inseq) ):
			char = string.lower( inseq[i] )

			for i in range( len(nucs) ):
				if char == nucs[i]:
					site += char
					break
				if i == len(nucs) - 1:
					self.message.insert( END,
					                     "Error: unknown bases in 'wildtype' site.\n" )
					return ""
		return site

	# package cutoff values into tk DoubleVars to allow user to set them
	def init_cutoffs( self ):

		cutlist = [ ( 'score', 0.0, 1.0 ),
		            (  'ddG', -3.0, 6.0 ),
		            ( 'HBsp', -1.0, 2.0 ) ]

		self.cutoffs = {}
		for cut in cutlist:
			self.cutoffs[ cut[0] ] = ( DoubleVar(), DoubleVar() )
			self.cutoffs[ cut[0] ][0].set( cut[1] )
			self.cutoffs[ cut[0] ][1].set( cut[2] )

### END OF APP CLASS ###

### BEGIN TKINTER WINDOW CLASSES ###
# all new windows (Toplevels) must be classes and belong here

# link destruction of this class to the emptying of the Genes_list that led to its creation?  Probably want keep multiple copies of Genes_list in App too
class GenesWindow( Toplevel ):

	def __init__( self, app, filename ):
		Toplevel.__init__( self, app.master, background = "#ddd" )

		self.genestext = Listbox( self, font = "Courier -12",
		                          height = 20, width = 80, background = "#fff" )
		genestextsb = Scrollbar( self, orient = VERTICAL )
		self.genestext.configure( yscrollcommand = genestextsb.set )
		genestextsb.configure( command = self.genestext.yview )

		self.search = Button( self, text = "Search using current matrix",
		                      highlightthickness = 0,
		                      command = lambda:
		                      self.genes_search( app, self.genes ) )

		self.genestext.grid( row = 0, column = 0, ipadx = 5, ipady = 5 )
		genestextsb.grid(    row = 0, column = 1, ipadx = 2, sticky = N+S )
		self.search.grid(    row = 1, column = 0, columnspan = 2, ipady = 3 )

		self.genes, self.bps = self.load_genes( app, filename )

		self.filename = filename

	def load_genes( self, app, filename ):

		app.message.insert( END, "opening gene file: %s...\n" % filename )

		files = [ filename ]
		genes = []

		# list of filenames to load
		if filename.split('.')[-1] == 'list':
			filepath = string.join( filename.split('/')[:-1], '/' )
			files = []
			file = open( filename, 'r' )
			for line in file:
				files.append( filepath + '/' + line.strip() )

		for filename in files:
			if os.path.isdir( filename ): continue
			if not os.path.exists( filename ):
				app.message.insert( END, '%s does not exist!' % filename )
				continue
			print filename
			file = open( filename, "r" )

			header = ''
			seq = ''

			dot = filename.split('.')
			file_ext = dot[-1]

			# 'special' format: text file cut-and-pasted from the Ensembl 'exon view.'  This forms a searchable segment (presumably w/ some intronic sequence flanking or connecting exons) for each run of sequence encountered between '.....' discontinuities.
			if file_ext == 'ens':
				# add any words that contain a,c,g,t,A,C,G,T to the word filter
				word_filter = \
					re.compile( '(upstream sequence|downstream sequence)' )
				# here note the inclusion of the dot (retains contig boundaries)
				genes.append( Gene( dot[0].split('/')[-1] ) )
				# read entire file into a string, and filter out some words
				segments = word_filter.sub( '', file.read() ).split('.')
				# let the Gene/Segment classes do the rest of the filtering
				for segment in segments:
					if segment == '' or segment [0:2] == '5\'': continue
					genes[-1].add_seg( segment )

			# FASTA-style
			else:
				for line in file:
					if line[0] == ">":
						if len(seq) != 0:
							genes.append( Gene( header, seq ) )
							header = ""
							seq = ''
						header = string.strip( ( line.split( app.delim.get() ) )[0] )
					else:
						line = line.strip() # remove spaces
						# append only nucleotide letters
						dna_filter.sub( '', line )
						seq += line

				# get the last one
				if len(seq) != 0: genes.append( Gene( header, seq ) )
			file.close()

		bps = 0
		for gene in genes:
			self.genestext.insert( END, gene )
			for seg in gene.segment_seqs:
				self.genestext.insert( END, str(seg) )
			bps += gene.length
			self.genestext.insert( END, '' )
		app.message.insert( END, "Gene file %s loaded: %i bp\n" % (filename,bps) )
		return genes, bps

	def genes_search( self, app, genes ):

		self.hitsarrays = []
		app.positions.reset_scores()

		if len(app.positions.list) == 0:
			app.status['text'] = "Error: no data exists for scoring!"
			return

		# WEIGHTING HAPPENS HERE ( never earlier! and never before scoremap! ). Weights are NEVER applied to original 'energy' values, only to 'scores', and only by using member functions of the Positions class

		# invert scoring matrix (for matrices in which positive scores are favorable)
		if app.invert_scores_flag.get() == 1:
			app.positions.invert_scores()

		# add weighted terms to every pocket
		if app.upweight_ends_flag.get() == 1:
			app.positions.upweight_ends()

		# add auxiliary matrix if present
		if app.__dict__.has_key( 'aux_matrix' ):
			app.positions.add_weight_matrix( app.aux_matrix )

		if app.wt_bonus_flag.get() == 1:
			app.positions.add_wt_bonuses( app.wildtype, app.wt_bonus.get() )
		if app.wt_best_flag.get() == 1:
			app.positions.make_wt_best( app.wildtype, app.wt_best.get() )

		# not all sorting terms correspond to scoring terms
		scoreterm = scorekey( app.sortkey.get() )

		app.message.insert( END, "searching genes...\n" )
		app.message.insert( END, "scoring with score term:  >>> %s <<<\n"
		                     % scoreterm )
		app.message.insert( END, "Progess: " )
		app.message.update_idletasks()

		self.search['state'] = DISABLED

		# keep these local to avoid memory leaks
		hits = HitManager( self.genestext, app.maxhits.get(), app.wildtype )
		# ScoreMap creation handles scoring priority setup here
		scoremap = ScoreMap( app.positions, scoreterm )

		sitelength = len(app.positions.list) + app.positions.pocketlength - 1
		count = 0
		count_report = self.bps / 10
		perc = 10
		for gene in genes:
			print gene.name
			for seg in gene.segment_seqs:
				self.genestext.insert( END, seg.name )
				self.genestext.see( END )
				self.genestext.update_idletasks()

				for position in range( len(seg.seq) - sitelength + 1 ):
					count += 1
					# pockets are all stored in lowercase!
					site = seg.seq[ position : position + sitelength ]

					score_fwd = site_score_fwd( site, app.positions,
					                            app.positions.pocketlength,
					                            scoremap, hits.fail, scoreterm )
					if score_fwd != "fail":
						hits.add_hit( Hit( seg.name, position, site, score_fwd, "fwd",
						                   diff_fwd( site, app.wildtype ) ) )

					score_rvs = site_score_rvs( site, app.positions,
					                            app.positions.pocketlength,
					                            scoremap, hits.fail, scoreterm )
					if score_rvs != "fail":
						hits.add_hit( Hit( seg.name, position, site, score_rvs, "rvs",
						                   diff_rvs( site, app.wildtype ) ) )

					if count == count_report:
						count = 0
						app.message.insert( END, '%i%s, ' % ( perc, "%" ) )
						self.genestext.see( END )
						perc += 10
						app.message.update_idletasks()

		app.message.insert( END, "done.\n" )

		# score and store the wildtype if it exists
		# skip optimizations and store energies by position
		if app.wildtype != "":
			score_wt, app.scores_wt = site_score_complete( app.wildtype,
			                          app.positions, app.positions.pocketlength,
			                          scoreterm )

			hits.add_wt_hit( Hit( "*** Wildtype Site ***", 0, app.wildtype,
			                      score_wt, " wt", 0 ) )

		hits.display( self.genestext, app.wildtype )

		logpath = self.filename.split("/")[-1]
		logpath = self.filename.split("\\")[-1]
		logpath = app.pdbroot + '.' + logpath + '.hits.' + scoreterm
		hits.write_file( logpath, app.wildtype )

		self.hitsarrays.append( HitsArrayWindow( self, app, hits, scoreterm ) )

		self.search['state'] = NORMAL

# inheritance from Toplevel is good practice, but I still see memory deallocation issues upon closure of the window.  Instances of this class may still be poorly managed?
class HitsArrayWindow( Toplevel ):

	def __init__( self, master, app, hits, scoreterm ):
		Toplevel.__init__( self, master )

		self.colorson = IntVar()
		self.colorson.set(1);
		self.hitstart = 0
		self.maxdisplay = app.maxhits_graphic
		if len(hits.hitlist) <= self.maxdisplay:
			self.maxdisplay = len(hits.hitlist)

		self.configure( background = "#ddd" )
		negative, positive = lookup_cutoffs( app.cutoffs, scoreterm )
		color_key = make_color_key( len(app.wildtype), app.colormode.get() )

		rownum = 0
		lab_term = Label( self, text = "TERM:", font = smallfixed, background='#ddd' )
		lab_term.grid( row = rownum, column = 0, sticky = W )
		lab2_term = Label( self, text = scoreterm, font = bigfixed, background='#ddd' )
		lab2_term.grid( row = rownum, column = 1 )

		# yay for lisp
		wt_score = reduce( lambda x, y: x + y, app.scores_wt )

		rownum += 1
		lab_wt = Label( self, text = "%6.1f  " % wt_score, background='#ddd' )
		lab_wt.grid( row = rownum, column = 0, sticky = S )
		array_wt = WTArray( self, app.wildtype, app.positions.seqmap,
		                    app.positions.pocketlength )
		array_wt.grid( row = rownum, column = 1 )
		lab2_wt = Label( self, text = "Wildtype", background='#ddd' )
		lab2_wt.grid( row = rownum, column = 2, sticky = S+W )

		# the hits
		self.hitarrays = HitArrays( self, app, hits, scoreterm, 0, self.maxdisplay,
		                            negative, positive, rownum )
		rownum = self.hitarrays.endrow + 1

		lab_key = Label( self, text = "KEY:", font = smallfixed, background='#ddd' )
		lab_key.grid( row = rownum, column = 0, sticky = E )
		array_key = KeyArray( self, len(app.wildtype), color_key,
		                      negative, positive )
		array_key.grid( row = rownum, column = 1  )
		lab_key = Label( self, text = "%4.1f to %4.1f" % (negative,positive), background='#ddd' )
		lab_key.grid( row = rownum, column = 2, sticky = W )

		rownum += 1
		self.prev = Button( self, text = "Prev", highlightthickness = 0,
		                    command = lambda: self.browse( 0, app, hits, scoreterm ),
		                    state = DISABLED )
		self.prev.grid( row = rownum, column = 0, sticky = E )

		self.rangelab = Label( self, background = "#ddd",
		                       text = "showing hits %i to %i" %
		                       ( self.hitstart, self.hitstart + self.maxdisplay ) )
		self.rangelab.grid( row = rownum, column = 1 )

		self.next = Button( self, text = "Next", highlightthickness = 0,
		               command = lambda: self.browse( 1, app, hits, scoreterm ) )
		self.next.grid( row = rownum, column = 2, sticky = W )
		rownum += 1

		colortoggle = Checkbutton( self, background = "#ddd",
		                           highlightthickness = 0,
		                           activebackground = "#ddd",
		                           text = "colors",
		                           variable = self.colorson )
		colortoggle.grid( row = rownum, column = 1 )

		if len(hits.hitlist) <= self.maxdisplay:
			self.next['state'] = DISABLED

	def browse( self, dir, app, hits, scoreterm ):
		if dir == 0: # higher in list
			self.hitstart -= self.maxdisplay
			if self.hitstart < 0: self.hitstart = 0
		else: # down in list
			self.hitstart += self.maxdisplay
			if self.hitstart >= len(hits.hitlist) - self.maxdisplay:
				self.hitstart = len(hits.hitlist) - self.maxdisplay

		self.prev['state'] = NORMAL
		self.next['state'] = NORMAL
		if self.hitstart >= len(hits.hitlist) - self.maxdisplay:
			self.next['state'] = DISABLED
		if self.hitstart == 0:
			self.prev['state'] = DISABLED

		self.rangelab['text'] = "showing hits %i to %i" % \
		                        ( self.hitstart, self.hitstart + self.maxdisplay )
		self.hitend = self.hitstart + self.maxdisplay
		self.hitarrays.set_arrays( app, hits, scoreterm, self.hitstart, self.maxdisplay,
		                           self.colorson.get() )

class HitArrays:

	def __init__( self, window, app, hits, scoreterm, hitstart, maxdisplay,
	              negative, positive, rownum ):

		self.hitlbls = []
		self.hitarys = []
		self.hittxts = []
		self.negative = negative
		self.positive = positive
		self.startrow = rownum

		i = 0
		while i < maxdisplay:
			rownum += 1

			# these objects must only be created once
			self.hitlbls.append( Label( window, font = smallfixed,
			                            background = "#ddd" ) )
			self.hitlbls[i].grid( row = rownum, column = 0, sticky = W )
			self.hitarys.append( SiteArray( window, app ) )
			self.hitarys[i].grid( row = rownum, column = 1 )
			self.hittxts.append( Label( window, font = smallfixed,
			                            background = "#ddd" ) )
			self.hittxts[i].grid( row = rownum, column = 2, sticky = W )
			i += 1

		self.endrow = rownum
		self.set_arrays( app, hits, scoreterm, hitstart, maxdisplay )

	def set_arrays( self, app, hits, scoreterm, hitstart, maxdisplay,
	                colorson = True ):

		hit_i = hitstart
		self.colors = []
		while hit_i < hitstart + maxdisplay:

			hit = hits.hitlist[hit_i]
			hitseq = hit.site
			if hitseq == "placeholder": continue # skip the dummy hit if encountered
			dir = hit.tag
			header = hit.gene
			maxbite = 15
			genepos = '  %s%7i%5s   ' % \
			          ( stringf_left(header,maxbite), hit.position + 1, dir )

			# rescore, skipping optimizations and storing energies by position
			# just translate rvs into fwd in this case
			if dir == "rvs":
				hitseq = rvs_comp_str( hitseq )

			score, scores = site_score_complete( hitseq, app.positions,
			                                     app.positions.pocketlength, scoreterm )

			# lowercase substitutions
			hitseq = hitseq.upper()
			hitseq = lowercase_substitutions( hitseq, app.wildtype )

			# make scorelist, color array
			colors = []
			if colorson == True:
				difs = avg_difs_to_single( app.scores_wt, scores, len(app.wildtype),
			  	                         app.positions.pocketlength )
				colors = make_colors( difs, hitseq, app.annotations, app.positions, self.positive,
			                           self.negative, app.colormode.get() )
			cell_i = hit_i - hitstart
			self.hitlbls[cell_i]['text'] = ( "%6.1f  " % score )
			self.hitarys[cell_i].set( hitseq, app.wildtype )
			self.hitarys[cell_i].set_colors( colors )
			self.hittxts[cell_i]['text'] = ( "%s" % genepos )
			hit_i += 1

class WTArray( Frame ):

	def __init__( self, master, wt, seqmap, pocketlength ):
		Frame.__init__( self, master, background = "#ddd" )

		# wildtype sequence with numerical position index
		indices = sorted( seqmap.keys(), as_int )
		# get max string length for layout purposes
		maxilen = 0
		for index in indices:
			ilen = len(index)
			if ilen > maxilen: maxilen = ilen
		self.cells = []
		ind_cells = []
		for i in range( len(wt) ):
			posindex = ""
			if i >= len(indices): # extrapolate terminal indices (if pocketlength > 1)
				posindex = str(int(indices[-1])+(i-len(indices))+1)
			else: posindex = indices[i]
			# justify indices for proper layout
			while ( len(posindex) < maxilen ):
				if posindex[0] == '-' or posindex[0] == ' ': posindex = posindex + ' '
				else: posindex = ' ' + posindex
			for j in range(len(posindex)):
				ind_cells.append( Cell( self, posindex[j] ) )
				ind_cells[-1].configure( background = "#ddd", font = smallfixed )
				ind_cells[-1].grid( row = j, column = i )

			self.cells.append( Cell( self ) )
			self.cells[i].grid( row = j+1, column = i )

		for i in range( len(self.cells) ):
			self.cells[i]['text'] = string.upper(wt[i])
			self.cells[i].configure( background = "#fff" )

class KeyArray( Frame ):

	def __init__( self, master, length, colors, negative, positive ):
		Frame.__init__( self, master )

		self.cells = []
		for i in range( length ):
			self.cells.append( Cell( self ) )
			self.cells[i].grid( row = 0, column = i+1 )

		zero = int( length*-1*float(negative)/(float(positive)-float(negative)) )

		for i in range( len(self.cells) ):
			symbol = " "
			if i == zero: symbol = "0"
			elif i == 0: symbol = "-"
			elif i == length - 1: symbol = "+"
			self.cells[i]['text'] = symbol
			self.cells[i].configure( background = colors[i] )

class SiteArray( Frame ):

	def __init__( self, master, app ):
		Frame.__init__( self, master )

		wt = app.wildtype
		self.cells = []
		for i in range( len(wt) ):
			self.cells.append( SiteCell( self, app, i ) )
			self.cells[i].grid( row = 0, column = i+1 )

	def set( self, hit, wt ):
		for i in range( len(self.cells) ):
			text = hit[i]
			self.cells[i]['text'] = text
			fontcol = "#000"
			#if hit[i].lower() != wt[i] and hit[i].upper() != wt[i]: fontcol = "#ccc" # ja temp
			self.cells[i].configure( foreground = fontcol )
#			if i < len(self.cells)-1:
#				self.cells[i].bind( "<Double-Button-1>", self.cells[i].load_pocket )

	def set_colors( self, colors ):
		if colors == []:
			for cell in self.cells: colors.append( "#3b3" )
		for i in range( len(self.cells) ):
			self.cells[i].configure( background = colors[i] )

class Cell( Label ):

	def __init__( self, master, text = "" ):
#		Label.__init__( self, master, text=text, font=bigfixed, width=10, height=10 )
		Label.__init__( self, master, text=text, font=bigfixed )

class SiteCell( Cell ):

	def __init__( self, master, app, index ):
		Cell.__init__( self, master )
		self.index = index
		self.app = app
		self.pocketseq = None

	# because it's a callback, no arguments are allowed (lame)
#	def load_pocket_pdb( self, event ):

class Titlemenu:

	def __init__( self, master, app ):

		titlemenu = Menu( master )
		filemenu = Menu( titlemenu, tearoff = 0 )

#		filemenu.add_command( label = "Load Rosetta data...",
#		                  command = lambda: app.loadfile( master ) )

#		filemenu.add_command( label = "Read pdbs from dir...",
#		                  command = lambda: app.loadpath( master ) )

#		filemenu.add_command( label = "Read best pdbs from dir...",
#		                  command = lambda: app.loadpath( master, 'best' ) )

		filemenu.add_command( label = "Load all...",
		                     command = lambda: Combo_match_dialog( master, app ) )

		filemenu.add_command( label = "Load sequence(s) to search...",
		                     command = lambda: app.loadseq( master ) )

		filemenu.add_command( label = "Load wildtype site sequence...",
#		                     command = lambda: File_match_dialog( master, 'wt', app ) )
		                     command = lambda: app.loadwt( master, app ) )

		filemenu.add_command( label = "Load PSSM file...",
		                      command = lambda: app.loadpssm( master ) )

		filemenu.add_command( label = "Load auxiliary matrix...",
		                      command = lambda: app.loadmatrix( master ) )

		filemenu.add_command( label = "Load lowercase annotations...",
		                      command = lambda: app.loadannotations( master ) )

		titlemenu.add_cascade( label = "Input", menu = filemenu )
		Seqmenu = Menu( titlemenu, tearoff = 0 )

#		Seqmenu.add_command( label = "Load auxiliary matrix...", command = lambda:
#		                     File_match_dialog( master, 'matrix', app ) )

		Seqmenu.add_command( label = "Search settings...",
		                     command = lambda: app.create_settings_dialog( master ))

		titlemenu.add_cascade( label = "Settings", menu = Seqmenu )

		titlemenu.add_command( label = "About",
		                       command = lambda: About( master ) )
		master.config( menu = titlemenu )

class Combo_match_dialog( Toplevel ):

	def __init__( self, master, app ):

		Toplevel.__init__( self, master, background = "#ddd" )
		lbls = {}
		self.mnus = {}
		row=0
#		for type in ['score','wt','matrix','seq']:
		for type in ['score','wt','seq','annotations']:
			lbls[ type ] = Label( self, text = string.capitalize(type),
			                      background = "#ddd" )
			lbls[ type ].grid( row = row, column = 0, sticky = W )
			self.mnus[ type ] = File_match_menu( self, type, app, 24 )
			self.mnus[ type ].grid( row = row, column = 1 )
			row+=1

		loadbut = Button( self, text = 'Load them',
                      highlightthickness = 0, command = self.load_all )
		loadbut.grid( row = row, columnspan = 2 )
		self.app = app

	def load_all( self ):
		for type, mnu in self.mnus.items():
			self.app.get_from_file( type, mnu.selection )
		self.destroy()

class File_match_dialog( Toplevel ):

	def __init__( self, master, match, app ):
		Toplevel.__init__( self, master, background = "#ddd" )

		self.match = match
		self.app = app
		lbl = Label( self, text = string.capitalize(match), background = "#ddd" )
		self.list = Listbox( self, font = "Courier -12", background = "#fff",
		                     highlightthickness = 0, width = 30 )
		listsb = Scrollbar( self, orient = VERTICAL )
		self.list.configure( yscrollcommand = listsb.set )
		self.list.bind( "<Double-Button-1>", self.load_file )
		listsb.configure( command = self.list.yview )
		lbl.grid( row = 0, column = 0 )
		self.list.grid( row = 1, column = 0, sticky = N+S )
		listsb.grid( row = 1, column = 1, sticky = N+S )

		fdbut = Button( self, text = '(It\'s not here!)',
		                highlightthickness = 0,
		                command = lambda: self.bigfd( app, match ) )

		fdbut.grid( row = 2, columnspan = 2, pady = 5 )

		files = matchlist( match )
		for file in files: self.list.insert( END, file )

	def bigfd( self, app, match ):
		fd = FileDialog.LoadFileDialog( self, title = "Load file..." )
		filename = fd.go()
		if filename != None:
			app.get_from_file( match, filename )
		self.destroy()

	def load_file( self, event ):
		self.app.get_from_file( self.match, self.frm.list.get( ACTIVE ) )
		self.destroy()

class File_match_menu( Menubutton ):

	def __init__( self, master, match, app, width_in, execute = False ):
		Menubutton.__init__( self, master, relief = RAISED, width = width_in )

		menu = Menu( self, tearoff = 0 )
		self.selection = ""
		files = matchlist( match )
		self.hard_addresses = {}
		index = 0
		for file in files:
			menu.add_command( label = file, command = lambda self = self, file = file:
			                  self.select( file, execute, app, match ) )

		if len(files) != 0:
			self['text'] = files[0]
			self.selection = files[0]
		else:
			self['text'] = 'no file'
		self['menu'] = menu

	def select( self, file, execute, app, match ):
		self.selection = file
		self['text'] = file
		if execute == True: app.get_from_file( match, file )

class About( Toplevel ):

	def __init__( self, master ):
		Toplevel.__init__( self, master )

		msg = "Written by Justin Ashworth\nBaker Lab, UW, 2006\nashwortj@u.washington.edu. Based on simple search algorithms written and developed by JA and Umut Ulge, Baker Lab, UW 2006.\n\nNote: Python is slow! This version is meant to offer interactive searching of relatively small amounts of sequence. For chromosome or genome-wide searches, ask JA about a much faster c++ version of the basic search algorithm.\nOr reimplement this in JAVA. ;)"
		aboutmsg = Message( self, width = 200, text = msg )
		closebutton = Button( self, text = "close", highlightthickness = 0,
		                      command = self.destroy )
		aboutmsg.grid()
		closebutton.grid( row = 1 )

class RClickDialog( Toplevel ):

	# a right-click bind will send this an 'event'
	def __init__( self, master, event ):
		Toplevel.__init__( self, master, background = "#ddd" )

		offset = master.geometry().split("+")
		offsetx = offset[1]
		offsety = offset[2]

		self.geometry( "%sx%s+%s+%s" % ( "60", "30",
		               str( event.x + int( offsetx ) ),
		               str( event.y + int( offsety ) ) ) )

		mbut = Menubutton( self, text = "Do...", relief = RAISED )

		mbut.m = Menu( mbut, tearoff = 0 )
		mbut.m.add_command( label = "Pymol", command = lambda:
		                    self.dial_pymol( self, master.tbox1.get( ACTIVE ) ) )
		mbut['menu'] = mbut.m
		mbut.grid()

	def dial_pymol( self, master, row ):

		keys = row.split()
		file = master.bestfiles[ master.pdbroot + "_" + keys[0] + "_" + keys[1] ]
		print file
		os.system( "pymol %s" % file )
		self.destroy()

### END OF GUI CLASES ###

### stand-alone functions ###

# pocket sorting functions, for the sole purpose of displaying the data
def sort_pockets( a, b, key, mode ):

	if mode == 'position' or mode == 'unique':
		# position indices are stored as strings, but sorted as ints
		if int(a.position) > int(b.position): return 1
		elif int(a.position) < int(b.position): return -1

	if mode == 'unique':
		if a.seq > b.seq: return 1
		elif a.seq < b.seq: return -1

	if a.energies[key] > b.energies[key]: return 1
	elif a.energies[key] < b.energies[key]: return -1
	return 0

def get_energies_from_pdbs( source, positions, bestfiles, term, mode = '' ):

	pdbs = {}
	print source

	if mode == 'best':
		for filename in os.listdir( source ):
			prefix = pdb_prefix( filename )
			if len(prefix) > 3:

				if positions.pdbroot == "": positions.pdbroot = prefix[0]
				prefix = prefix[0] + "_" + prefix[1] + "_" + prefix[2]
				if not pdbs.has_key( prefix ):
					pdbs[ prefix ] = []
				pdbs[ prefix ].append( filename )

		index = -1
		for prefix, pdblist in sorted( pdbs.items() ):
			print prefix, ":", len(pdblist), "model(s)"
			bestfile, bkt, ddg, hbb = get_pdb_energies_list( source, pdblist, term )
			bestfiles[ prefix ] = bestfile
			ids = prefix.split('_')
			position = ids[1]
			seq = ids[2]

			if positions.pocketlength == 0: positions.pocketlength = len(seq)
			if not positions.seqmap.has_key( position ):
				index += 1
				positions.seqmap[ position ] = index
				positions.list.append( {} )

			positions.list[index][seq] = Pocket( position, seq, bkt, ddg, hbb, bestfile )

	else: # store all pdbs
		index = -1
		for filename in sorted( glob.glob('*_*pdb') ):
			print filename
			prefix_list = pdb_prefix( filename )
			prefix = prefix_list[0]
			i = 1
			while i < len(prefix_list):
				prefix += '_'+prefix_list[i]
				if i == 2: break
				i += 1
			position = prefix_list[1]
			seq = prefix_list[2]
			E = get_pdb_energies( filename )
			bestfiles[ prefix ] = filename
			if not positions.seqmap.has_key( position ):
				index += 1
				positions.seqmap[ position ] = index
				positions.list.append( {} )
			positions.list[index][filename] = \
				Pocket( position, seq, E['bk_tot'], E['DDG'], E['base-spec'], filename )

def get_energies_from_file( source, positions, filenames ):

	file = open( source, "r" )
	index = -1
	for line in file:
		fields = line.split()
		position = fields[0]
		seq = fields[1]
		bkt = float(fields[2])
		ddg = float(fields[3])
		hbb = float(fields[4])
		filename = fields[5]

		if positions.pdbroot == "": positions.pdbroot = filename.split('_')[0]
		prefix = positions.pdbroot + "_" + position + "_" + seq
		filenames[ prefix ] = filename

		if positions.pocketlength == 0: positions.pocketlength = len(seq)

		if not positions.seqmap.has_key( position ):
			index += 1
			positions.seqmap[ position ] = index
			positions.list.append( {} )

		positions.list[ index ][ seq ] = Pocket( position, seq, bkt, ddg, hbb, filename )
	file.close()

def get_energies_from_pssm( source, positions ):

	file = open( source, 'r' )
	index = -1
	for line in file:
		# default base column keys (should be a,c,g,t)
		basekeys = bases
		if line[0] == '#pos': # read column header for keys, if it exists
			basekeys = line[1:-1] # file specifies column keys
			for base in basekeys: base = string.lower(base) # ensure lower-case
		elif line[0] == '#': continue # comments
		line = line.split()
		position = line[0]
		# fill positions
		if not positions.seqmap.has_key(position):
			index += 1
			positions.seqmap[position] = index
			positions.list.append( {} )
		# fill scores for each column(base)
		for i in range( len(basekeys) ):
			key = basekeys[i]
			score = float( line[i+1] )
			if options.invert:
				score = -1.0 * score
			positions.list[index][key] = Pocket( position, key, score, 0., 0., 'pssm' )
	file.close()
	# pssm assumed to work on single-basepair level
	positions.pocketlength = 1
	# labelling
	positions.pdbroot = source.split('/')[-1].split('.')[0]

# not all sorting terms correspond to scoring terms
def scorekey( key ):

	if key == '' or key == 'seq': return 'score'
	else: return key

def matchlist( match ):
	files = []
	matches = [ re.compile( match ) ]
	extras = []

	if match == 'score': extras = ['score','energ','pssm']
	if match == 'wt': extras = ['wildtype','native']
	elif match == 'matrix': extras = ['cen','middle','bonus','pssm']
	elif match == 'seq': extras = ['dna','fasta','gene','.ens']
	elif match == 'annotations': extras = ['ann']
	for extra in extras:
		matches.append( re.compile( extra ) )

	for file in sorted( os.listdir('.') ):
		if re_pdb.search( file ): continue
		if re_hit.search( file ): continue
		for match in matches:
			if match.search( file ):
				files.append( file )
				break
	return files

# these ones are optimized for speed (scoring priority and early rejection)
def site_score_fwd( site, positions, pocketlength, scoremap, fail, scoreterm ):

	score = 0.
	for i in range( len(site) - pocketlength + 1 ):
		position = scoremap.best_cases[i][0]
		seq = site[ position : position + pocketlength ]
#		print i, position, pocketlength, seq
		score += positions.list[ position ][ seq ].scores[ scoreterm ]
		if score + scoremap.best_cases[i][1] > fail:
#			print "failed site at index %i" % i
			return "fail"
#	print score, fail
	return score

def site_score_rvs( site, positions, pocketlength, scoremap, fail, scoreterm ):

	score = 0.
	for i in range( len(site) - pocketlength + 1 ):
		pos_rvs = len(site) - scoremap.best_cases[i][0] - pocketlength
		seq_rvs = ""
		j = pocketlength - 1
		while j >= 0:
			seq_rvs += comp[ site[pos_rvs+j] ]
			j -= 1
		score += \
		 positions.list[scoremap.best_cases[i][0]][seq_rvs].scores[ scoreterm ]
#		print "rvs", i, pos_rvs, seq_rvs, score
		if score + scoremap.best_cases[i][1] > fail:
#			print "failed rvs site at index %i" % i
			return "fail"
#	print score, fail
	return score

# more comfortable ones (for non-repetitive calling and full score info)
# also gives you a list
def site_score_complete( site, positions, pocketlength, scoreterm ):

	score = 0.
	scores = []
	for i in range( len(site) - pocketlength + 1 ):
		seq = site[ i : i + pocketlength ]
		score_i = positions.list[ i ][ seq ].scores[ scoreterm ]
		scores.append( score_i )
		score += score_i
	return score, scores

# No need for a reverse version of site_score_complete()!
# Just invert your sequence for reverse hits!

def load_simple_matrix( source ):

	matrix = {}
	print "opening simple weight matrix file %s\n" % source
	file = open( source, 'r' )
	for line in file:
		fields = line.split(':')
		if fields[0][0] == '#': continue
		position = fields[0].strip()
		if not matrix.has_key( position ): matrix[ position ] = []
		weights = fields[1].split('|')
		matrix[ position ] = tuple_weights( weights )
	return matrix

def tuple_weights( weights ):
	new = []
	for weight in weights:
		weight = weight.split(',')
		weight = ( weight[0].strip(), weight[1].strip() )
		new.append( weight )
	return new

def load_pocket_matrix( source ):

	matrix = []
	print "opening pocket weight matrix file %s\n" % source
	file = open( source, "r" )
	for line in file:
		fields = line.split()
		position = fields[0]
		seq = fields[1]
		energy = float(fields[2])
		matrix.append( Pocket( position, seq, energy, 0, 0, "" ) )
	return matrix

def avg_difs_to_single( wt_scores, scores, seqlength, pocketlength ):

	assert( seqlength == len(wt_scores) + pocketlength - 1 )
	assert( seqlength == len(scores) + pocketlength - 1 )

	# init lists
	difs = []
	for i in range( seqlength ):
		difs.append( [] )
	# accumulate scores
	for i in range( len(wt_scores) ):
		for j in range( pocketlength ):
			difs[i+j].append( scores[i] - wt_scores[i] )
#			print i, scores[i], wt_scores[i], difs[i]
	# average
	for i in range( len(difs) ):
		tot = 0
		for val in difs[i]:
			tot += val
		difs[i] = float(tot/float(len(difs[i])))
	return difs

def write_raw( positions, name ):
	datafile = open( name, "w" )
	for position in positions.list:
		for pocket in sorted( position.values() ):
			datafile.write( str(pocket) + '\n' )
	datafile.close()
	print 'score file \'%s\' written.\n' % name

def center_in_master( m_geom, width, height ):

	m = geometries( m_geom )
	x = m['x'] + m['w']/2 - width/2
	y = m['y'] + m['h']/2 - height/2

	s_geom = "%sx%s+%s+%s" % ( str(width), str(height), str(x), str(y) )
	return s_geom

# hash from Tkinter geometry string
def geometries( geometry ):

	geoms = {}
	terms = geometry.split("+")
	wh = terms[0].split("x")
	geoms['w'] = int( wh[0] )
	geoms['h'] = int( wh[1] )
	geoms['x'] = int( terms[1] )
	geoms['y'] = int( terms[2] )

	return geoms

def lookup_cutoffs( cutoffs, scoreterm ):

	negative = cutoffs[ scoreterm ][0].get()
	positive = cutoffs[ scoreterm ][1].get()
	return negative, positive

### init (main) ###

# script mode (no GUI)
if options.script:
	# old version deleted: to be (re)written
	pass

else:
	root = Tk()
	bigfixed = tkFont.Font( family = "Courier", size = "16", weight = "bold" )
	bigfixedlight = tkFont.Font( family = "Courier", size = "24" )
	smallfixed = tkFont.Font( family = "Courier", size = "12" )
#	bigfixed = tkFont.Font( family = "Courier", weight = "bold" )
#	bigfixedlight = tkFont.Font( family = "Courier" )
#	smallfixed = tkFont.Font( family = "Courier" )
	root.title( "PySSM" )
	root.wm_resizable(0,0)
	app = App( root )
	root.mainloop()
