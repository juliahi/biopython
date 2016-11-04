from galaxy import eggs

import pkg_resources
pkg_resources.require( "bx-python" )

import logging, os, sys, time, sets, tempfile, shutil
import data
from galaxy import util
from galaxy.datatypes.sniff import *
from cgi import escape
import urllib
from bx.intervals.io import *
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.tabular import Tabular

class Gaf( Tabular ):
    """Tab delimited data in GAF format"""
    file_ext = "gaf"

    #MetadataElement( name="columns", default=17, desc="Number of columns", readonly=True )
    #MetadataElement( name="comment_lines", default=0, desc="Number of comment lines", readonly=False, optional=True, no_value=0 )

    #MetadataElement( name="version", default="2.0", desc="GAF version", readonly=True )
    
    def __init__(self, **kwd):
        Tabular.__init__(self, **kwd)
        #self.do_something_else()


    def init_meta( self, dataset, copy_from=None ):
        Tabular.init_meta( self, dataset, copy_from=copy_from )
        self.column_names = ['DB',  'DB Object ID', 'DB Object Symbol', 'Qualifier',  'GO ID',
                             'DB:Reference (|DB:Reference)', 'Evidence Code', 'With (or) From', 
                             'Aspect', 'DB Object', 'DB Object Synonym (|Synonym)', 'DB Object Type', 
                             'Taxon(|taxon)', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID', 
                             ]


    def set_meta(  self, dataset, **kwd ):
        comment_lines = 0
        columns = 0
        headers = get_headers( dataset.file_name, '\t' )   
        for hdr in headers:
            if len( hdr ) > 1 and hdr[0]:
                if hdr[0].startswith( '!' ):
                    comment_lines += 1
                    #if hdr[0].startswith('!gaf-version:'): dataset.metadata.version = hdr[0][len('!gaf-version:'):].strip()
                else:
                    columns = max(columns, len(hdr))
                        
        dataset.metadata.comment_lines = comment_lines
        dataset.metadata.columns = columns
        dataset.metadata.column_names = self.column_names[:columns]
        
        
    def sniff( self, filename ):
        with open( filename ) as handle:
            for line in handle:
                if line.startswith('!gaf-version:'):
                    headers = get_headers( filename, '\t' )[:100]
                    for hdr in headers:
                        if len( hdr ) > 1 and hdr[0] and not hdr[0].startswith( '!' ):
                            if not 15 <= len(hdr) <= 17:
                                return False
                    return True
                elif not line.startswith('!'):
                    return False
                            
        return False
