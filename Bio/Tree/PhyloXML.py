# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes corresponding to phyloXML elements.
"""
__docformat__ = "epytext en"

import re
import warnings

from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import BaseTree
from BaseTree import trim_str


class PhyloXMLWarning(Warning):
    """Warning for non-compliance with the phyloXML specification."""
    pass


def check_str(text, testfunc):
    """Check a string using testfunc, and warn if there's no match."""
    if text is not None and not testfunc(text):
        warnings.warn("String %s doesn't match the given regexp" % text,
                      PhyloXMLWarning, stacklevel=2)


# Core elements

class PhyloElement(BaseTree.TreeElement):
    """Base class for all PhyloXML objects."""
    def __init__(self, **kwargs):
        """Set all keyword arguments as instance attributes.
        """
        self.__dict__.update(kwargs)

    def __str__(self):
        """Show the class name and an identifying attribute."""
        s = self.__class__.__name__
        if hasattr(self, 'name') and self.name:
            return '%s %s' % (s, trim_str(self.name))
        if hasattr(self, 'value') and self.value:
            return '%s %s' % (s, trim_str(self.value))
        if hasattr(self, 'id') and self.id:
            return '%s %s' % (s, self.id)
        return s


class Phyloxml(PhyloElement):
    """Root node of the PhyloXML document.

    Contains an arbitrary number of Phylogeny elements, possibly followed by
    elements from other namespaces.

    @param attributes: (XML namespace definitions)
    @param phylogenies: list of phylogenetic trees
    @param other: list of arbitrary non-phyloXML elements, if any
    """
    def __init__(self, attributes, phylogenies=None, other=None):
        self.attributes = attributes
        self.phylogenies = phylogenies or []
        self.other = other or []

    def __getitem__(self, index):
        """Get a phylogeny by index or name."""
        if isinstance(index, int) or isinstance(index, slice):
            return self.phylogenies[index]
        if not isinstance(index, basestring):
            raise KeyError, "can't use %s as an index" % type(index)
        for tree in self.phylogenies:
            if tree.name == index:
                return tree
        else:
            raise KeyError, "no phylogeny found with name " + repr(index)

    def __iter__(self):
        """Iterate through the phylogenetic trees in this object."""
        return iter(self.phylogenies)

    def __len__(self):
        """Number of phylogenetic trees in this object."""
        return len(self.phylogenies)


class Other(PhyloElement):
    """Container for non-phyloXML elements in the tree.

    Usually, an Other object will have either a 'value' or a non-empty list
    of 'children', but not both. This is not enforced here, though.

    @param tag: local tag for the XML node
    @param namespace: XML namespace for the node -- should not be the default
        phyloXML namespace.
    @param attributes: string attributes on the XML node
    @param value: text contained directly within this XML node
    @param children: list of child nodes, if any (also Other instances)
    """
    def __init__(self, tag, namespace=None, attributes=None, value=None,
            children=None):
        self.tag = tag
        self.namespace = namespace
        self.attributes = attributes
        self.value = value
        self.children = children or []

    def __iter__(self):
        """Iterate through the children of this object (if any)."""
        return iter(self.children)


class Phylogeny(PhyloElement, BaseTree.Tree):
    """A phylogenetic tree.

    @param rooted: True if this tree is rooted
    @param rerootable: True if this tree is rerootable
    @param branch_length_unit: unit for branch_length values on clades
    @type type: str

    @param name: string identifier for this tree, not required to be unique
    @param id: unique identifier for this tree (type Id)
    @param description: plain-text description
    @param date: date for the root node of this tree (type Date)
    @param confidences: list of Confidence objects for this tree
    @param clade: the root node/clade of this tree
    @param clade_relations: list of CladeRelation objects
    @param sequence_relations: list of SequenceRelation objects
    @param properties: list of Property objects
    @param other: list of non-phyloXML elements (type Other)
    """
    def __init__(self, rooted,
            rerootable=None, branch_length_unit=None, type=None,
            # Child nodes
            name=None, id=None, description=None, date=None, clade=None,
            # Collections
            confidences=None, clade_relations=None, sequence_relations=None,
            properties=None, other=None,
            ):
        assert isinstance(rooted, bool)
        PhyloElement.__init__(self, rerootable=rerootable,
                branch_length_unit=branch_length_unit, type=type,
                rooted=rooted, name=name, id=id, description=description,
                date=date, clade=clade,
                confidences=confidences or [],
                clade_relations=clade_relations or [],
                sequence_relations=sequence_relations or [],
                properties=properties or [],
                other=other or [],
                )

    def to_phyloxml(self, **kwargs):
        """Create a new PhyloXML object containing just this phylogeny."""
        return Phyloxml(kwargs, phylogenies=[self])

    @property
    def confidence(self):
        """Equivalent to self.confidences[0] if there is only 1 value.

        See also: Clade.confidence, Clade.taxonomy
        """
        if len(self.confidences) == 0:
            raise RuntimeError("Phylogeny().confidences is empty")
        if len(self.confidences) > 1:
            raise RuntimeError("more than 1 confidence value available; "
                               "use Phylogeny().confidences")
        return self.confidences[0]


class Clade(PhyloElement, BaseTree.Node, BaseTree.Tree):
    """Describes a branch of the current phylogenetic tree.

    Used recursively, describes the topology of a phylogenetic tree.

    Both 'color' and 'width' elements apply for the whole clade unless
    overwritten in-sub clades.

    @param branch_length: parent branch length of this clade
    @param id_source: link other elements to a clade (on the xml-level)

    @param name: short string label for this clade
    @param confidences: list of Confidence objects, used to indicate the
        support for a clade/parent branch.
    @param width: branch width for this clade (including parent branch)
    @param color: color used for graphical display of this clade
    @param node_id: unique identifier for the root node of this clade
    @param taxonomies: list of Taxonomy objects
    @param sequences: list of Sequence objects
    @param events: describe such events as gene-duplications at the root
        node/parent branch of this clade
    @param binary_characters: a BinaryCharacters object
    @param distributions: list of Distribution objects
    @param date: a date for the root node of this clade (type Date)
    @param references: list of Reference objects
    @param properties: list of Property objects
    @param clades: list of sub-clades (type Clade)
    @param other: list of non-phyloXML objects

    @param left_idx: precomputed values for certain tree operations
    @param right_idx: precomputed values for certain tree operations
    """
    def __init__(self, parent=None,
            # Attributes
            branch_length=None, id_source=None,
            # Child nodes
            name=None, width=None, color=None, node_id=None, events=None,
            binary_characters=None, date=None,
            # Collections
            confidences=None, taxonomies=None, sequences=None,
            distributions=None, references=None, properties=None, clades=None,
            other=None,
            # BaseTree.Node
            left_idx=None, right_idx=None,
            ):
        PhyloElement.__init__(self, parent=parent, id_source=id_source,
                name=name, branch_length=branch_length, width=width,
                color=color, node_id=node_id, events=events,
                binary_characters=binary_characters, date=date,
                left_idx=left_idx, right_idx=right_idx,
                confidences=confidences or [],
                taxonomies=taxonomies or [],
                sequences=sequences or [],
                distributions=distributions or [],
                references=references or [],
                properties=properties or [],
                clades=clades or [],
                other=other or [],
                )

    def to_phylogeny(self, **kwargs):
        """Create a new phylogeny containing just this clade."""
        # ENH: preserve some attributes of the parent phylogeny
        return Phylogeny(clade=self, **kwargs)

    # Mimic BaseTree.Node
    @property
    def label(self):
        return str(self)

    # Mimic BaseTree.Tree
    @property
    def rooted(self):
        if self.parent is not None:
            return self.parent.rooted

    # Shortcuts for list attributes that are usually only 1 item
    # XXX should these raise RuntimeError, AttributeError or IndexError?
    @property
    def confidence(self):
        if len(self.confidences) == 0:
            raise RuntimeError("Clade().confidences is empty")
        if len(self.confidences) > 1:
            raise RuntimeError("more than 1 confidence value available; "
                               "use Clade().confidences")
        return self.confidences[0]

    @property
    def taxonomy(self):
        if len(self.taxonomies) == 0:
            raise RuntimeError("Clade().taxonomies is empty")
        if len(self.taxonomies) > 1:
            raise RuntimeError("more than 1 taxonomy value available; "
                               "use Clade().taxonomies")
        return self.taxonomies[0]

    # Sequence-type behavior methods

    def __getitem__(self, index):
        """Get a sub-clade by index (integer or slice)."""
        if isinstance(index, int) or isinstance(index, slice):
            return self.clades[index]
        ref = self
        for idx in index:
            ref = ref.clades[idx]
        return ref

    def __iter__(self):
        """Iterate through the clades (sub-nodes) within this clade."""
        return iter(self.clades)

    def __len__(self):
        """Number of clades directy under this element."""
        return len(self.clades)


# PhyloXML-specific complex types

class Accession(PhyloElement):
    """Captures the local part in a sequence identifier.

    Example: In 'UniProtKB:P17304', the value of Accession is 'P17304'  and the
    'source' attribute is 'UniProtKB'.
    """
    def __init__(self, value, source):
        self.value = value
        self.source = source


class Annotation(PhyloElement):
    """The annotation of a molecular sequence.

    It is recommended to annotate by using the optional 'ref' attribute (some
    examples of acceptable values for the ref attribute: 'GO:0008270',
    'KEGG:Tetrachloroethene degradation', 'EC:1.1.1.1').

    @type ref: str
    @param source: plain-text source for this annotation
    @param evidence: describe evidence as free text (e.g. 'experimental')
    @type type: str

    @param desc: free text description
    @param confidence: state the type and value of support (type Confidence)
    @param properties: list of typed and referenced annotations from external
        resources
    @type uri: Uri
    """
    re_ref = re.compile(r'[a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+')

    def __init__(self, 
            # Attributes
            ref=None, source=None, evidence=None, type=None,
            # Child nodes
            desc=None, confidence=None, uri=None,
            # Collection
            properties=None):
        check_str(ref, self.re_ref.match)
        PhyloElement.__init__(self, ref=ref, source=source, evidence=evidence,
                type=type, desc=desc, confidence=confidence, uri=uri,
                properties=properties or [])


class BinaryCharacters(PhyloElement):
    """The names and/or counts of binary characters present, gained, and lost
    at the root of a clade. 
    """
    def __init__(self,
            # Attributes
            type=None, gained_count=None, lost_count=None, present_count=None,
            absent_count=None,
            # Child nodes (flattened into collections)
            gained=None, lost=None, present=None, absent=None):
        PhyloElement.__init__(self,
                type=type, gained_count=gained_count, lost_count=lost_count,
                present_count=present_count, absent_count=absent_count,
                gained=gained or [],
                lost=lost or [],
                present=present or [],
                absent=absent or [])


class BranchColor(PhyloElement):
    """Indicates the color of a clade when rendered graphically.

    The color applies to the whole clade unless overwritten by the color(s) of
    sub-clades.

    Color values should be unsigned bytes, or integers from 0 to 255.
    """
    def __init__(self, red, green, blue):
        assert isinstance(red, int)
        assert isinstance(green, int)
        assert isinstance(blue, int)
        self.red = red
        self.green = green
        self.blue = blue

    def to_rgb(self):
        """Return a 24-bit hexadecimal RGB representation of this color.

        The returned string is suitable for use in HTML/CSS.

        Example:

            >>> bc = BranchColor(12, 200, 100)
            >>> bc.to_rgb()
            '0cc864'
        """
        return hex(self.red * (16**4)
                + self.green * (16**2)
                + self.blue)[2:].zfill(6)


class CladeRelation(PhyloElement):
    """Expresses a typed relationship between two clades.

    For example, this could be used to describe multiple parents of a clade.

    @type id_ref_0: str
    @type id_ref_1: str
    @type distance: str
    @type type: str

    @type confidence: Confidence
    """
    def __init__(self, type, id_ref_0, id_ref_1,
            distance=None, confidence=None):
        PhyloElement.__init__(self, distance=distance, type=type,
                id_ref_0=id_ref_0, id_ref_1=id_ref_1, confidence=confidence)


class Confidence(PhyloElement):
    """A general purpose confidence element.

    For example, this can be used to express the bootstrap support value of a
    clade (in which case the 'type' attribute is 'bootstrap').

    @type value: float
    @type type: str
    """
    def __init__(self, value, type):
        self.value = value
        self.type = type


class Date(PhyloElement):
    """A date associated with a clade/node.

    Its value can be numerical by using the 'value' element and/or free text
    with the 'desc' element' (e.g. 'Silurian'). If a numerical value is used, it
    is recommended to employ the 'unit' attribute.

    @param unit: type of numerical value (e.g. 'mya' for 'million years ago')

    @type value: float
    @param desc: plain-text description of the date
    @param minimum: lower bound on the date value
    @param maximum: upper bound on the date value
    """
    def __init__(self, value=None, unit=None, desc=None, 
            minimum=None, maximum=None):
        PhyloElement.__init__(self, value=value, unit=unit, desc=desc, 
                minimum=minimum, maximum=maximum)

    def __str__(self):
        """Show the class name and the human-readable date."""
        s = self.__class__.__name__
        if self.unit and self.value is not None:
            return '%s %s %s' % (s, self.value, self.unit)
        if self.desc is not None:
            return '%s %s' % (s, self.desc)
        return s


class Distribution(PhyloElement):
    """Geographic distribution of the items of a clade (species, sequences).

    Intended for phylogeographic applications.

    The location can be described either by free text in the 'desc' element
    and/or by the coordinates of one or more 'Points' (similar to the 'Point'
    element in Google's KML format) or by 'Polygons'.
    """
    def __init__(self, desc=None, points=None, polygons=None):
        PhyloElement.__init__(self, desc=desc,
                points=points or [],
                polygons=polygons or [])


class DomainArchitecture(PhyloElement):
    """Domain architecture of a protein.

    @param length: total length of the protein sequence (type int)
    @param domains: list of ProteinDomain objects
    """
    def __init__(self, length=None, domains=None):
        # assert len(domains)
        PhyloElement.__init__(self, length=length, domains=domains)


class Events(PhyloElement):
    """Events at the root node of a clade (e.g. one gene duplication).

    All attributes are set to None by default, but this object can also be
    treated as a dictionary, in which case None values are treated as missing
    keys and deleting a key resets that attribute's value back to None.
    """
    ok_type = set(('transfer', 'fusion', 'speciation_or_duplication', 'other',
                    'mixed', 'unassigned'))

    def __init__(self, type=None, duplications=None, speciations=None,
            losses=None, confidence=None):
        check_str(type, self.ok_type.__contains__)
        PhyloElement.__init__(self, type=type, duplications=duplications,
                speciations=speciations, losses=losses, confidence=confidence)

    def iteritems(self):
        return ((k, v) for k, v in self.__dict__.iteritems() if v is not None)

    def iterkeys(self):
        return (k for k, v in self.__dict__.iteritems() if v is not None)

    def itervalues(self):
        return (v for v in self.__dict__.itervalues() if v is not None)

    def items(self):
        return list(self.iteritems())

    def keys(self):
        return list(self.iterkeys())

    def values(self):
        return list(self.itervalues())

    def __len__(self):
        return len(self.values())

    def __getitem__(self, key):
        if not hasattr(self, key):
            raise KeyError(key)
        val = getattr(self, key)
        if val is None:
            raise KeyError("%s has not been set in this object" % repr(key))
        return val

    def __setitem__(self, key, val):
        setattr(self, key, val)

    def __delitem__(self, key):
        setattr(self, key, None)

    def __iter__(self):
        return iter(self.iterkeys())

    def __contains__(self, key):
        return (hasattr(self, key) and getattr(self, key) is not None)


class Id(PhyloElement):
    """A general purpose identifier element.

    Allows to indicate the provider (or authority) of an identifier, e.g. NCBI,
    along with the value itself.
    """
    def __init__(self, value, provider=None):
        PhyloElement.__init__(self, provider=provider, value=value)


class MolSeq(PhyloElement):
    """Store a molecular sequence.

    @param value: the sequence, as a string
    @param is_aligned: True is mol_seq is aligned (usu. meaning gaps are
        introduced and all aligned seqs are the same length)
    """
    re_value = re.compile(r'[a-zA-Z\.\-\?\*_]+')

    def __init__(self, value, is_aligned=None):
        check_str(value, self.re_value.match)
        self.value = value
        self.is_aligned = is_aligned

    def __str__(self):
        return self.value


class Point(PhyloElement):
    """Geographic coordinates of a point, with an optional altitude.

    Used by element 'Distribution'.

    @param geodetic_datum: indicate the geodetic datum (also called 'map
        datum'). For example, Google's KML uses 'WGS84'. (required)
    @param lat: latitude
    @param long: longitude
    @param alt: altitude
    @param alt_unit: unit for the altitude (e.g. 'meter')
    """
    def __init__(self, geodetic_datum, lat, long, alt=None, alt_unit=None):
        PhyloElement.__init__(self, geodetic_datum=geodetic_datum,
                lat=lat, long=long, alt=alt, alt_unit=alt_unit)


class Polygon(PhyloElement):
    """A polygon defined by a list of 'Points' (used by element 'Distribution').

    @param points: list of 3 or more points representing vertices.
    """
    def __init__(self, points=None):
        self.points = points or []


class Property(PhyloElement):
    """A typed and referenced property from an external resources.

    Can be attached to 'Phylogeny', 'Clade', and 'Annotation' objects.

    @param ref: reference to an external resource, e.g. "NOAA:depth"

    @param unit: the unit of the property, e.g. "METRIC:m" (optional)

    @param datatype: indicates the type of a property and is limited to
        xsd-datatypes (e.g. 'xsd:string', 'xsd:boolean', 'xsd:integer',
        'xsd:decimal', 'xsd:float', 'xsd:double', 'xsd:date', 'xsd:anyURI').

    @param applies_to: indicates the item to which a property applies to (e.g.
        'node' for the parent node of a clade, 'parent_branch' for the parent
        branch of a clade, or just 'clade').

    @param id_ref: allows to attached a property specifically to one element
        (on the xml-level). (optional)

    @type value: str
    """
    re_ref = re.compile(r'[a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+')
    ok_applies_to = set(('phylogeny', 'clade', 'node', 'annotation',
                         'parent_branch', 'other'))
    ok_datatype = set(('xsd:string', 'xsd:boolean', 'xsd:decimal', 'xsd:float',
        'xsd:double', 'xsd:duration', 'xsd:dateTime', 'xsd:time', 'xsd:date',
        'xsd:gYearMonth', 'xsd:gYear', 'xsd:gMonthDay', 'xsd:gDay',
        'xsd:gMonth', 'xsd:hexBinary', 'xsd:base64Binary', 'xsd:anyURI',
        'xsd:normalizedString', 'xsd:token', 'xsd:integer',
        'xsd:nonPositiveInteger', 'xsd:negativeInteger', 'xsd:long', 'xsd:int',
        'xsd:short', 'xsd:byte', 'xsd:nonNegativeInteger', 'xsd:unsignedLong',
        'xsd:unsignedInt', 'xsd:unsignedShort', 'xsd:unsignedByte',
        'xsd:positiveInteger'))

    def __init__(self, value, ref, applies_to, datatype,
            unit=None, id_ref=None):
        check_str(ref, self.re_ref.match)
        check_str(applies_to, self.ok_applies_to.__contains__)
        check_str(datatype, self.ok_datatype.__contains__)
        check_str(unit, self.re_ref.match)
        PhyloElement.__init__(self, unit=unit, id_ref=id_ref, value=value,
                ref=ref, applies_to=applies_to, datatype=datatype)


class ProteinDomain(PhyloElement):
    """Represents an individual domain in a domain architecture.

    The locations use 0-based indexing, as most Python objects including
    SeqFeature do, rather than the usual biological convention starting at 1.
    This means the start and end attributes can be used directly as slice
    indexes on Seq objects.

    @param start: start of the domain on the sequence, using 0-based indexing
    @type start: non-negative integer
    @param end: end of the domain on the sequence
    @type end: non-negative integer
    @param confidence: can be used to store e.g. E-values. (type float)
    @param id: unique identifier/name
    """
    # TODO: confirm that 'start' counts from 1, not 0
    def __init__(self, value, start, end, confidence=None, id=None):
        PhyloElement.__init__(self, value=value, start=start, end=end,
                confidence=confidence, id=id)

    @classmethod
    def from_seqfeature(cls, feat):
        return ProteinDomain(feat.id,
                feat.location.nofuzzy_start,
                feat.location.nofuzzy_end,
                confidence=feat.qualifiers.get('confidence'))

    def to_seqfeature(self):
        feat = SeqFeature(location=FeatureLocation(self.start, self.end),
                          id=self.value)
        if hasattr(self, 'confidence'):
            feat.qualifiers['confidence'] = self.confidence
        return feat


class Reference(PhyloElement):
    """Literature reference for a clade.

    It is recommended to use the 'doi' attribute instead of the free text
    'desc' element whenever possible.
    """
    re_doi = re.compile(r'[a-zA-Z0-9_\.]+/[a-zA-Z0-9_\.]+')

    def __init__(self, doi=None, desc=None):
        check_str(doi, self.re_doi.match)
        self.doi = doi
        self.desc = desc


class Sequence(PhyloElement):
    """A molecular sequence (Protein, DNA, RNA) associated with a node.

    One intended use for 'id_ref' is to link a sequence to a taxonomy (via the
    taxonomy's 'id_source') in case of multiple sequences and taxonomies per
    node. 

    @param type: type of sequence ('dna', 'rna', or 'protein').
    @type id_ref: str
    @type id_source: str

    @param symbol: short  symbol of the sequence, e.g. 'ACTM' (max. 10 chars)
    @type accession: Accession
    @param name: full name of the sequence, e.g. 'muscle Actin'
    @param location: location of a sequence on a genome/chromosome.
    @param mol_seq: the actual sequence, as a string
    @type uri: Uri
    @param annotations: list of Annotation objects
    @param domain_architecture: protein domains on this sequence (type
        DomainArchitecture)
    @param other: list of non-phyloXML elements (type Other)
    """
    re_symbol = re.compile(r'\S{1,10}')
    ok_type = set(('rna', 'dna', 'protein'))

    def __init__(self, 
            # Attributes
            type=None, id_ref=None, id_source=None,
            # Child nodes
            symbol=None, accession=None, name=None, location=None,
            mol_seq=None, uri=None, domain_architecture=None,
            # Collections
            annotations=None, other=None,
            ):
        check_str(type, self.ok_type.__contains__)
        check_str(symbol, self.re_symbol.match)
        PhyloElement.__init__(self, type=type, id_ref=id_ref,
                id_source=id_source, symbol=symbol, accession=accession,
                name=name, location=location, mol_seq=mol_seq, uri=uri,
                domain_architecture=domain_architecture,
                annotations=annotations or [],
                other=other or [],
                )

    @classmethod
    def from_seqrecord(cls, record):
        kwargs = {
                'accession': Accession('', record.id),
                'symbol': record.name,
                'name': record.description,
                'mol_seq': str(record.seq),
                }
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            kwargs['type'] = 'dna'
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            kwargs['type'] = 'rna'
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            kwargs['type'] = 'protein'

        # Unpack record.annotations
        annot_attrib = {}
        annot_conf = None
        annot_prop = None
        annot_uri = None
        for key in ('ref', 'source', 'evidence', 'type'):
            if key in record.annotations:
                annot_attrib[key] = record.annotations[key]
        if 'confidence' in record.annotations:
            # NB: record.annotations['confidence'] = [value, type]
            annot_conf = Confidence(*record.annotations['confidence'])
        if 'properties' in record.annotations:
            # NB: record.annotations['properties'] = {...}
            annot_props = [Property(**prop)
                           for prop in record.annotations['properties']]
        if 'uri' in record.annotations:
            # NB: record.annotations['uri'] = {...}
            annot_uri = Uri(**record.annotations['uri'])
        kwargs['annotations'] = [Annotation(annot_attrib, {
            'desc': record.annotations.get('desc', None),
            'confidence': annot_conf,
            'properties': [annot_prop],
            'uri': annot_uri,
            })]

        # Unpack record.features
        if record.features:
            kwargs['domain_architecture'] = DomainArchitecture(
                    length=len(record.seq),
                    domains=[ProteinDomain.from_seqfeature(feat)
                             for feat in record.features])

        # Not handled:
        # attributes: id_ref, id_source
        # kwargs['location'] = None
        # kwargs['uri'] = None -- redundant here?
        return Sequence(**kwargs)

    def to_seqrecord(self):
        alphabets = {'dna': Alphabet.generic_dna,
                     'rna': Alphabet.generic_rna,
                     'protein': Alphabet.generic_protein}
        seqrec = SeqRecord(
                Seq(self.mol_seq,
                    alphabets.get(self.type, Alphabet.generic_alphabet)),
                id=str(self.accession),
                name=self.symbol,
                description=self.name,
                # dbxrefs=None,
                # features=None,
                )
        # TODO: repack seqrec.annotations
        return seqrec


class SequenceRelation(PhyloElement):
    """Express a typed relationship between two sequences.

    For example, this could be used to describe an orthology (in which case
    attribute 'type' is 'orthology'). 

    @param id_ref_0: first sequence reference identifier
    @param id_ref_1: second sequence reference identifier
    @param distance: distance between the two sequences (type float)
    @param type: describe the type of relationship

    @type confidence: Confidence
    """
    ok_type = set(('orthology', 'one_to_one_orthology', 'super_orthology',
        'paralogy', 'ultra_paralogy', 'xenology', 'unknown', 'other'))

    def __init__(self, type, id_ref_0, id_ref_1,
            distance=None, confidence=None):
        check_str(type, self.ok_type.__contains__)
        PhyloElement.__init__(self, distance=distance, type=type,
                id_ref_0=id_ref_0, id_ref_1=id_ref_1, confidence=confidence)


class Taxonomy(PhyloElement):
    """Describe taxonomic information for a clade.

    @param id_source: link other elements to a taxonomy (on the XML level)

    @param id: unique identifier of a taxon, e.g. Id('6500',
        provider='ncbi_taxonomy') for the California sea hare
    @param code: store UniProt/Swiss-Prot style organism codes, e.g. 'APLCA'
        for the California sea hare 'Aplysia californica' (restricted string)
    @param scientific_name: the standard scientific name for this organism,
        e.g. 'Aplysia californica' for the California sea hare
    @param authority: keep the authority, such as 'J. G. Cooper, 1863',
        associated with the 'scientific_name'
    @param common_names: list of common names for this organism
    @param synonyms: ???
    @param rank: taxonomic rank (restricted string)
    @type uri: Uri
    @param other: list of non-phyloXML elements (type Other)
    """
    re_code = re.compile(r'[a-zA-Z0-9_]{2,10}')
    ok_rank = set(('domain', 'kingdom', 'subkingdom', 'branch', 'infrakingdom',
        'superphylum', 'phylum', 'subphylum', 'infraphylum', 'microphylum',
        'superdivision', 'division', 'subdivision', 'infradivision',
        'superclass', 'class', 'subclass', 'infraclass', 'superlegion',
        'legion', 'sublegion', 'infralegion', 'supercohort', 'cohort',
        'subcohort', 'infracohort', 'superorder', 'order', 'suborder',
        'superfamily', 'family', 'subfamily', 'supertribe', 'tribe', 'subtribe',
        'infratribe', 'genus', 'subgenus', 'superspecies', 'species',
        'subspecies', 'variety', 'subvariety', 'form', 'subform', 'cultivar',
        'unknown', 'other'))

    def __init__(self, 
            # Attributes
            id_source=None,
            # Child nodes
            id=None, code=None, scientific_name=None, authority=None,
            rank=None, uri=None,
            # Collections
            common_names=None, synonyms=None, other=None,
            ):
        check_str(code, self.re_code.match)
        check_str(rank, self.ok_rank.__contains__)
        PhyloElement.__init__(self, id_source=id_source, id=id, code=code,
                scientific_name=scientific_name, authority=authority,
                rank=rank, uri=uri,
                common_names=common_names or [],
                synonyms=synonyms or [],
                other=other or [],
                )

    def __str__(self):
        """Show the class name and an identifying attribute."""
        s = self.__class__.__name__
        if self.code is not None:
            return '%s %s' % (s, self.code)
        if self.scientific_name is not None:
            return '%s %s' % (s, self.scientific_name)
        if self.rank is not None:
            return '%s %s' % (s, self.rank)
        if self.id is not None:
            return '%s %s' % (s, self.id)
        return s


class Uri(PhyloElement):
    """A uniform resource identifier.

    In general, this is expected to be an URL (for example, to link to an image
    on a website, in which case the 'type' attribute might be 'image' and 'desc'
    might be 'image of a California sea hare').
    """
    def __init__(self, value, desc=None, type=None):
        PhyloElement.__init__(self, value=value, desc=desc, type=type)
