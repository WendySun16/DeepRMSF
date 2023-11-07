# --- UCSF Chimera Copyright ---
# Copyright (c) 2000-2011 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  This notice must be embedded in or
# attached to all copies, including partial copies, of the
# software or any revisions or derivations thereof.
# --- UCSF Chimera Copyright ---
#
# $Id: copyright 34705 2011-10-19 23:37:43Z pett $

#TODO: nucleotides is a _lot_ of work
#TODO: then: clipping (check against #18620)
#TODO: then: 2D arrows
#TODO: then: volumes
#TODO: then: sequences: intrinsic; finish headers (settings); letter coloring(?)

import replyobj
pipes_string = " Pipes and Planks"
surface_export_warned = False
label_2d_rgba_warned = False
label_2d_size_warned = False
label_2d_style_warned = False
label_2d_font_name_warned = False
label_2d_multiline_warned = False

def export_cx(filename):
	from OpenSave import osOpen
	f = osOpen(filename, 'w')
	write_utility_functions(f)

	# 3D models
	write_prolog(f)
	by_id = {}
	metal_pbg_lookup = {}
	msms_models = []
	from chimera import openModels, viewer
	for m in openModels.list():
		if m.__class__.__name__ == 'Molecule':
			write_molecule(f, m)
			pbg = m.metalComplexGroup(create=False)
			if pbg:
				metal_pbg_lookup[pbg] = m
			by_id.setdefault(m.id, []).append(m)
		elif m.__class__.__name__ == 'MSMSModel':
			msms_models.append(m)
			by_id.setdefault(m.molecule, []).append(m)
		else:
			if not m.name.endswith(pipes_string):
				replyobj.warning("ChimeraX export does not yet support %s models" % m.__class__.__name__)
			continue
	global_groups = write_pseudobonds(f, metal_pbg_lookup)
	for msms_model in msms_models:
		write_msms_model(f, msms_model)
	write_selection(f)

	# 2D models
	from Ilabel import LabelsModel
	labels_2d_model = LabelsModel(create=False)
	if labels_2d_model:
		write_2d_labels(f, labels_2d_model)
	from Ilabel.ColorKey import getKeyModel
	key = getKeyModel(create=False)
	if key:
		write_color_key(f, key)

	# scene parameters
	cam = viewer.camera
	bg = viewer.background if viewer.backgroundMethod == viewer.Solid else None
	write_epilog(f, by_id, global_groups, cam.eyePos(0), cam.fieldOfView, viewer.windowSize, bg)

	# alignments
	from chimera.extension import manager
	from MultAlignViewer.MAViewer import MAViewer
	for mav in [inst for inst in manager.instances if isinstance(inst, MAViewer)]:
		write_alignment(f, mav)
	f.close()

def pickled(data):
	from SimpleSession.save import pickled as pickled2
	#for i in range(len(data)):
	#	pickled2(data[i])
	return "p" + pickled2(data)[2:-1] + ", encoding='latin1')"  # change the initial 'cPickle' to 'pickle'

def color_val(color):
	return None if color is None else color.rgba()

def simple_offset(off):
	return None if off is None else off.data()

def write_utility_functions(f):
	print>>f, color_converter_definition
	print>>f, restore_model_definition
	print>>f, make_structure_definition
	print>>f, make_molecular_surface_definition
	print>>f, restore_coordset_definition
	print>>f, restore_res_definition
	print>>f, restore_atom_definition
	print>>f, restore_bond_definition
	print>>f, make_pseudobonds_definition
	print>>f, restore_selection_definition
	print>>f, restore_camera_definition
	print>>f, restore_window_size_definition
	print>>f, make_alignment_definition
	print>>f, restore_seq_definition
	print>>f, restore_2d_labels_definition
	print>>f, restore_color_key_definition

def write_prolog(f):
	print>>f, "from chimerax.core.commands import run"
	print>>f, "run(session, 'close session')"
	print>>f, "import pickle, base64"
	print>>f, "global_atom_map = {}"
	print>>f, "structure_map = {}"
	print>>f, "residue_map = {}"
	print>>f, "label_data = []"
	print>>f, "font_mapping = {'Serif': 'Times', 'Sans Serif': 'Helvetica', 'Fixed': 'Courier'}"

def write_pseudobonds(f, metal_pbg_lookup):
	from chimera import PseudoBondMgr
	mgr = PseudoBondMgr.mgr()
	global_groups = []
	data = []
	from StructMeasure.DistMonitor import monitoredGroups, precision, _pref
	for grp in mgr.pseudoBondGroups:
		if grp.category.startswith('internal-chain'):
			continue
		if not grp.pseudoBonds:
			continue
		if grp in metal_pbg_lookup:
			parent, category = "m%d" % id(metal_pbg_lookup[grp]), "metal coordination bonds"
		else:
			parent, category = None, grp.category
			if category != "missing segments":
				global_groups.append(grp)
		if grp in monitoredGroups:
			monitored = (precision(), _pref["show units"])
		else:
			monitored = False
		data.append([parent, category, monitored, [bond_data(pb) for pb in grp.pseudoBonds]])
		print>>f, "pbg%d = make_pseudobonds(%s)" % (id(grp), pickled(data))
	return global_groups

def write_selection(f):
	from chimera import selection
	sel_atoms = ["a%d" % id(a) for a in selection.currentAtoms()]
	sel_bonds = [("a%d" % id(b.atoms[0]), "a%d" % id(b.atoms[1])) for b in selection.currentBonds()]
	sel_pbs = [("pbg%d" % id(pb.pseudoBondGroup), "a%d" % id(pb.atoms[0]), "a%d" % id(pb.atoms[1]))
		for pb in selection.currentPseudobonds()]
	print>>f, "restore_selection(*%s)" % pickled((sel_atoms, sel_bonds, sel_pbs))

def write_2d_labels(f, model):
	data = []
	global label_2d_multiline_warned
	from chimera import OGLFont
	for label in model.labels:
		if len(label.lines) > 1 and not label_2d_multiline_warned:
			replyobj.warning("ChimeraX export does not support multi-line 2D labels; creating multiple"
				" single-line labels instead")
			label_2d_multiline_warned = True
		for i, line in enumerate(label.lines):
			if not line:
				continue
			rgba, size, style, font_name, text = line_data(line)
			if i > 0:
				from chimera import viewer
				w, h = viewer.windowSize
				pos = (label.pos[0],  label.pos[1] - (i * size) / float(h))
			else:
				pos = label.pos
			data.append([pos, label.background, label.margin, label.outline, label.opacity, rgba, size,
				style & OGLFont.bold, style & OGLFont.italic, font_name, text])
	print>>f, "restore_2d_labels(%s)" % pickled(data)

def write_color_key(f, key):
	from chimera import OGLFont
	keywords = {}
	pos = key.getKeyPosition()
	if pos is None:
		keywords['display'] = False
	else:
		ll, ur = pos
		keywords['pos'] = (min([ll[0], ur[0]]), min([ll[1], ur[1]]))
		keywords['size'] = (abs(ur[0] - ll[0]), abs(ll[1] - ur[1]))
	keywords['font_size'] = key.getFontSize()
	keywords['bold'] = bool(key.getFontStyle() & OGLFont.bold)
	keywords['italic'] = bool(key.getFontStyle() & OGLFont.italic)
	keywords['color_treatment'] = key.getColorTreatment()
	keywords['justification'] = key.getJustification()
	keywords['label_side'] = key.getLabelSide()
	keywords['numeric_label_spacing'] = key.getNumLabelSpacing()
	keywords['label_rgba'] = key.getLabelColor()
	# ChimeraX internally adds 5 to the offset
	keywords['label_offset'] = key.getLabelOffset() - 5
	keywords['font'] = key.getFontTypeface()
	keywords['border'] = key.getBorderColor() != None
	if keywords['border']:
		keywords['border_rgba'] = key.getBorderColor()
	keywords['border_width'] = key.getBorderWidth()
	keywords['ticks'] = key.getTickMarks()
	keywords['tick_length'] = key.getTickLength()
	keywords['tick_thickness'] = key.getTickThickness()
	keywords['rgbas_and_labels'] = key.getRgbasAndLabels()
	print>>f, "restore_color_key(%s)" % pickled(keywords)

def line_data(line):
	global label_2d_rgba_warned, label_2d_size_warned, label_2d_style_warned, label_2d_font_name_warned

	rgbas = set([c.rgba for c in line])
	if len(rgbas) > 1 and not label_2d_rgba_warned:
		replyobj.warning("ChimeraX export does not support mixed colors within a single 2D label line")
		label_2d_rgba_warned = True
	rgba = rgbas.pop()

	sizes = set([c.size for c in line])
	if len(sizes) > 1 and not label_2d_size_warned:
		replyobj.warning("ChimeraX export does not support mixed font sizes within a single 2D label line")
		label_2d_size_warned = True
	size = sizes.pop()

	styles = set([c.style for c in line])
	if len(styles) > 1 and not label_2d_style_warned:
		replyobj.warning("ChimeraX export does not support mixed font styles within a single 2D label line")
		label_2d_style_warned = True
	style = styles.pop()

	font_names = set([c.fontName for c in line])
	if len(font_names) > 1 and not label_2d_font_name_warned:
		replyobj.warning("ChimeraX export does not support mixed fonts within a single 2D label line")
		label_2d_font_name_warned = True
	font_name = font_names.pop()

	return (rgba, size, style, font_name, "".join([unicode(c) for c in line]))

def write_epilog(f, by_id, global_groups, eye_pos, field_of_view, window_size, bg_color):
	for id_info, models in by_id.items():
		model_list = ", ".join(["m%d" % id(m) for m in models])
		if isinstance(id_info, int):
			print>>f, "session.models.add([%s], minimum_id=%d, _from_session=True)" % (
				", ".join(["m%d" % id(m) for m in models]), id_info+1)
		else:
			print>>f, "m%d.add([%s])" % (id(id_info), model_list)
		for model in models:
			if model.__class__.__name__ == 'Molecule':
				print>>f, """
if m%d.num_coordsets > 1:
	from chimerax.core.commands import run
	run(session, "coordset slider %%s" %% m%d.atomspec, log=False)
""" % (id(model), id(model))
	print>>f, "session.models.add([%s], _from_session=True)" % ", ".join(
		["pbg%d" % id(pbg) for pbg in global_groups])
	print>>f, "restore_camera(%s, %s)" % (eye_pos, field_of_view)
	print>>f, "restore_window_size(%d, %d)" % window_size
	print>>f, """
from chimerax.label.label3d import label as label3d
from chimerax.core.objects import Objects
from chimerax.atomic import Atoms
for obj, obj_type, label, color, offset in label_data:
	if obj_type == "atoms":
		objs = Objects(atoms=Atoms([obj]))
	else:
		objs = Objects(atoms=obj.atoms)
	label3d(session, objects=objs, object_type=obj_type, text=label, color=color, offset=offset)
"""
	print>>f, "session.main_view.background_color = (%g, %g, %g, %g)" % (
		(0.0, 0.0, 0.0, 1.0) if bg_color is None else bg_color.rgba())

def write_molecule(f, m):
	tube = False
	from chimera import openModels
	for other_m in openModels.list():
		if other_m.id == m.id and other_m.name.startswith(m.name) and other_m.name.endswith(pipes_string):
			tube = True
			break
	coord_sets = m.coordSets.values()
	coord_sets.sort(lambda cs1, cs2: cmp(cs1.id, cs2.id))
	data = [model_data(m), m.autochain, m.ballScale, color_val(m.color), m.lineType, m.lineWidth,
		m.lowerCaseChains, m.mol2comments, m.mol2data, m.pdbHeaders, m.pdbVersion, m.pointSize,
		m.ribbonHidesMainchain, color_val(m.ribbonInsideColor), m.silhouette, m.stickScale,
		color_val(m.surfaceColor or m.color), m.surfaceOpacity, tube, m.vdwDensity, m.wireStipple,
		[coordset_data(cs) for cs in coord_sets], [res_data(r) for r in m.residues],
		[bond_data(b) for b in m.bonds]]
	print>>f, "m%d = structure_map['m%d'] = make_structure(%s)" % (id(m), id(m), pickled(data))

def write_msms_model(f, m):
	color_mode = m.colorMode
	global surface_export_warned
	if color_mode == 2 and not surface_export_warned:
		replyobj.warning("ChimeraX export does not support per-vertex surface coloring")
		surface_export_warned = True
	data = [model_data(m), m.probeRadius, m.density, color_val(m.color), color_mode,
		["a%d" % id(a) for a in m.atomMap]]
	print>>f, "m%d = make_molecular_surface(%s)" % (id(m), pickled(data))

def model_data(m):
	return [m.name, m.display, m.openState.xform.getOpenGLMatrix()]

def coordset_data(cs):
	return [cs.id, cs.xyzArray()]

def res_data(r):
	try:
		label_coord = r.labelCoord().data()
	except:
		# throws exception if nothing displayed
		label_coord = (0.0,0.0,0.0)
	return ["r%d" % id(r), r.currentLabelOffset().data(), color_val(r.fillColor), r.fillDisplay,
		r.fillMode, r.id.chainId, r.id.insertionCode, r.id.position, r.isHelix, r.isHet, r.isStrand,
		r.label, color_val(r.labelColor or r.ribbonColor or r.molecule.color),
		label_coord, simple_offset(r.labelOffset), color_val(r.ribbonColor or r.molecule.color),
		r.ribbonDisplay, r.ribbonDrawMode, r.ssId, r.type, [atom_data(a) for a in r.atoms]]

def atom_data(a):
	return ["a%d" % id(a), a.altLoc, a.anisoU, a.bfactor if a.haveBfactor else None, a.coord().data(),
		a.coordIndex, a.currentLabelOffset().data(), a.defaultRadius, a.display, a.drawMode,
		a.element.number, a.hide, a.idatmIsExplicit, a.idatmType, a.label,
		color_val(a.labelColor or a.color or a.molecule.color), a.labelCoord().data(),
		simple_offset(a.labelOffset), a.name[:4], a.occupancy if a.haveOccupancy else None, a.radius,
		a.serialNumber, color_val(a.shownColor()), a.surfaceCategory, color_val(a.surfaceColor or
		a.molecule.color), a.surfaceDisplay, a.surfaceOpacity, a.vdw, color_val(a.vdwColor)]

def bond_data(b):
	"""Also used for pseudobonds"""
	return [["a%d" % id(a) for a in b.atoms], [color_val(b.color), b.currentLabelOffset().data(), b.display,
		b.drawMode, b.halfbond, b.label, color_val(b.labelColor), b.labelCoord().data(),
		simple_offset(b.labelOffset), b.radius]]

def write_alignment(f, mav):
	assoc_info = {}
	for mol, seq in mav.associations.items():
		match_map = seq.matchMaps[mol]
		info = {}
		for key, value in match_map.items():
			if type(key) == int:
				info[key] = "r%d" % id(value)
		assoc_info[mav.seqs.index(seq)] = info

	from MultAlignViewer.RegionBrowser import SEL_REGION_NAME
	rb = mav.regionBrowser
	regions = [(r, None) for r in rb.regions if r.name != SEL_REGION_NAME]
	for seq, seq_regions in rb.sequenceRegions.items():
		regions.extend([(r, seq) for r in seq_regions])
	region_data = []
	for r, seq in regions:
		region_data.append([r.name, region_blocks(mav, r), r.shown, r.interiorRGBA, r.borderRGBA,
			r.highlighted, r.coverGaps, None if seq is None else mav.seqs.index(seq)])
	data = [mav.title, [seq_data(seq) for seq in mav.seqs], hdr_data(mav), assoc_info, region_data]
	print>>f, "make_alignment(%s)" % pickled(data)

def region_blocks(mav, region):
	blocks = []
	for block in region.blocks:
		line1, line2, i1, i2 = block
		if line1 not in mav.seqs and line2 not in mav.seqs:
			continue
		seq1 = line1 if line1 in mav.seqs else mav.seqs[0]
		seq2 = line2 if line2 in mav.seqs else mav.seqs[0]
		blocks.append((mav.seqs.index(seq1), mav.seqs.index(seq2), i1, i2))
	return blocks

def seq_data(seq):
	return [seq.name, str(seq)]

def hdr_data(mav):
	data = { hdr.name: [] for hdr in mav.headers(shownOnly=True) }
	return data

restore_model_definition = """
def restore_model_data(m, data):
	name, display, matrix = data
	m.name = name
	m.display = display
	from chimerax.geometry import Place
	from numpy import array, float64, transpose
	reformatted = []
	for r in range(3):
		vector = []
		reformatted.append(vector)
		for c in range(4):
			vector.append(matrix[r + 4*c])
	m.scene_position = Place(array(reformatted, float64))
"""

make_structure_definition = """
def make_structure(structure_data):
	model_data, auto_chain, ball_scale, color, line_type, line_width, lower_case_chains, mol2_comments, \\
		mol2_data, pdb_headers, pdb_version, point_size, ribbon_hides_backbone, ribbon_inside_color, \\
		silhouette, stick_scale, surface_color, surface_opacity, tube, vdw_density, wire_stipple, \\
		coordset_data, res_data, bond_data = structure_data
	from chimerax.atomic import AtomicStructure
	s = AtomicStructure(session, auto_style=False)
	for key, value in pdb_headers.items():
		s.set_metadata_entry(key, value)
	restore_model_data(s, model_data)
	cs_id_lookup = restore_coordset_data(s, coordset_data)
	if tube:
		# s.ribbon_mode_strand is actively buggy
		s.ribbon_mode_helix = s.RIBBON_MODE_ARC
	atom_map = restore_res_data(s, res_data, ribbon_hides_backbone, tube)
	restore_bond_data(s, bond_data, atom_map, stick_scale)
	global_atom_map.update(atom_map)
	s.pdb_version = pdb_version
	s._set_chain_descriptions(session)
	from chimerax.pdb.pdb import set_logging_info
	set_logging_info(s)
	return s
"""

make_molecular_surface_definition = """
def make_molecular_surface(surface_data):
	model_data, probe_radius, density, color, color_mode, atoms = surface_data

	from chimerax.atomic import Atoms
	enclose_atoms = Atoms([global_atom_map[a_id] for a_id in atoms])
	show_atoms = Atoms([a for a in enclose_atoms if a._c2cx_surface_display])
	from math import sqrt
	grid = density / sqrt(3.0)
	from chimerax.atomic.molsurf import MolecularSurface
	s = MolecularSurface(session, enclose_atoms, show_atoms, probe_radius, grid, None, None,
		"name restored later", rgba_to_color(color).uint8x4(), None, True)
	restore_model_data(s, model_data)
	if s.name.startswith("MSMS "):
		s.name = s.name[5:]
	s.calculate_surface_geometry()
	s.auto_update = True
	if color_mode == 1:
		# per atom
		from numpy import array
		s.color_atom_patches(enclose_atoms, None, array([a._c2cx_surface_color for a in enclose_atoms]))
	elif color_mode == 2:
		session.logger.warning("Cannot restore custom surface colors")
	return s
"""

restore_coordset_definition = """
def restore_coordset_data(s, coordset_data):
	xyzs = []
	id_lookup = {}
	for i, data in enumerate(coordset_data):
		cs_id, cs_xyzs = data
		id_lookup[cs_id] = i
		xyzs.append(cs_xyzs)
	from numpy import array
	s.add_coordsets(array(xyzs))
	return id_lookup
"""

restore_res_definition = """
def restore_res_data(s, res_data, ribbon_hides_backbone, tube):
	global residue_map
	atom_map = {}
	for restore_id, current_label_offset, fill_color, fill_display, fill_mode, chain_id, \\
			insertion_code, position, is_helix, is_het, is_strand, label, label_color, \\
			label_coord, label_offset, ribbon_color, ribbon_display, ribbon_draw_mode, \\
			ss_id, type, atom_data in res_data:
		r = s.new_residue(type, chain_id, position, insertion_code)
		residue_map[restore_id] = r
		r.ss_type = r.SS_HELIX if is_helix else (r.SS_STRAND if is_strand else r.SS_COIL)
		r.ss_id = ss_id
		r.ribbon_display = ribbon_display or tube
		r.ribbon_color = rgba_to_color(ribbon_color).uint8x4()
		r.ribbon_hide_backbone = ribbon_hides_backbone
		if label:
			global label_data
			label_data.append((r, "residues", label, rgba_to_color(label_color), label_offset))
		atom_map.update(restore_atom_data(r, atom_data))
	xsection = s.ribbon_xs_mgr.STYLE_ROUND if ribbon_draw_mode == 2 else s.ribbon_xs_mgr.STYLE_SQUARE
	s.ribbon_xs_mgr.set_coil_style(xsection)
	s.ribbon_xs_mgr.set_helix_style(xsection)
	s.ribbon_xs_mgr.set_sheet_style(xsection)
	return atom_map
"""

restore_atom_definition = """
def restore_atom_data(r, atom_data):
	atom_map = {}
	draw_mode_map = { 0: 2, 1: 0, 2: 2, 3: 1 }
	for restore_id, alt_loc, aniso_u, bfactor, coord, coord_index, current_label_offset, \\
			default_radius, display, draw_mode, element_number, hide, idatm_is_explicit, idatm_type, label, \\
			label_color, label_coord, label_offset, name, occupancy, radius, serial_number, \\
			shown_color, surface_category, surface_color, surface_display, surface_opacity, vdw, vdw_color \\
			in atom_data:
		a = r.find_atom(name)
		if alt_loc != '' and a:
			a.set_alt_loc(alt_loc, True)
			a.coord = coord
		else:
			a = r.structure.new_atom(name, element_number)
			r.add_atom(a)
			if alt_loc != '':
				a.set_alt_loc(alt_loc, True)
			a.coord_index = coord_index
			a.display = display
			a.draw_mode = draw_mode_map[draw_mode]
			a.hide = a.HIDE_RIBBON if hide else 0
			if idatm_is_explicit:
				a.idatm_type = idatm_type
			if label:
				global label_data
				label_data.append((a, "atoms", label, rgba_to_color(label_color), label_offset))
			if default_radius != radius:
				a.radius = radius
			a.color = rgba_to_color(shown_color).uint8x4()
		if occupancy is not None:
			a.occupancy = occupancy
		if bfactor is not None:
			a.bfactor = bfactor
		a.serial_number = serial_number
		a.aniso_u6 = tuple((aniso_u[x][y]
			for x,y in [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)])) if aniso_u else None
		a._c2cx_surface_display = surface_display
		scolor = rgba_to_color(surface_color)
		if surface_opacity >= 0.0:
			scolor.rgba[-1] = surface_opacity
		a._c2cx_surface_color = scolor.uint8x4()
		atom_map[restore_id] = a
	return atom_map
"""

restore_bond_definition = """
def restore_bond_data(s, bond_data, atom_map, stick_scale):
	for atom_ids, bond_attrs in bond_data:
		# bonds between alt locs in Chimera1 can be duplicative in ChimeraX...
		a1, a2 = [atom_map[ident] for ident in atom_ids]
		if a1 not in a2.neighbors:
			b = s.new_bond(a1, a2)
			if stick_scale != 1.0:
				b.radius *= stick_scale
"""

make_pseudobonds_definition = """
def make_pseudobonds(pseudobond_data):
	for parent, category, monitored, pb_data in pseudobond_data:
		if category == "missing segments":
			for atom_ids, pb_attrs in pb_data:
				a1, a2 = [global_atom_map[ident] for ident in atom_ids]
				pbg = a1.structure.pseudobond_group("missing structure")
				pb = pbg.new_pseudobond(a1, a2)
				if a1.connects_to(a2):
					a1.structure.delete_bond(a1.bonds[a1.neighbors.index(a2)])
			continue
		if parent is None:
			pbg = session.pb_manager.get_group(category)
		else:
			pbg = structure_map[parent].pseudobond_group(category)
		if monitored:
			precision, show_units = monitored
			monitor = pbg.session.pb_dist_monitor
			monitor.add_group(pbg, session_restore=True)
			monitor.decimal_places = precision
			monitor.show_units = show_units
		for atom_ids, pb_attrs in pb_data:
			a1, a2 = [global_atom_map[ident] for ident in atom_ids]
			pb = pbg.new_pseudobond(a1, a2)
	return pbg
"""

color_converter_definition = """
def rgba_to_color(rgba):
	from chimerax.core.colors import Color
	return Color(rgba)
"""

restore_selection_definition = """
def restore_selection(atom_data, bond_data, pb_data):
	session.selection.clear()
	for atom_id in atom_data:
		global_atom_map[atom_id].selected = True
	for id1, id2 in bond_data:
		a1 = global_atom_map[id1]
		a2 = global_atom_map[id2]
		for b in a1.bonds:
			if b.other_atom(a1) == a2:
				b.selected = True
				break
	for pbg_id, atom_id1, atom_id2 in pb_data:
		pbg = eval(pbg_id)
		pb_atoms = set([global_atom_map[aid] for aid in (atom_id1, atom_id2)])
		for pb in pbg.pseudobonds:
			if set(pb.atoms) == pb_atoms:
				pb.selected = True
"""

restore_camera_definition = """
def restore_camera(eye_pos, field_of_view):
	from chimerax.geometry import translation
	session.view.camera.position = translation(eye_pos)
	session.view.camera.field_of_view = field_of_view
	# as a nicety, multiply all matrices by the inverse of the
	# lowest-ID model's matrix, so that newly opened models
	# are in correct relative position to that model at least
	models = session.models[:]
	if models:
		models.sort()
		inverse = models[0].position.inverse(is_orthonormal=True)
		for m in models:
			m.position = m.position * inverse
		session.view.camera.position = inverse * session.view.camera.position
"""

restore_window_size_definition = """
def restore_window_size(w,h):
	from chimerax.graphics.windowsize import window_size
	window_size(session, w, h)
"""

make_alignment_definition = """
def make_alignment(alignment_data):
	title, seq_data, hdr_data, assoc_info, region_data = alignment_data
	seqs = [restore_seq(sd) for sd in seq_data]
	aln = session.alignments.new_alignment(seqs, title, auto_associate=False)
	for hdr in aln.headers:
		hdr.shown = hdr.name in hdr_data
	global residue_map
	from chimerax.atomic import SeqMatchMap
	for seq_index, map_info in assoc_info.items():
		if not map_info:
			continue
		align_seq = aln.seqs[seq_index]
		match_map = SeqMatchMap(align_seq, residue_map[list(map_info.values())[0]].chain)
		for index, res_id in map_info.items():
			match_map.match(residue_map[res_id], index)
		aln.prematched_assoc_structure(match_map, False, False)
	if aln.viewers:
		viewer = aln.viewers[0]
		for name, blocks, shown, fill, outline, highlighted, cover_gaps, seq_index in region_data:
			viewer.new_region(name=name, blocks=blocks, fill=fill, outline=outline, select=highlighted,
				cover_gaps=cover_gaps, session_restore=True, shown=shown,
				sequence=(None if seq_index is None else aln.seqs[seq_index]))
"""

restore_seq_definition = """
def restore_seq(seq_data):
	name, sequence = seq_data
	from chimerax.atomic import Sequence
	return Sequence(name=name, characters=sequence)
"""

restore_2d_labels_definition = """
def restore_2d_labels(labels_data):
	from chimerax.label.label2d import Label
	from chimerax.core.colors import Color
	used = set()
	for pos, background, margin, outline, opacity, rgba, size, bold, italic, font, text in labels_data:
		if background is not None:
			background = Color(background)
		font = font_mapping.get(font, font)
		label_name = text
		label_num = 1
		while label_name in used:
			label_name = "%s (%d)" % (label_name, label_num)
			label_num += 1
		Label(session, label_name, text=text, xpos=pos[0], ypos=pos[1],
			color=Color(rgba).uint8x4(), size=size, font=font, bold=bold, italic=italic,
			background=background, outline_width=outline, margin=margin)
"""

restore_color_key_definition = """
def restore_color_key(key_info):
	from chimerax.color_key.model import get_model
	key = get_model(session)
	for k, v in key_info.items():
		if k == "font":
			v = font_mapping.get(v, v)
		setattr(key, k, v)
"""
